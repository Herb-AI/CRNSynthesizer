struct GenericAtomConnections <: AbstractLocalConstraint
    path::Vector{Int}
    grammar_data::GrammarData
end

import Base.==
function ==(a::GenericAtomConnections, b::GenericAtomConnections)
    return a.path == b.path
end

import Base.hash
function hash(a::GenericAtomConnections, h::UInt)
    return hash(a.path, h)
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::GenericAtomConnections, path::Vector{Int}
)::Bool
    # Check if the path is a child of the atom node
    if length(path) <= length(constraint.path)
        return false
    end

    if path[1:length(constraint.path)] != constraint.path
        return false
    end

    # Schedule for atom nodes and relevant holes
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    # Only schedule for types that can change bonds or atoms
    return type == :atom ||
           type == :bond ||
           (
               node isa Hole && (
               type == :chain ||
               type == :structure ||
               type == :branch ||
               type == :branches ||
               type == :molecule_list ||
               type == :molecule ||
               type == :fragment_X_entry ||
               type == :fragment_X_exit
           )
           )
end

function HerbConstraints.propagate!(solver::Solver, constraint::GenericAtomConnections)
    propagate_atoms!(solver, constraint, constraint.path)
end

function get_relevant_bonds(solver::GenericSolver, path::Vector{Int})
    # Get the relevant bonds for the current node
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    @match type begin
        :ringbond => begin
            return [push!(copy(path), 1)], []
        end
        :ringbonds => begin
            if !isfilled(node)
                return Vector{Vector{Int}}(), [path]
            end

            rule = HerbCore.get_rule(node)
            @match solver.grammar.rules[rule] begin
                :(ringbond *
                ringbonds) => begin
                    ringbond_bonds,
                    ringbond_holes = get_relevant_bonds(
                        solver, push!(copy(path), 1)
                    )
                    ringbonds_bonds,
                    ringbonds_holes = get_relevant_bonds(
                        solver, push!(copy(path), 2)
                    )
                    return vcat(ringbond_bonds, ringbonds_bonds),
                    vcat(ringbond_holes, ringbonds_holes)
                end
                :("") => return Vector{Vector{Int}}(), Vector{Vector{Int}}()
                x => throw("Unknown ringbond rule: $x")
            end
        end

        :branch => begin
            return [push!(copy(path), 1)], []
        end
        :branches => begin
            if !isfilled(node)
                return Vector{Vector{Int}}(), [path]
            end

            rule = HerbCore.get_rule(node)
            @match solver.grammar.rules[rule] begin
                :(branch *
                branches) => begin
                    branch_bonds,
                    branch_holes = get_relevant_bonds(
                        solver, push!(copy(path), 1)
                    )
                    branches_bonds,
                    branches_holes = get_relevant_bonds(
                        solver, push!(copy(path), 2)
                    )
                    return vcat(branch_bonds, branches_bonds),
                    vcat(branch_holes, branches_holes)
                end
                :("") => return Vector{Vector{Int}}(), Vector{Vector{Int}}()
                x => throw("Unknown branch rule: $x")
            end
        end

        x => throw("Unknown node type: $x")
    end

    return Vector{Vector{Int}}(), Vector{Vector{Int}}()
end

# Function that searches for atoms in the tree and propagates the atom constraint on them
function propagate_atoms!(
        solver::Solver,
        constraint::GenericAtomConnections,
        path::Vector{Int};
        bond_paths::Vector{Vector{Int}} = Vector{Vector{Int}}(),
        holes::Vector{Vector{Int}} = Vector{Vector{Int}}()
)
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    # Check if the node has children we can propagate to
    if isuniform(node)
        if !isfilled(node)
            return
        end
        rule = HerbCore.get_rule(node)
        if type == :fragment_X_entry
            for i in eachindex(node.children)
                propagate_atoms!(solver, constraint, push!(copy(path), i))
            end
            return
        elseif type == :fragment_X_exit
            if rule == :("(-" * chain * ")")
                propagate_atoms!(
                    solver, constraint, push!(copy(path), 1), bond_paths = [copy(path)])
            elseif rule == :("(-" * fragment_X_entry * ")")
                propagate_atoms!(solver, constraint, push!(copy(path), 1))
            end
            return
        end

        @match solver.grammar.rules[rule] begin
            # Main molecule rule with one child
            :chain => propagate_atoms!(solver, constraint, push!(copy(path), 1))
            :(from_SMILES(chain)) => propagate_atoms!(
                solver, constraint, push!(copy(path), 1))

            # Chain option with two children
            :(atom *
            ringbonds) => begin
                ring_bonds, ring_holes = get_relevant_bonds(solver, push!(copy(path), 2))
                bonds::Vector{Vector{Int}} = vcat(ring_bonds, bond_paths)
                propagate_atom!(solver, constraint, push!(copy(path), 1), bonds, holes)
            end

            # Chain option with three children
            :(SMILES_combine_chain(bond,
            structure,
            chain)) => begin
                # The structule will have the bond and previous collected bond_paths and holes
                bond_paths = push!(copy(bond_paths), push!(copy(path), 1))
                propagate_atoms!(
                    solver,
                    constraint,
                    push!(copy(path), 2),
                    bond_paths = bond_paths,
                    holes = holes
                )

                # The chain will only have the bond be relevant
                propagate_atoms!(
                    solver,
                    constraint,
                    push!(copy(path), 3),
                    bond_paths = [push!(copy(path), 1)]
                )
            end

            # Chain option with three children
            :(structure * bond *
            chain) => begin
                # The structure will have the bond and previous collected bond_paths and holes
                bond_paths = push!(copy(bond_paths), push!(copy(path), 2))
                propagate_atoms!(
                    solver,
                    constraint,
                    push!(copy(path), 1),
                    bond_paths = bond_paths,
                    holes = holes
                )

                # The chain will only have the bond be relevant
                propagate_atoms!(
                    solver,
                    constraint,
                    push!(copy(path), 3),
                    bond_paths = [push!(copy(path), 2)]
                )
            end

            :(structure * "-" * fragment_X_entry) => begin
                bond_paths = push!(copy(bond_paths), push!(copy(path), 2))
                propagate_atoms!(
                    solver, constraint, push!(copy(path), 1),
                    bond_paths = bond_paths, holes = holes
                )
                propagate_atoms!(solver, constraint, push!(copy(path), 2))
            end

            # Structure option with three children
            :(atom * ringbonds *
            branches) => begin
                # TODO: Check if the following propagate_atoms! is correct
                # propagate_atoms!(solver, constraint, push!(copy(path), 3), bond_paths=bond_paths, holes=holes)

                ring_bonds, ring_holes = get_relevant_bonds(solver, push!(copy(path), 2))
                branch_bonds,
                branch_holes = get_relevant_bonds(
                    solver, push!(copy(path), 3)
                )
                bonds::Vector{Vector{Int}} = vcat(ring_bonds, branch_bonds, bond_paths)
                holes::Vector{Vector{Int}} = vcat(ring_holes, branch_holes, holes)
                propagate_atom!(solver, constraint, push!(copy(path), 1), bonds, holes)
            end

            x => throw("Unknown rule type: $x")
        end
    end
end

function propagate_atom!(
        solver::Solver,
        constraint::GenericAtomConnections,
        atom_path::Vector{Int},
        bond_paths::Vector{Vector{Int}},
        holes::Vector{Vector{Int}}
)
    # Get the possible atoms the current node can be
    atom = get_node_at_location(solver, atom_path)
    rules = get_rules(atom)

    # println("tree: ", solver.state.tree)

    # Get the max and min number of connections for the atom
    max_connections = 0
    min_connections = 1000 # TODO: Look into max int value
    for rule in rules
        connections = grammar_to_atom_connections(constraint.grammar_data, rule)
        max_connections = max(max_connections, connections)
        min_connections = min(min_connections, connections)
    end

    length_bonds = length(bond_paths)
    length_holes = length(holes)

    # println()
    # println("min: ", min_connections, "  max: ", max_connections)
    # println("bonds: ", length_bonds, "  holes: ", length_holes)
    # println("bonds: ", bond_paths)
    # println("holes: ", holes)

    if length_holes == 0 && length_bonds < min_connections
        HerbConstraints.set_infeasible!(solver)
        # println("infeasible")
        return nothing
    end

    if length_bonds + length_holes == 0
        HerbConstraints.set_infeasible!(solver)
        # println("infeasible")
        return nothing
    end

    if length_bonds > max_connections
        HerbConstraints.set_infeasible!(solver)
        return nothing
    end

    if length_bonds == max_connections
        for hole in holes
            terminals = findall(==(:("")), solver.grammar.rules)
            c_remove_all_but!(solver, hole, terminals, false)
        end
    end
end
