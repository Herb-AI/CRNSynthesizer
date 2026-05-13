mutable struct LocalUniformAtomConnections <: AbstractLocalConstraint
    path::Vector{Int}
    bond_paths::Vector{Vector{Int}}
    grammar_data::GrammarData
end

import Base.==
function ==(a::LocalUniformAtomConnections, b::LocalUniformAtomConnections)
    return a.path == b.path
end

import Base.hash
function hash(a::LocalUniformAtomConnections, h::UInt)
    return hash(a.path, h) + hash(a.bond_paths, h)
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::LocalUniformAtomConnections, path::Vector{Int}
)::Bool
    if path == constraint.path || path in constraint.bond_paths
        return true
    end
    return false
end

function HerbConstraints.propagate!(solver::Solver, constraint::LocalUniformAtomConnections)
    # Get the possible atoms the current node can be
    atom = get_node_at_location(solver, constraint.path)
    rules = get_rules(atom)

    # Get the max and min number of connections for the atom
    max_connections = 0
    min_connections = 1000 # TODO: Look into max int value
    for rule in rules
        connections = grammar_to_atom_connections(constraint.grammar_data, rule)
        max_connections = max(max_connections, connections)
        min_connections = min(min_connections, connections)
    end

    # Get the max and min number of connections from the bonds
    max_bonds = 0
    min_bonds = 0
    for bond_path in constraint.bond_paths
        bond = get_node_at_location(solver, bond_path)
        bond_rules = get_rules(bond)
        # println("bond rules: ", bond_rules)

        local_max = 0
        local_min = 1000 # TODO: Look into max int value
        for bond_rule in bond_rules
            connections = grammar_to_bond_connections(constraint.grammar_data, bond_rule)
            local_max = max(local_max, connections)
            local_min = min(local_min, connections)
        end
        max_bonds += local_max
        min_bonds += local_min
    end

    # If the ranges don't match, set the constraint to infeasible
    if min_bonds > max_connections || min_connections > max_bonds
        HerbConstraints.set_infeasible!(solver)
        return nothing
    end

    # println()
    # println("max_connections: ", max_connections, " min_connections: ", min_connections)
    # println("max_bonds: ", max_bonds, " min_bonds: ", min_bonds)
    # println(length(constraint.bond_paths), " bonds: ", constraint.bond_paths)

    # TODO: find a better way to do this
    # Update all the ranges
    available_atoms = Int[]
    for rule in rules
        connections = grammar_to_atom_connections(constraint.grammar_data, rule)
        if connections >= min_bonds && connections <= max_bonds
            push!(available_atoms, rule)
        end
    end
    c_remove_all_but!(solver, constraint.path, available_atoms, false)
end

function get_relevant_bonds(solver, path)
    grammar = solver.grammar
    node = get_node_at_location(solver, path)

    @match get_node_type(grammar, node) begin
        :ringbonds => begin
            rule = grammar.rules[HerbCore.get_rule(node)]

            if rule == :("")
                return Vector{Vector{Int}}()
            elseif rule == :(ringbond * ringbonds)
                ringbond = get_relevant_bonds(solver, push!(copy(path), 1))
                ringbonds = get_relevant_bonds(solver, push!(copy(path), 2))
                return vcat(ringbond, ringbonds)
            else
                throw("Unknown ringbond rule: $rule")
            end
        end

        :ringbond => begin
            bond_path = push!(copy(path), 1)
            return [bond_path]
        end

        :branches => begin
            rule = grammar.rules[HerbCore.get_rule(node)]

            if rule == :("")
                return []
            elseif rule == :(branch * branches)
                branch = get_relevant_bonds(solver, push!(copy(path), 1))
                branches = get_relevant_bonds(solver, push!(copy(path), 2))
                return vcat(branch, branches)
            else
                throw("Unknown branch rule: $rule")
            end
        end

        :branch => begin
            bond_path = push!(copy(path), 1)
            return [bond_path]
        end

        unknown => begin
            throw("Unknown node type: $unknown")
        end
    end
end

function post_atom_constraints!(
        solver::Solver,
        path::Vector{Int},
        grammar_data::GrammarData;
        relevant_bonds::Vector{Vector{Int}} = Vector{Vector{Int}}()
)
    grammar = solver.grammar
    node = get_node_at_location(solver, path)

    @match get_node_type(grammar, node) begin
        :molecule => begin
            post_atom_constraints!(solver, push!(copy(path), 1), grammar_data)
        end

        :chain => begin
            rule = grammar.rules[HerbCore.get_rule(node)]

            if rule == :(atom * ringbonds)
                bonds = get_relevant_bonds(solver, push!(copy(path), 2))
                bonds = vcat(bonds, relevant_bonds)

                atom_constraint = LocalUniformAtomConnections(
                    push!(copy(path), 1), bonds, grammar_data
                )
                HerbConstraints.post!(solver, atom_constraint)
            elseif rule == :(SMILES_combine_chain(bond, structure, chain))
                bond_path = push!(copy(path), 1)
                bonds = [relevant_bonds..., bond_path]

                post_atom_constraints!(
                    solver, push!(copy(path), 2), grammar_data, relevant_bonds = bonds
                )
                post_atom_constraints!(
                    solver,
                    push!(copy(path), 3),
                    grammar_data,
                    relevant_bonds = [bond_path]
                )
            elseif rule == :(structure * bond * chain)
                bond_path = push!(copy(path), 2)
                bonds = [relevant_bonds..., bond_path]

                post_atom_constraints!(
                    solver, push!(copy(path), 1), grammar_data, relevant_bonds = bonds
                )
                post_atom_constraints!(
                    solver,
                    push!(copy(path), 3),
                    grammar_data,
                    relevant_bonds = [bond_path]
                )
            elseif rule == :(structure * "-" * fragment_X_entry)
                bond_path = push!(copy(path), 2)
                bonds = [relevant_bonds..., bond_path]

                post_atom_constraints!(
                    solver, push!(copy(path), 1), grammar_data, relevant_bonds = bonds
                )
                post_atom_constraints!(solver, push!(copy(path), 2), grammar_data)
            else
                throw("Unknown chain rule: $rule")
            end
        end

        :structure => begin
            ringbonds = get_relevant_bonds(solver, push!(copy(path), 2))
            branches = get_relevant_bonds(solver, push!(copy(path), 3))
            bonds = vcat(ringbonds, branches, relevant_bonds)

            atom_constraint = LocalUniformAtomConnections(
                push!(copy(path), 1), bonds, grammar_data
            )
            HerbConstraints.post!(solver, atom_constraint)

            post_atom_constraints!(solver, push!(copy(path), 3), grammar_data)
        end

        :fragment_X_exit => begin
            rule = grammar.rules[HerbCore.get_rule(node)]

            if rule == :("(-" * chain * ")")
                post_atom_constraints!(solver, push!(copy(path), 1), grammar_data,
                    relevant_bonds = [copy(path)])
            elseif rule == :("(" * special_bond * fragment_X_entry * ")")
                post_atom_constraints!(solver, push!(copy(path), 2), grammar_data)
            end
        end

        :fragment_X_entry || :starting_fragment => begin
            for i in eachindex(node.children)
                post_atom_constraints!(solver, push!(copy(path), i), grammar_data)
            end
        end

        :branches => begin
            rule = grammar.rules[HerbCore.get_rule(node)]

            if rule == :("")
                return nothing
            elseif rule == :(branch * branches)
                post_atom_constraints!(solver, push!(copy(path), 1), grammar_data)
                post_atom_constraints!(solver, push!(copy(path), 2), grammar_data)
            else
                throw("Unknown branch rule: $rule")
            end
        end

        :branch => begin
            bond_path = push!(copy(path), 1)
            post_atom_constraints!(
                solver, push!(copy(path), 2), grammar_data, relevant_bonds = [bond_path]
            )
        end

        unknown => begin
            println(unknown)
            throw("Unknown node type: $unknown")
        end
    end
end
