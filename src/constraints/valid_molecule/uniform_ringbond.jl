mutable struct localUniformRingbonds <: AbstractLocalConstraint
    path::Vector{Int}
    ringbond_paths::Vector{Tuple{Vector{Int}, Vector{Int}}}
    connected_groups::Vector{Vector{Vector{Int}}}
    grammar_data::GrammarData
end

import Base.==
function ==(a::localUniformRingbonds, b::localUniformRingbonds)
    return a.path == b.path && a.ringbonds == b.ringbonds
end

import Base.hash
function hash(a::localUniformRingbonds, h::UInt)
    return hash(a.path, h) + hash(a.ringbond_paths, h)
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::localUniformRingbonds, path::Vector{Int}
)::Bool
    # Check if the update was in a ringbond
    for ringbond_path in [x[1] for x in constraint.ringbond_paths]
        if path[1:(end - 1)] == ringbond_path
            return true
        end
    end

    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    # If the update was in a fragment that can contain ringbonds,
    # we need to schedule and recompute the incompatible groups
    if type == :fragment_X_entry || type == :starting_fragment
        constraint.ringbond_paths, constraint.connected_groups, _ = get_ringbond_paths(
            solver, constraint.path)
        return true
    end

    return false
end

function get_digit_path(solver::Solver, ringbond_path::Vector{Int})::Vector{Int}
    node = get_node_at_location(solver, ringbond_path)
    type = get_node_type(solver.grammar, node)

    if type == :ringbond
        return push!(copy(ringbond_path), 2)
    end
    return push!(copy(ringbond_path), 1)
end

function get_ringbond_bond_info(
        solver::Solver, grammar_data::GrammarData, ringbond_path::Vector{Int})
    ringbond_node = get_node_at_location(solver, ringbond_path)
    ringbond_type = get_node_type(solver.grammar, ringbond_node)

    if ringbond_type == :ringbond
        bond_path = push!(copy(ringbond_path), 1)
        bond_node = get_node_at_location(solver, bond_path)
        bond_rules = get_rules(bond_node)
        return bond_path, bond_rules
    end
    rule = solver.grammar.rules[HerbCore.get_rule(ringbond_node)]
    bond_rules = [bond_to_grammar(grammar_data, rule.args[2])]
    return ringbond_path, bond_rules
end

function HerbConstraints.propagate!(solver::Solver, constraint::localUniformRingbonds)

    # TODO: make dynamic; TEMP: restrict max ringbonds digit to 9
    max_ringbond_grammar = 2
    # Constraint the max digit of the ringbonds
    max_ringbonds = length(constraint.ringbond_paths) ÷ 2
    for (i, (ringbond_path, _)) in enumerate(constraint.ringbond_paths)
        rule = digit_to_grammar(
            constraint.grammar_data, min(max_ringbonds, max_ringbond_grammar)
        )
        grammar_i = digit_to_grammar(constraint.grammar_data, min(i, max_ringbond_grammar))
        node = get_node_at_location(solver, get_digit_path(solver, ringbond_path))
        if node isa StateHole
            c_remove_above!(
                solver, get_digit_path(solver, ringbond_path),
                min(rule, grammar_i); fix_point = false
            )
        end
    end

    # Get the filled ringbond digits
    filled_ringbonds = Dict{Int, Tuple}()
    for (ringbond_path, _) in constraint.ringbond_paths
        digit = get_node_at_location(solver, get_digit_path(solver, ringbond_path))
        if isfilled(digit)
            rule = HerbCore.get_rule(digit)
            if haskey(filled_ringbonds, rule)
                if !isnothing(filled_ringbonds[rule][2])
                    HerbConstraints.set_infeasible!(solver)
                end

                (first_rb, _) = filled_ringbonds[rule]
                filled_ringbonds[rule] = (first_rb, ringbond_path)
            else
                filled_ringbonds[rule] = (ringbond_path, nothing)
            end
        end
    end

    # If the ringbond digit is the same, they both should have the same bond
    for (_, pair) in collect(filled_ringbonds)
        if !isnothing(pair[2])
            # Get the bond of the first ringbond
            left_path, left_rules = get_ringbond_bond_info(
                solver, constraint.grammar_data, pair[1])

            # Get the bond of the second ringbond
            right_path, right_rules = get_ringbond_bond_info(
                solver, constraint.grammar_data, pair[2])

            # Make their domains equal (only if the path wasn't mocked)
            if left_path != pair[1]
                c_remove_all_but!(solver, left_path, right_rules, false)
            end
            if right_path != pair[2]
                c_remove_all_but!(solver, right_path, left_rules, false)
            end
        end
    end

    # Check that a group should not have the same digits
    for group in constraint.connected_groups
        seen = Set{Int}()
        # println("Group: ", group)
        # Filter empty paths from the group
        group = filter(path -> !isempty(path), group)
        for ringbond_path in group
            # println("Ringbond path: ", ringbond_path)
            digit = get_node_at_location(solver, get_digit_path(solver, ringbond_path))
            # println("Digit: ", digit)
            if isfilled(digit)
                rule = HerbCore.get_rule(digit)
                if rule in seen
                    HerbConstraints.set_infeasible!(solver)
                    # println("Infeasible: ", rule)
                end
                push!(seen, rule)
            end
        end
        # println("Seen: ", seen)
    end

    # Check that there are no two ringbonds between the same two atoms
    filled_atoms = Dict{Int, Tuple}()
    for (ringbond_path, atom_path) in constraint.ringbond_paths
        digit = get_node_at_location(solver, get_digit_path(solver, ringbond_path))
        if isfilled(digit)
            rule = HerbCore.get_rule(digit)
            if haskey(filled_atoms, rule)
                if !isnothing(filled_atoms[rule][2])
                    HerbConstraints.set_infeasible!(solver)
                end

                (first_rb, _) = filled_atoms[rule]
                filled_atoms[rule] = (first_rb, atom_path)
            else
                filled_atoms[rule] = (atom_path, nothing)
            end
        end
    end
    # println("Filled atoms: ", filled_atoms)
    filled_atom_pairs = Dict{Tuple, Int}()
    for (rule, (atom1, atom2)) in collect(filled_atoms)
        if !isnothing(atom2)
            highest_atom = max(atom1, atom2)
            lowest_atom = min(atom1, atom2)
            if haskey(filled_atom_pairs, (highest_atom, lowest_atom))
                HerbConstraints.set_infeasible!(solver)
                # println("Infeasible: ", rule, " ", (highest_atom, lowest_atom), " ", filled_atom_pairs[(highest_atom, lowest_atom)])
            else
                filled_atom_pairs[(highest_atom, lowest_atom)] = rule
            end
        end
    end
    # println("Filled pairs: ", filled_atom_pairs)
end

function is_single_atom(smiles::String)::Bool
    return count(==(']'), smiles) == 1
end

function postprocess_grouped_fragment_connection_points(
        solver::UniformSolver, path::Vector{Int},
        forbidden_group::Vector{Vector{Int}}, connects_to_single_atom::Bool,
        all_ringbonds::Vector{Tuple{Vector{Int}, Vector{Int}}},
        all_forbidden::Vector{Vector{Vector{Int}}},
        atom_ringbonds::Vector{Tuple{Vector{Int}, Vector{Int}}},
        branch_children::Vector{Int})::Vector{Vector{Int}}
    forbidden_for_branches = [x[1] for x in atom_ringbonds]
    atom_forbidden = connects_to_single_atom ?
                     vcat(forbidden_group, forbidden_for_branches) :
                     forbidden_for_branches
    push!(all_forbidden, atom_forbidden)
    append!(all_ringbonds, atom_ringbonds)
    for child_id in branch_children
        branch_ringbonds, branch_forbidden,
        _ = get_ringbond_paths(
            solver,
            push!(copy(path), child_id),
            forbidden_group = forbidden_for_branches
        )
        append!(all_ringbonds, branch_ringbonds)
        append!(all_forbidden, branch_forbidden)
    end
    return forbidden_for_branches
end

function get_ringbond_paths(
        solver::UniformSolver,
        path::Vector{Int};
        forbidden_group::Vector{Vector{Int}} = Vector{Vector{Int}}(),
        atom_path::Vector{Int} = Vector{Int}()
)::Tuple{
        Vector{Tuple{Vector{Int}, Vector{Int}}}, Vector{Vector{Vector{Int}}}, Vector{Vector{Int}}
}

    # Get the node at the specified path
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    @match type begin
        :molecule => return get_ringbond_paths(solver, push!(copy(path), 1))
        :chain => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :(SMILES_combine_chain(bond,
                structure,
                chain)) => begin
                    structure_ringbonds, structure_forbidden,
                    ringbond_group = get_ringbond_paths(
                        solver,
                        push!(copy(path), 2),
                        forbidden_group = forbidden_group
                    )
                    chain_ringbonds, chain_forbidden,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 3), forbidden_group = ringbond_group
                    )
                    return vcat(structure_ringbonds, chain_ringbonds),
                    vcat(structure_forbidden, chain_forbidden),
                    []
                end
                :(structure * bond *
                chain) => begin
                    structure_ringbonds, structure_forbidden,
                    ringbond_group = get_ringbond_paths(
                        solver,
                        push!(copy(path), 1),
                        forbidden_group = forbidden_group
                    )
                    chain_ringbonds, chain_forbidden,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 3), forbidden_group = ringbond_group
                    )
                    return vcat(structure_ringbonds, chain_ringbonds),
                    vcat(structure_forbidden, chain_forbidden),
                    []
                end
                :(atom *
                ringbonds) => begin
                    ringbonds_paths, _,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 2), atom_path = push!(copy(path), 1)
                    )
                    forbidden_group = vcat(forbidden_group, [x[1] for x in ringbonds_paths])
                    ringbonds_group = copy([x[1] for x in ringbonds_paths])
                    return ringbonds_paths, [forbidden_group], ringbonds_group
                end
                :(structure * "-" *
                fragment_X_entry) => begin
                    structure_ringbonds, structure_forbidden,
                    ringbond_group = get_ringbond_paths(
                        solver,
                        push!(copy(path), 1),
                        forbidden_group = forbidden_group
                    )
                    fragment_ringbonds, fragment_forbidden,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 2), forbidden_group = ringbond_group
                    )
                    return vcat(structure_ringbonds, fragment_ringbonds),
                    vcat(structure_forbidden, fragment_forbidden),
                    []
                end
                _ => throw("Unknown chain rule: $rule")
            end
        end
        :structure => begin
            ringbond_paths, _,
            _ = get_ringbond_paths(
                solver, push!(copy(path), 2), atom_path = push!(copy(path), 1)
            )
            forbidden_group = vcat(forbidden_group, [x[1] for x in ringbond_paths])
            branches_ringbonds, branches_forbidden,
            _ = get_ringbond_paths(
                solver,
                push!(copy(path), 3),
                forbidden_group = [x[1] for x in ringbond_paths]
            )
            return vcat(ringbond_paths, branches_ringbonds),
            vcat([forbidden_group], branches_forbidden),
            [x[1] for x in ringbond_paths]
        end
        :ringbonds => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :("") => return [], [], []
                :(ringbond *
                ringbonds) => begin
                    ringbond_path = push!(copy(path), 1)
                    ringbond_paths, _,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 2), atom_path = atom_path
                    )
                    return vcat([(ringbond_path, atom_path)], ringbond_paths), [], []
                end
                _ => throw("Unknown ringbonds rule: $rule")
            end
        end
        :branches => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :("") => return [], [], []
                :(branch *
                branches) => begin
                    branch_ringbonds, branch_forbidden,
                    _ = get_ringbond_paths(
                        solver,
                        push!(copy(path), 1),
                        forbidden_group = forbidden_group
                    )
                    branches_ringbonds, branches_forbidden,
                    _ = get_ringbond_paths(
                        solver,
                        push!(copy(path), 2),
                        forbidden_group = forbidden_group
                    )
                    return vcat(branch_ringbonds, branches_ringbonds),
                    vcat(branch_forbidden, branches_forbidden),
                    []
                end
                _ => throw("Unknown branches rule: $rule")
            end
        end
        :branch => begin
            return get_ringbond_paths(
                solver, push!(copy(path), 2), forbidden_group = forbidden_group
            )
        end
        :fragment_X_entry || :starting_fragment => begin
            if !isfilled(node)
                all_ringbonds = Vector{Tuple{Vector{Int}, Vector{Int}}}()
                all_forbidden = Vector{Vector{Vector{Int}}}()
                for child_id in eachindex(node.children)
                    child_ringbonds, child_forbidden,
                    _ = get_ringbond_paths(solver, push!(copy(path), child_id))
                    append!(all_ringbonds, child_ringbonds)
                    append!(all_forbidden, child_forbidden)
                end
                return all_ringbonds, all_forbidden, []
            end
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            connects_to_single_atom = false
            all_ringbonds = Vector{Tuple{Vector{Int}, Vector{Int}}}()
            all_forbidden = Vector{Vector{Vector{Int}}}()
            atom_ringbonds = Vector{Tuple{Vector{Int}, Vector{Int}}}()
            branch_children = Vector{Int}()
            virtual_atom_count = 0
            children_count = 0

            args = rule isa Expr ? rule.args[2:end] : [rule]
            for arg in args
                if arg isa String
                    forbidden_group = postprocess_grouped_fragment_connection_points(
                        solver, path, forbidden_group, connects_to_single_atom,
                        all_ringbonds, all_forbidden, atom_ringbonds, branch_children)
                    atom_ringbonds = Vector{Tuple{Vector{Int}, Vector{Int}}}()
                    branch_children = Vector{Int}()
                    connects_to_single_atom = is_single_atom(arg)
                    virtual_atom_count += 1
                else
                    children_count += 1
                    child_node = get_node_at_location(
                        solver, push!(copy(path), children_count))
                    if !isfilled(child_node)
                        throw("Unfilled child node at path: $(push!(copy(path), children_count))")
                        continue
                    end
                    child_rule = solver.grammar.rules[HerbCore.get_rule(child_node)]
                    @match child_rule begin
                        :("-" * digit) => begin
                            push!(atom_ringbonds,
                                (push!(copy(path), children_count),
                                    push!(copy(path), virtual_atom_count)))
                        end
                        _ => push!(branch_children, children_count)
                    end
                end
            end

            if !isempty(atom_ringbonds) || !isempty(branch_children)
                postprocess_grouped_fragment_connection_points(
                    solver, path, forbidden_group, connects_to_single_atom,
                    all_ringbonds, all_forbidden, atom_ringbonds, branch_children)
            end

            return all_ringbonds, all_forbidden, []
        end
        :fragment_X_exit => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :("(-" * chain * ")") => begin
                    return get_ringbond_paths(
                        solver, push!(copy(path), 1), forbidden_group = forbidden_group)
                end
                # Push temporary ringbond path to check for parity of ringbonds in the molecule
                # The correct path will be updated when fragment_X_entry/starting_fragment is filled
                :("-" * digit) => begin
                    return [(copy(path), [-1])], [], []
                end
                _ => return get_ringbond_paths(
                    solver, push!(copy(path), 2), forbidden_group = forbidden_group)
            end
        end
        _ => throw("Unknown node type: $type")
    end
end

function post_ringbond_constraints!(
        solver::Solver, path::Vector{Int}, grammar_data::GrammarData
)
    # Get the ringbonds underneath the current path
    ringbond_paths, incompatible_groups, _ = get_ringbond_paths(solver, path)

    # Check if there are an even number of ringbonds
    if length(ringbond_paths) % 2 != 0
        HerbConstraints.set_infeasible!(solver)
    end

    # Post the local constraint with the incompatible pairs
    HerbConstraints.post!(
        solver,
        localUniformRingbonds(path, ringbond_paths, incompatible_groups, grammar_data)
    )
end
