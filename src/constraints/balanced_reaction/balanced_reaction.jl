struct BalancedReaction <: AbstractGrammarConstraint
    complete_grammar::Bool

    function BalancedReaction(; complete_grammar::Bool = true)
        return new(complete_grammar)
    end
end

function HerbCore.is_domain_valid(
        constraint::BalancedReaction, grammar::ContextSensitiveGrammar
)
    # TODO: Check if the grammar is valid for balanced reactions
    return true
end

function HerbCore.update_rule_indices!(constraint::BalancedReaction, n_rules::Integer)
    return nothing
end

function HerbCore.update_rule_indices!(
        constraint::BalancedReaction,
        n_rules::Integer,
        mapping::AbstractDict{<:Integer, <:Integer},
        constraints::Vector{<:AbstractConstraint}
)
    return nothing
end

struct AtomCounts
    atom_counts::OrderedDict{String, Int}
    hash::UInt
end
function AtomCounts(atom_counts::OrderedDict{String, Int})
    return AtomCounts(atom_counts, hash(atom_counts))
end

function Base.show(io::IO, a::AtomCounts)
    atom_counts = a.atom_counts
    atom_counts_str = join(
        [string(k) * (v == 1 ? "" : string(v)) for (k, v) in atom_counts], " + "
    )
    print(io, atom_counts_str)
end

import Base.==
function ==(a::AtomCounts, b::AtomCounts)
    return a.hash == b.hash
end

import Base.hash
function hash(a::AtomCounts, h::UInt)
    return a.hash
end

struct ReactionPossibility
    rules::Vector{Int}
    atom_counts::AtomCounts
    hash::UInt
end
function ReactionPossibility(rules::Vector{Int}, atom_counts::OrderedDict{String, Int})
    atom_counts = AtomCounts(atom_counts)
    return ReactionPossibility(rules, atom_counts, hash(rules, hash(atom_counts)))
end
function ReactionPossibility(rules::Vector{Int}, atom_counts::AtomCounts)
    return ReactionPossibility(rules, atom_counts, hash(rules, hash(atom_counts)))
end

function Base.show(io::IO, a::ReactionPossibility)
    rules_str = join([string(r) for r in a.rules], ", ")
    atom_counts_str = string(a.atom_counts)
    print(io, "Rules: [", rules_str, "], Atom Counts: ", atom_counts_str)
end

import Base.==
function ==(a::ReactionPossibility, b::ReactionPossibility)
    return a.hash == b.hash
end

import Base.hash
function hash(a::ReactionPossibility, h::UInt)
    return a.hash
end

struct LocalBalancedReaction <: AbstractLocalConstraint
    path::Vector{Int}
    input_molecule_paths::Vector{Vector{Int}}
    output_molecule_paths::Vector{Vector{Int}}
    rule_to_atoms::Dict{Int, Dict{String, Int}}
    left_possibilities::Vector{ReactionPossibility}
    right_possibilities::Vector{ReactionPossibility}
end

import Base.==
function ==(a::LocalBalancedReaction, b::LocalBalancedReaction)
    return a.path == b.path &&
           a.input_molecule_paths == b.input_molecule_paths &&
           a.output_molecule_paths == b.output_molecule_paths
end

import Base.hash
function hash(a::LocalBalancedReaction, h::UInt)
    return hash(a.path, h) +
           hash(a.input_molecule_paths, h) +
           hash(a.output_molecule_paths, h)
end

function get_atom_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :atom
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_atom_paths(solver, next_path))
    end
    return result
end

function get_molecule_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :molecule
        return [path]
    end
    if type == :required_molecule
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_molecule_paths(solver, next_path))
    end
    return result
end

function get_prefilled_atoms(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    if node isa Hole
        return OrderedDict{String, Int}(), OrderedDict{String, Int}()
    end
    rule = HerbCore.get_rule(node)
    rule_expr = solver.grammar.rules[rule]

    @match rule_expr begin
        :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules))) =>
            begin
                # Safety check: ensure children are fully instantiated before extracting prefilled atoms
                if !isfilled(node.children[2]) || !isfilled(node.children[4])
                    return OrderedDict{String, Int}(), OrderedDict{String, Int}()
                end

                input_rule = HerbCore.get_rule(node.children[2])
                inputs_vec = solver.grammar.rules[input_rule]

                output_rule = HerbCore.get_rule(node.children[4])
                outputs_vec = solver.grammar.rules[output_rule]

                left_atoms = OrderedDict{String, Int}()
                for mol in inputs_vec
                    left_atoms = mergewith(+, left_atoms, count_atoms(mol))
                end

                right_atoms = OrderedDict{String, Int}()
                for mol in outputs_vec
                    right_atoms = mergewith(+, right_atoms, count_atoms(mol))
                end

                return left_atoms, right_atoms
            end
        _ => return OrderedDict{String, Int}(), OrderedDict{String, Int}()
    end
end

function post_reaction_constraints!(solver::Solver, reaction_paths::Vector{Vector{Int}})
    reaction_constraints = []
    for (i, path) in enumerate(reaction_paths)
        node = get_node_at_location(solver, path)
        rule = HerbCore.get_rule(node)
        rule_expr = solver.grammar.rules[rule]

        is_partial = @match rule_expr begin
            :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules))) =>
                true
            _ => false
        end

        input_paths = get_molecule_paths(solver, push!(copy(path), 1))
        output_paths = get_molecule_paths(solver, push!(copy(path), is_partial ? 3 : 2))

        rule_to_atoms = Dict{Int, Dict{String, Int}}()
        molecule_rules = findall(solver.grammar.domains[:molecule])

        for rule in molecule_rules
            rule_to_atoms[rule] = count_atoms(solver.grammar.rules[rule])
        end

        prefilled_left, prefilled_right = get_prefilled_atoms(solver, path)

        possible_left = get_possibilities(
            solver, input_paths; rule_counts = rule_to_atoms, prefilled_atoms = prefilled_left
        )
        possible_right = get_possibilities(
            solver, output_paths; rule_counts = rule_to_atoms, prefilled_atoms = prefilled_right
        )

        left_length = length(possible_left)
        right_length = length(possible_right)

        # Do a pre pass with the intersection
        l = [x.atom_counts for x in possible_left]
        r = [x.atom_counts for x in possible_right]
        intersection = intersect(l, r)
        filter!(x -> x.atom_counts in intersection, possible_left)
        filter!(x -> x.atom_counts in intersection, possible_right)

        if left_length == 0 || right_length == 0
            HerbConstraints.set_infeasible!(solver)
            return nothing
        end

        balanced_reaction = LocalBalancedReaction(
            path, input_paths, output_paths, rule_to_atoms, possible_left, possible_right
        )
        HerbConstraints.post!(solver, balanced_reaction)
        push!(reaction_constraints, balanced_reaction)
    end

    return reaction_constraints
end

function HerbConstraints.on_new_node(
        solver::Solver, constraint::BalancedReaction, path::Vector{Int}
)
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    if constraint.complete_grammar
        if solver isa GenericSolver
            if type == :reaction
                # TODO: check the starting vectors for the input and output paths
                HerbConstraints.post!(solver, LocalGenericBalancedReaction(path))
            end
        else
            # Uniform solver
            node = get_node_at_location(solver, path)
            type = get_node_type(solver.grammar, node)
            if type == :reaction
                input_paths = get_atom_paths(solver, push!(copy(path), 1))
                output_paths = get_atom_paths(solver, push!(copy(path), 2))
                HerbConstraints.post!(
                    solver, LocalUniformBalancedReaction(path, input_paths, output_paths)
                )
            end
        end
    else
        if solver isa UniformSolver && type == :reaction
            node = get_node_at_location(solver, path)
            rule = HerbCore.get_rule(node)
            rule_expr = solver.grammar.rules[rule]
            is_partial = @match rule_expr begin
                :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules))) =>
                    true
                _ => false
            end

            input_paths = get_molecule_paths(solver, push!(copy(path), 1))
            output_paths = get_molecule_paths(solver, push!(copy(path), is_partial ? 3 : 2))

            rule_to_atoms = Dict{Int, Dict{String, Int}}()
            molecule_rules = findall(solver.grammar.domains[:molecule])

            for rule in molecule_rules
                rule_to_atoms[rule] = count_atoms(solver.grammar.rules[rule])
            end

            prefilled_left, prefilled_right = get_prefilled_atoms(solver, path)

            possible_left = get_possibilities(
                solver, input_paths; rule_counts = rule_to_atoms, prefilled_atoms = prefilled_left
            )
            possible_right = get_possibilities(
                solver, output_paths; rule_counts = rule_to_atoms,
                prefilled_atoms = prefilled_right
            )

            left_length = length(possible_left)
            right_length = length(possible_right)

            # Do a pre pass with the intersection
            l = [x.atom_counts for x in possible_left]
            r = [x.atom_counts for x in possible_right]
            intersection = intersect(l, r)
            filter!(x -> x.atom_counts in intersection, possible_left)
            filter!(x -> x.atom_counts in intersection, possible_right)

            if left_length == 0 || right_length == 0
                HerbConstraints.set_infeasible!(solver)
                return nothing
            end

            HerbConstraints.post!(
                solver,
                LocalBalancedReaction(
                    path,
                    input_paths,
                    output_paths,
                    rule_to_atoms,
                    possible_left,
                    possible_right
                )
            )
        end
    end
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::LocalBalancedReaction, path::Vector{Int}
)
    if path in constraint.input_molecule_paths || path in constraint.output_molecule_paths
        # push!(constraint.updated_paths, path)
        return true
    end

    return false
end

function count_atoms(solver::Solver, paths::Vector{Vector{Int}})::OrderedDict{String, Int}
    result = OrderedDict{String, Int}()
    for path in paths
        node = get_node_at_location(solver, path)
        type = get_node_type(solver.grammar, node)

        if type == :molecule
            if isfilled(node)
                rule = solver.grammar.rules[HerbCore.get_rule(node)]
                result = mergewith(+, result, count_atoms(rule))
            else
                return nothing
            end
        end
    end

    return result
end

function get_possibilities(
        solver::Solver,
        paths::Vector{Vector{Int}};
        rule_counts::Dict{Int, Dict{String, Int}} = Dict{Int, Dict{String, Int}}(),
        prefilled_atoms::OrderedDict{String, Int} = OrderedDict{String, Int}()
)

    # Build atom indices dynamically based on the elements and charges present in the grammar
    all_keys = OrderedSet{String}()
    for (rule, counts) in rule_counts
        union!(all_keys, keys(counts))
    end
    if isempty(all_keys)
        molecule_rules = findall(solver.grammar.domains[:molecule])
        for rule in molecule_rules
            union!(all_keys, keys(count_atoms(solver.grammar.rules[rule])))
        end
    end
    union!(all_keys, keys(prefilled_atoms))

    sorted_atoms = sort(collect(all_keys))
    atom_indices = OrderedDict{String, Int}(atom => i for (i, atom) in enumerate(sorted_atoms))
    n_atoms = max(length(sorted_atoms), 1)

    # Initialize with prefilled atom counts
    initial_counts = zeros(Int, n_atoms)
    for (atom, count) in prefilled_atoms
        initial_counts[atom_indices[atom]] = count
    end
    current_options = [(initial_counts, Int[])]

    # Process each path
    for path in paths
        node = get_node_at_location(solver, path)
        rules = get_rules(node)

        new_options = Vector{Tuple{Vector{Int}, Vector{Int}}}()

        for (atom_counts, rules_used) in current_options
            for rule_ind in rules
                # Enforce non-decreasing order
                if isempty(rules_used) || rule_ind >= rules_used[end]
                    # Get or calculate atom counts for this rule
                    if !haskey(rule_counts, rule_ind)
                        count = count_atoms(solver.grammar.rules[rule_ind])
                        rule_counts[rule_ind] = count
                    else
                        count = rule_counts[rule_ind]
                    end

                    # Merge atom counts
                    merged_count::Vector{Int} = copy(atom_counts)
                    for (atom, count) in count
                        merged_count[atom_indices[atom]] += count
                    end

                    merged_rules::Vector{Int} = vcat(rules_used, rule_ind)

                    push!(new_options, (merged_count, merged_rules))
                end
            end
        end

        current_options = new_options
    end

    # Convert to final result format
    results = ReactionPossibility[]

    for (counts, rules) in current_options
        atom_counts = OrderedDict{String, Int}()
        for (i, count) in enumerate(counts)
            if count > 0
                atom_counts[sorted_atoms[i]] = count
            end
        end
        push!(results, ReactionPossibility(rules, atom_counts))
    end

    return results
end

function HerbConstraints.propagate!(solver::Solver, constraint::LocalBalancedReaction)
    if solver isa GenericSolver
        return nothing
    end

    possible_left = copy(constraint.left_possibilities)
    for (index, path) in enumerate(constraint.input_molecule_paths)
        node = get_node_at_location(solver, path)
        rules = OrderedSet{Int}(get_rules(node))
        filter!(x -> x.rules[index] in rules, possible_left)
    end

    possible_right = copy(constraint.right_possibilities)
    for (index, path) in enumerate(constraint.output_molecule_paths)
        node = get_node_at_location(solver, path)
        rules = OrderedSet{Int}(get_rules(node))
        filter!(x -> x.rules[index] in rules, possible_right)
    end

    l = Set([x.atom_counts for x in possible_left])
    r = Set([x.atom_counts for x in possible_right])
    intersection = intersect(l, r)

    if isempty(intersection)
        HerbConstraints.set_infeasible!(solver)
        return nothing
    end

    filter!(x -> x.atom_counts in intersection, possible_left)
    filter!(x -> x.atom_counts in intersection, possible_right)

    if isempty(possible_left) || isempty(possible_right)
        HerbConstraints.set_infeasible!(solver)
        return nothing
    end

    fixed_left = Int[]
    for (i, path) in enumerate(constraint.input_molecule_paths)
        possibilities = [x.rules[i] for x in possible_left]
        # println("Possibilities for input ", i, ": ", possibilities)
        c_remove_all_but!(solver, path, possibilities, false)
        if !isfeasible(solver)
            return nothing
        end

        node = get_node_at_location(solver, path)
        rules = get_rules(node)
        if length(possibilities) == 1 || length(rules) == 1
            push!(fixed_left, possibilities[1])
        end
    end

    fixed_right = Int[]
    for (i, path) in enumerate(constraint.output_molecule_paths)
        possibilities = [x.rules[i] for x in possible_right]
        # println("Possibilities for output ", i, ": ", possibilities)
        c_remove_all_but!(solver, path, possibilities, false)

        node = get_node_at_location(solver, path)
        rules = get_rules(node)
        if length(possibilities) == 1 || length(rules) == 1
            push!(fixed_right, possibilities[1])
        end
    end

    # any item that is fixed in the right must not be on the left
    for path in constraint.input_molecule_paths
        c_remove!(solver, path, fixed_right, false)
    end

    # any item that is fixed in the left must not be on the right
    for path in constraint.output_molecule_paths
        c_remove!(solver, path, fixed_left, false)
    end

    if (length(possible_left) == 1)
        for (i, path) in enumerate(constraint.input_molecule_paths)
            c_remove_all_but!(solver, path, [possible_left[1].rules[i]], false)
        end
    end

    if (length(possible_right) == 1)
        for (i, path) in enumerate(constraint.output_molecule_paths)
            c_remove_all_but!(solver, path, [possible_right[1].rules[i]], false)
        end
    end

    if !isfeasible(solver)
        return nothing
    end
end

function is_valid(candidate::Reaction, constraint::BalancedReaction)
    input_counts = Dict{String, Int}()
    output_counts = Dict{String, Int}()

    for (num, molecule) in candidate.inputs
        for (atom, count) in count_atoms(molecule)
            input_counts[atom] = get(input_counts, atom, 0) + count * num
        end
    end

    for (num, molecule) in candidate.outputs
        for (atom, count) in count_atoms(molecule)
            output_counts[atom] = get(output_counts, atom, 0) + count * num
        end
    end

    filter!(p -> p.second != 0, input_counts)
    filter!(p -> p.second != 0, output_counts)

    if input_counts != output_counts
        return false
    end

    return true
end
