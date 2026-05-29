struct ContainsReactions <: AbstractGrammarConstraint
    required_molecules::OrderedDict{RequiredMolecule, Vector{Int}}
end

function ContainsReactions(
        grammar::ContextSensitiveGrammar,
        problem_required_molecules::Vector{RequiredMolecule},
        reactions::Vector{Reaction}
)
    reaction_to_rule = Dict{Any, Int}()
    for (i, rule) in enumerate(grammar.rules)
        reaction_to_rule[rule] = i
    end

    required_molecules::OrderedDict{
        RequiredMolecule, Vector{Int}} = OrderedDict{
        RequiredMolecule, Vector{Int}
    }()
    for required_molecule in problem_required_molecules
        required_molecules[required_molecule] = []

        for reaction in reactions
            if required_molecule.position == INPUT
                if any(input -> input[2] == required_molecule.molecule, reaction.inputs)
                    push!(
                        required_molecules[required_molecule],
                        reaction_to_rule[:($reaction)]
                    )
                end
            elseif required_molecule.position == OUTPUT
                if any(output -> output[2] == required_molecule.molecule, reaction.outputs)
                    push!(
                        required_molecules[required_molecule],
                        reaction_to_rule[:($reaction)]
                    )
                end
            else
                if any(input -> input[2] == required_molecule.molecule, reaction.inputs) ||
                   any(output -> output[2] == required_molecule.molecule, reaction.outputs)
                    push!(
                        required_molecules[required_molecule],
                        reaction_to_rule[:($reaction)]
                    )
                end
            end
        end
    end

    return ContainsReactions(required_molecules)
end

function HerbCore.is_domain_valid(
        constraint::ContainsReactions, grammar::ContextSensitiveGrammar
)
    # TODO: Implement the logic to check if the domain of the ContainsReactions constraint is valid
    return true
end

struct LocalContainsReactions <: AbstractLocalConstraint
    path::Vector{Int}
    required_molecules::OrderedDict{RequiredMolecule, Vector{Int}}
end

function HerbConstraints.on_new_node(
        solver::Solver, constraint::ContainsReactions, path::Vector{Int}
)
    if length(path) == 0 && solver isa UniformSolver

        # for (required_molecule, indices) in constraint.required_molecules
        #     println(indices)
        # end

        HerbConstraints.post!(
            solver, LocalContainsReactions(path, constraint.required_molecules)
        )
    end
end

function HerbConstraints.propagate!(solver::Solver, constraint::LocalContainsReactions)

    # Check if there are holes
    found = falses(length(constraint.required_molecules))
    found,
    holes = check_contains_reactions(
        solver, constraint, push!(copy(constraint.path), 1), found
    )

    # Check if the constraint is satisfied
    if all(found)
        HerbConstraints.deactivate!(solver, constraint)
    elseif holes == 0
        HerbConstraints.set_infeasible!(solver)
        return nothing
    elseif holes == 1
        # Get all the missing required molecule indices
        required_indices = nothing
        for (found_ind,
            (required_molecule, indices)) in enumerate(constraint.required_molecules)
            if !found[found_ind]
                if isnothing(required_indices)
                    required_indices = Set(indices)
                else
                    required_indices = intersect(required_indices, Set(indices))
                end
            end
        end
        if length(required_indices) == 0
            # If there are no matching indices, the constraint is infeasible
            HerbConstraints.set_infeasible!(solver)
            return nothing
        end
    end

    # println("Found: ", found)

end

function check_contains_reactions(
        solver::Solver, constraint::LocalContainsReactions, path::Vector{Int}, found, holes = 0
)
    node = get_node_at_location(solver, path)
    rule = solver.grammar.rules[get_rule(node)]

    @match rule begin
        :(vcat(reaction,
            reaction_list)) => begin
            node = get_node_at_location(solver, push!(copy(path), 1))
            if !isfilled(node)
                holes += 1
            else
                for (found_ind,
                    (required_molecule, indices)) in enumerate(constraint.required_molecules)
                    if get_rule(node) in indices
                        found[found_ind] = true
                    end
                end
            end
            return check_contains_reactions(
                solver, constraint, push!(copy(path), 2), found, holes
            )
        end

        :([reaction]) => begin
            node = get_node_at_location(solver, push!(copy(path), 1))
            if !isfilled(node)
                holes += 1
            else
                for (found_ind,
                    (required_molecule, indices)) in enumerate(constraint.required_molecules)
                    if get_rule(node) in indices
                        found[found_ind] = true
                    end
                end
            end
            return found, holes
        end

        :(Vector{Reaction}()) => begin
            return found, holes
        end

        _ => throw("Unexpected rule: $rule")
    end
end

function is_valid(candidate::ReactionNetwork, constraint::ContainsReactions)
    required_molecules = keys(constraint.required_molecules)
    found = falses(length(required_molecules))

    for reaction in candidate.reactions
        for (_, molecule) in reaction.inputs
            for (i, required_molecule) in enumerate(required_molecules)
                if (
                    required_molecule.position == INPUT ||
                    required_molecule.position == UNKNOWN
                ) && molecule == required_molecule.molecule
                    found[i] = true
                end
            end
        end

        for (_, molecule) in reaction.outputs
            for (i, required_molecule) in enumerate(required_molecules)
                if (
                    required_molecule.position == OUTPUT ||
                    required_molecule.position == UNKNOWN
                ) && molecule == required_molecule.molecule
                    found[i] = true
                end
            end
        end
    end

    return all(found)
end
