function network_grammar(; settings::SynthesizerSettings = SynthesizerSettings())
    grammar = @csgrammar begin
        network = ReactionNetwork(reaction_list)
        reaction_list = vcat(reaction, reaction_list)
        reaction_list = Vector{Reaction}()
    end

    # Forbids reactions to occur more than once in the network
    addconstraint!(grammar, Forbidden(@c_rulenode 2{a, 2{a, c}}))

    if !(
        haskey(settings.options, :disable_ordered_reaction_list) &&
        settings.options[:disable_ordered_reaction_list]
    )
        # Makes the reaction list ordered to break symmetries
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end

    return grammar
end

function network_grammar(
        atoms::Vector{String};
        problem::ProblemDefinition = ProblemDefinition(),
        settings::SynthesizerSettings = SynthesizerSettings()
)
    grammar = network_grammar()
    merge_grammars!(grammar, reaction_grammar(; settings = settings))
    merge_grammars!(grammar, SMILES_grammar(atoms))

    # required = problem.required_molecules
    # required_rules = Vector{Tuple{Int, ReactionPosition}}()
    # for (i, required_molecule) in enumerate(required)
    #     molecule = required_molecule.molecule
    #     add_rule!(grammar, :(required_molecule = $molecule))
    #     rule = findfirst(==(:($molecule)), grammar.rules)
    #     push!(required_rules, (rule, required_molecule.position))
    # end

    # if !(haskey(settings.options, :disable_contains_molecules) && settings.options[:disable_contains_molecules])
    #     add_rule!(grammar, :(molecule_list = vcat([required_molecule], molecule_list)))
    #     addconstraint!(grammar, ContainsMolecules(required_rules))
    # end

    return grammar
end

function network_grammar(
        molecules::Vector{Molecule};
        problem::ProblemDefinition = ProblemDefinition(),
        required_molecules::Vector{Molecule} = Vector{Molecule}(),
        settings::SynthesizerSettings = SynthesizerSettings(),
        check_required::Bool = true
)
    grammar = network_grammar(; settings = settings)
    merge_grammars!(
        grammar, reaction_grammar(; add_required_rule = true, complete_grammar = false)
    )

    for molecule in molecules
        add_rule!(grammar, :(molecule = $molecule))
    end

    required = problem.required_molecules
    required_rules = Vector{Tuple{Int, ReactionPosition}}()
    req_molecules = Vector{Tuple{Molecule, ReactionPosition}}()
    for (i, required_molecule) in enumerate(required)
        molecule = required_molecule.molecule
        add_rule!(grammar, :(required_molecule = $molecule))
        rule = findfirst(==(:($molecule)), grammar.rules)
        push!(required_rules, (rule, required_molecule.position))
        push!(req_molecules, (required_molecule.molecule, required_molecule.position))
    end


    for molecule in required_molecules
        add_rule!(grammar, :(required_molecule = $molecule))
        rule = findfirst(==(:($molecule)), grammar.rules)
        push!(required_rules, (rule, UNKNOWN))
        push!(req_molecules, (molecule, UNKNOWN))
    end

    if length(required) == 0 && length(required_molecules) == 0 && check_required
        addconstraint!(grammar, Forbidden(@c_rulenode 7))
    elseif !(
        haskey(settings.options, :disable_contains_molecules) &&
        settings.options[:disable_contains_molecules]
    )
        addconstraint!(grammar, ContainsMolecules(req_molecules, required_rules))
    end

    return grammar
end

function network_grammar(
        reactions::Vector{Reaction};
        problem::ProblemDefinition = ProblemDefinition(),
        settings::SynthesizerSettings = SynthesizerSettings()
)
    grammar = network_grammar(; settings = settings)

    for reaction in reactions[1:(end - 1)]
        c_add_rule!(grammar, :(reaction = $reaction))
    end
    reaction = reactions[end]
    c_add_rule!(grammar, :(reaction = $reaction); recalculate = true)

    reaction_to_rule = Dict{Any, Int}()
    for (i, rule) in enumerate(grammar.rules)
        reaction_to_rule[rule] = i
    end

    required_molecules::Dict{
        RequiredMolecule, Vector{Int}} = Dict{
        RequiredMolecule, Vector{Int}
    }()
    for required_molecule in problem.required_molecules
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

    if !(
        haskey(settings.options, :disable_contains_molecules) &&
        settings.options[:disable_contains_molecules]
    )
        addconstraint!(grammar, ContainsReactions(required_molecules))
    end

    return grammar
end

function network_grammar(
        problem::ProblemDefinition; settings::SynthesizerSettings = SynthesizerSettings()
)
    return network_grammar(get_atoms(problem); problem = problem, settings = settings)
end
