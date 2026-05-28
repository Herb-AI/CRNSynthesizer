function reaction_grammar(;
        settings::SynthesizerSettings = SynthesizerSettings(),
        add_required_rule::Bool = false,
        complete_grammar::Bool = true
)
    grammar = @csgrammar begin
        reaction = Reaction(molecule_list, molecule_list)
        molecule_list = vcat([molecule], molecule_list)
        molecule_list = Vector{Molecule}()
    end

    if add_required_rule
        add_rule!(grammar, :(molecule_list = vcat([required_molecule], molecule_list)))
    end

    # Makes sure that there is at least one molecule in the reaction
    addconstraint!(grammar, Forbidden(@c_rulenode 1{3, c}))
    addconstraint!(grammar, Forbidden(@c_rulenode 1{c, 3}))

    # Forbids reactions to be the same on both sides
    addconstraint!(grammar, Forbidden(@c_rulenode 1{a, a}))

    # Makes the molecule list ordered to break symmetries
    if !(
        haskey(settings.options, :disable_ordered_molecule_list) &&
        settings.options[:disable_ordered_molecule_list]
    )
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))

        if add_required_rule
            # Makes sure that the required molecule is always at the front of the list
            addconstraint!(grammar, Ordered((@c_rulenode 4{a, 4{b, c}}), [:a, :b]))
            addconstraint!(grammar, Ordered((@c_rulenode 2{a, 4{b, c}}), [:a, :b]))
            addconstraint!(grammar, Ordered((@c_rulenode 4{a, 2{b, c}}), [:a, :b]))
        end
    end

    # Make sure the reactions are atom balanced
    if !(
        haskey(settings.options, :disable_balanced_reaction) &&
        settings.options[:disable_balanced_reaction]
    )
        addconstraint!(grammar, BalancedReaction(; complete_grammar = complete_grammar))
    end

    return grammar
end

function partial_reaction_grammar(partial_reaction::Reaction)
    inputs = Vector{Molecule}()
    for (count, mol) in partial_reaction.inputs
        for _ in 1:count
            push!(inputs, mol)
        end
    end
    outputs = Vector{Molecule}()
    for (count, mol) in partial_reaction.outputs
        for _ in 1:count
            push!(outputs, mol)
        end
    end

    grammar = @csgrammar begin
        reaction = Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules))
    end
    grammar = add_rule!(grammar, :(input_molecules = $inputs))
    grammar = add_rule!(grammar, :(output_molecules = $outputs))

    return grammar
end

function reaction_grammar(
        molecules::Vector{Molecule}; settings::SynthesizerSettings = SynthesizerSettings(),
        partial_reaction::Union{Nothing, Reaction} = nothing
)
    grammar = reaction_grammar(; settings = settings, complete_grammar = false)

    for molecule in molecules
        add_rule!(grammar, :(molecule = $molecule))
    end

    if !isnothing(partial_reaction)
        new_grammar = partial_reaction_grammar(partial_reaction)
        merge_grammars!(new_grammar, grammar)
        grammar = new_grammar

        # Add Forbidden constraints for empty input/output boundary cases after merging
        rule_expr = :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules)))
        idx = findfirst(r -> r == rule_expr, grammar.rules)
        empty_list_idx = findfirst(r -> r == :(Vector{Molecule}()), grammar.rules)
        if !isnothing(idx)
            addconstraint!(grammar, Forbidden(RuleNode(idx, [VarNode(:a), VarNode(:b), VarNode(:a), VarNode(:c)])))
            if !isnothing(empty_list_idx)
                if isempty(partial_reaction.inputs)
                    addconstraint!(grammar, Forbidden(RuleNode(idx, [RuleNode(empty_list_idx), VarNode(:a), VarNode(:b), VarNode(:c)])))
                end
                if isempty(partial_reaction.outputs)
                    addconstraint!(grammar, Forbidden(RuleNode(idx, [VarNode(:a), VarNode(:b), RuleNode(empty_list_idx), VarNode(:c)])))
                end
            end
        end
    end

    return grammar
end

function reaction_grammar(
        atoms::Vector{String}; settings::SynthesizerSettings = SynthesizerSettings()
)
    grammar = reaction_grammar(; settings = settings, complete_grammar = true)
    merge_grammars!(grammar, SMILES_grammar(atoms; settings = settings))

    return grammar
end
