function synthesize_molecules(
        atom_valences::OrderedDict{String, Int};
        settings::SynthesizerSettings = SynthesizerSettings(),
        starting_element::Union{Symbol, AbstractRuleNode} = :molecule,
        fragment_rules::OrderedDict{Int, Vector{Expr}} = OrderedDict{Int, Vector{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[]
)::Vector{Molecule}
    grammar = SMILES_grammar(
        atom_valences; settings = settings, fragment_rules = fragment_rules,
        starting_fragments = starting_fragments)
    iterator = get_iterator(settings, grammar, starting_element)

    candidates = Vector{Molecule}()
    start_time = time()
    for program in iterator
        molecule = interpret_molecule(program, grammar)
        push!(candidates, molecule)

        if check_stop_condition(settings, start_time, candidates, molecule)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end
