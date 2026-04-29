function synthesize_molecules(
        atoms::Vector{Atom}, settings::SynthesizerSettings, starting_element::Union{Symbol, AbstractRuleNode} = :molecule
)::Vector{Molecule}
    grammar = SMILES_grammar(atoms; settings = settings)
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
