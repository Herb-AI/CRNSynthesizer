function synthesize_molecules(
        atoms::Vector{Atom}, settings::SynthesizerSettings,
        starting_element::Union{Symbol, AbstractRuleNode} = :molecule; fragment_rules::Dict{
            Int, Vector{Expr}} = Dict{Int, Vector{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[],
        known_molecules::Vector{Molecule} = Molecule[]
)::Vector{Molecule}
    grammar = !haskey(settings.options, :use_bottom_up) ||
              !settings.options[:use_bottom_up] ?
              SMILES_grammar(atoms; settings = settings, fragment_rules = fragment_rules,
        starting_fragments = starting_fragments) :
              bu_SMILES_grammar(
        atoms; settings = settings, fragment_rules = fragment_rules,
        starting_fragments = starting_fragments)
    iterator = get_iterator(settings, grammar, starting_element)

    candidates = Vector{Molecule}()
    start_time = time()
    for program in iterator
        if !haskey(settings.options, :use_bottom_up) || !settings.options[:use_bottom_up]
            molecule = interpret_molecule(program, grammar)
            push!(candidates, molecule)
        else
            molecule = bu_interpret_molecule(program, grammar)
            if is_valid(molecule)
                push!(candidates, molecule)
            end
        end

        if check_stop_condition(settings, start_time, candidates, molecule)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end
