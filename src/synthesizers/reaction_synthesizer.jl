
function synthesize_reactions(
        atoms::Vector{String}, settings::SynthesizerSettings
)::Vector{Reaction}
    grammar = reaction_grammar(atoms; settings = settings)
    iterator = get_iterator(settings, grammar, :reaction)

    candidates = Vector{Reaction}()
    start_time = time()
    for program in iterator
        reaction = interpret_reaction(program, grammar)

        push!(candidates, reaction)

        if check_stop_condition(settings, start_time, candidates, reaction)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_reactions(
        molecules::Vector{Molecule}, settings::SynthesizerSettings;
        known_molecules::Vector{Molecule} = Molecule[]
)::Vector{Reaction}
    metric_name = get(settings.options, :similarity_metric, :none)
    sorted_molecules = (isempty(known_molecules) || metric_name == :none) ? molecules :
                       sort_molecules_by_similarity(molecules, known_molecules; metric_name=metric_name)
    grammar = reaction_grammar(sorted_molecules; settings = settings)
    iterator = get_iterator(settings, grammar, :reaction)

    candidates = Vector{Reaction}()
    start_time = time()
    for program in iterator
        reaction = interpret_reaction(program, grammar)

        push!(candidates, reaction)

        if check_stop_condition(settings, start_time, candidates, reaction)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_reactions(
        problem::ProblemDefinition,
        molecule_settings::SynthesizerSettings,
        reaction_settings::SynthesizerSettings;
        initial_molecules_count::Int = 10,
        fragment_rules::Dict{Int, Set{Expr}} = Dict{Int, Set{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[]
)
    start_time = time()

    molecules = OrderedSet{Molecule}()
    metric_name = get(reaction_settings.options, :similarity_metric, :none)

    for molecule in problem.known_molecules
        push!(molecules, molecule)
    end
    check_stop_condition(
        molecule_settings, start_time, molecules, nothing; check_all_candidates = true
    )

    molecule_grammar = SMILES_grammar(problem.atom_valences; settings = molecule_settings, fragment_rules = fragment_rules, starting_fragments = starting_fragments)
    molecule_iterator = get_iterator(molecule_settings, molecule_grammar, :molecule)

    reactions = OrderedSet{Reaction}()
    molecule_score_cache = Dict{Molecule, Float64}()

    stop_condition = false
    for molecule_program in molecule_iterator
        molecule = interpret_molecule(molecule_program, molecule_grammar)

        if !(molecule in problem.known_molecules)
            push!(molecules, molecule)
        end

        if length(molecules) < initial_molecules_count &&
           !check_stop_condition(molecule_settings, start_time, molecules, molecule)
            continue
        end

        reactions_grammar = reaction_grammar(
            sort_molecules_by_similarity(collect(molecules), problem.known_molecules; metric_name=metric_name, cache=molecule_score_cache); settings = reaction_settings
        )
        reaction_iterator = get_iterator(reaction_settings, reactions_grammar, :reaction)
        for reaction_program in reaction_iterator
            reaction = interpret_reaction(reaction_program, reactions_grammar)
            if reaction in reactions
                continue
            end
            push!(reactions, reaction)

            if check_stop_condition(reaction_settings, start_time, reactions, reaction)
                stop_condition = true
                break
            end
        end

        if stop_condition ||
           check_stop_condition(molecule_settings, start_time, molecules, molecule)
            break
        end
    end

    return collect(reactions), molecules
end

