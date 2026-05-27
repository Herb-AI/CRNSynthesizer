
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
