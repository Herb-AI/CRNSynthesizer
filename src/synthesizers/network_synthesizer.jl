function synthesize_networks(
        atoms::Vector{String}, settings::SynthesizerSettings;
        problem = ProblemDefinition()
)::Vector{ReactionNetwork}
    grammar = network_grammar(atoms; problem = problem, settings = settings)
    iterator = get_iterator(settings, grammar, :network)
    interpreter = x -> interpret_network(x, grammar)

    candidates = Vector{ReactionNetwork}()
    find_programs!(iterator, settings, interpreter, candidates)

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_networks(
        molecules::Vector{Molecule},
        settings::SynthesizerSettings;
        problem = ProblemDefinition(),
        job = nothing
)::Vector{ReactionNetwork}
    molecules = filter(molecules) do molecule
        return !(molecule in problem.known_molecules)
    end
    metric_name = get(settings.options, :similarity_metric, :none)
    grammar = network_grammar(sort_molecules_by_similarity(molecules, problem.known_molecules; metric_name=metric_name); problem = problem, settings = settings)

    iterator = get_iterator(settings, grammar, :network)
    interpreter = x -> interpret_network(x, grammar)

    candidates = Vector{ReactionNetwork}()
    find_programs!(iterator, settings, interpreter, candidates)

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_networks(
        reactions::Vector{Reaction}, settings::SynthesizerSettings;
        problem = ProblemDefinition()
)::Vector{ReactionNetwork}
    metric_name = get(settings.options, :similarity_metric, :none)
    grammar = network_grammar(sort_reactions_by_similarity(reactions, problem.known_molecules; metric_name=metric_name),; problem = problem, settings = settings)
    iterator = get_iterator(settings, grammar, :network)
    interpreter = x -> interpret_network(x, grammar)

    candidates = Vector{ReactionNetwork}()
    find_programs!(iterator, settings, interpreter, candidates)

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_networks(
        problem::ProblemDefinition, settings::SynthesizerSettings
)::Vector{ReactionNetwork}
    return synthesize_networks(get_atoms(problem), settings; problem = problem)
end

function synthesize_networks(
        problem::ProblemDefinition,
        molecule_settings::SynthesizerSettings,
        network_settings::SynthesizerSettings;
        initial_molecules_count::Int = 100,
        pbar = nothing,
        fragment_rules::Dict{Int, Set{Expr}} = Dict{Int, Set{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[]
)
    molecule_job = nothing
    network_job = nothing
    if !isnothing(pbar)
        molecule_job = addjob!(pbar)
        network_job = addjob!(pbar)
    end

    start_time = time()
    molecules = Set{Molecule}()
    metric_name = get(network_settings.options, :similarity_metric, :none)
    molecule_grammar = SMILES_grammar(get_atoms(problem); fragment_rules = fragment_rules, starting_fragments = starting_fragments)
    molecule_iterator = get_iterator(molecule_settings, molecule_grammar, :molecule)

    networks = Vector{ReactionNetwork}()

    stop_condition = false
    for molecule_program in molecule_iterator
        molecule = interpret_molecule(molecule_program, molecule_grammar)

        if !(molecule in problem.known_molecules)
            push!(molecules, molecule)

            if !isnothing(pbar)
                update!(molecule_job)
                yield()
            end
        end

        if length(molecules) < initial_molecules_count &&
           !check_stop_condition(molecule_settings, start_time, molecules, molecule)
            continue
        end

        networks_grammar = network_grammar(
            sort_molecules_by_similarity(collect(molecules), problem.known_molecules; metric_name=metric_name); problem = problem, settings = network_settings
        )
        network_iterator = get_iterator(network_settings, networks_grammar, :network)
        for network_program in network_iterator
            network = interpret_network(network_program, networks_grammar)

            push!(networks, network)

            if !isnothing(pbar)
                update!(network_job)
                yield()
            end

            if check_stop_condition(network_settings, start_time, networks, network)
                stop_condition = true
                break
            end
        end

        if stop_condition ||
           check_stop_condition(molecule_settings, start_time, molecules, molecule)
            break
        end
    end

    if !isnothing(pbar)
        removejob!(pbar, molecule_job)
        removejob!(pbar, network_job)
    end

    if haskey(network_settings.options, :unique_candidates) &&
       network_settings.options[:unique_candidates]
        networks = unique(networks)
    end

    return networks, molecules
end

function synthesize_networks_2(
        problem::ProblemDefinition,
        reaction_settings::SynthesizerSettings,
        network_settings::SynthesizerSettings;
        initial_reactions_count::Int = 100,
        pbar = nothing
)
    reaction_job = nothing
    network_job = nothing
    if !isnothing(pbar)
        reaction_job = addjob!(pbar)
        network_job = addjob!(pbar)
    end

    start_time = time()

    reactions = Set{Reaction}()
    metric_name = get(network_settings.options, :similarity_metric, :simpson)
    combine_method = get(network_settings.options, :similarity_combine, :multiplicative)

    reactions_grammar = reaction_grammar(get_atoms(problem); settings = reaction_settings)
    reaction_iterator = get_iterator(reaction_settings, reactions_grammar, :reaction)

    networks = Vector{ReactionNetwork}()

    stop_condition = false
    for reaction_program in reaction_iterator
        reaction = interpret_reaction(reaction_program, reactions_grammar)
        if reaction in reactions
            continue
        end
        push!(reactions, reaction)

        if !isnothing(pbar)
            update!(reaction_job)
            yield()
        end
        if length(reactions) < initial_reactions_count &&
           !check_stop_condition(reaction_settings, start_time, reactions, reaction)
            continue
        end

        networks_grammar = network_grammar(
            sort_reactions_by_similarity(collect(reactions), problem.known_molecules; metric_name=metric_name, combine=combine_method); problem = problem, settings = network_settings
        )
        network_iterator = get_iterator(network_settings, networks_grammar, :network)
        for network_program in network_iterator
            network = interpret_network(network_program, networks_grammar)

            push!(networks, network)

            if !isnothing(pbar)
                update!(network_job)
                yield()
            end

            if check_stop_condition(network_settings, start_time, networks, network)
                stop_condition = true
                break
            end
        end

        if stop_condition ||
           check_stop_condition(reaction_settings, start_time, reactions, reaction)
            break
        end
    end

    if !isnothing(pbar)
        removejob!(pbar, reaction_job)
        removejob!(pbar, network_job)
    end

    if haskey(network_settings.options, :unique_candidates) &&
       network_settings.options[:unique_candidates]
        networks = unique(networks)
    end

    return networks, reactions
end

function synthesize_networks(
        problem::ProblemDefinition,
        molecule_settings::SynthesizerSettings,
        reaction_settings::SynthesizerSettings,
        network_settings::SynthesizerSettings;
        initial_molecules_count::Int = 10,
        initial_reactions_count::Int = 1000,
        pbar = nothing,
        fragment_rules::Dict{Int, Set{Expr}} = Dict{Int, Set{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[]
)
    start_time = time()

    molecules = Set{Molecule}()
    metric_name = get(network_settings.options, :similarity_metric, :none)
    combine_method = get(network_settings.options, :similarity_combine, :multiplicative)

    for molecule in problem.known_molecules
        push!(molecules, molecule)
    end
    for molecule in get_molecules(problem)
        push!(molecules, molecule)
    end
    check_stop_condition(
        molecule_settings, start_time, molecules, nothing; check_all_candidates = true
    )

    molecule_grammar = SMILES_grammar(get_atoms(problem); settings = molecule_settings, fragment_rules = fragment_rules, starting_fragments = starting_fragments)
    molecule_iterator = get_iterator(molecule_settings, molecule_grammar, :molecule)

    reactions = Set{Reaction}()

    networks = Vector{ReactionNetwork}()

    stop_condition = false
    for molecule_program in molecule_iterator
        molecule = interpret_molecule(molecule_program, molecule_grammar)
        push!(molecules, molecule)

        if length(molecules) < initial_molecules_count &&
           !check_stop_condition(molecule_settings, start_time, molecules, molecule)
            continue
        end

        reactions_grammar = reaction_grammar(
            sort_molecules_by_similarity(collect(molecules), problem.known_molecules; metric_name=metric_name); settings = reaction_settings)
        reaction_iterator = get_iterator(reaction_settings, reactions_grammar, :reaction)
        for reaction_program in reaction_iterator
            reaction = interpret_reaction(reaction_program, reactions_grammar)
            if reaction in reactions
                continue
            end
            push!(reactions, reaction)

            # molecules = Set(synthesize_molecules(get_atoms(problem), molecule_settings))
            # for molecule in get_molecules(problem)
            #     push!(molecules, molecule)
            # end
            # reactions = synthesize_reactions(unique(molecules), reaction_settings)

            if length(reactions) < initial_reactions_count &&
               !check_stop_condition(reaction_settings, start_time, reactions, reaction)
                # if length(reactions) < initial_reactions_count
                continue
            end

            # println("Found $(length(reactions)) reactions so far.")

            networks_grammar = network_grammar(
                sort_reactions_by_similarity(collect(reactions), problem.known_molecules; metric_name=metric_name, combine=combine_method); problem = problem, settings = network_settings
            )
            network_iterator = get_iterator(network_settings, networks_grammar, :network)
            for network_program in network_iterator
                network = interpret_network(network_program, networks_grammar)

                push!(networks, network)

                if check_stop_condition(network_settings, start_time, networks, network)
                    stop_condition = true
                    break
                end
            end

            if stop_condition ||
               check_stop_condition(reaction_settings, start_time, reactions, reaction)
                break
            end
        end

        if stop_condition ||
           check_stop_condition(molecule_settings, start_time, molecules, molecule)
            break
        end
    end

    if haskey(network_settings.options, :unique_candidates) &&
       network_settings.options[:unique_candidates]
        networks = unique(networks)
    end

    # TEMP
    return networks, reactions, molecules
    # return networks
end

# Similarity and sorting functions moved to similarity.jl
