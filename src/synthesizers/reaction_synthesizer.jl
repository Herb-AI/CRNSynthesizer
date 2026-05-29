
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

function get_reaction_molecules(reaction::Reaction)::Vector{Molecule}
    mols = Molecule[]
    for (count, mol) in reaction.inputs
        push!(mols, mol)
    end
    for (count, mol) in reaction.outputs
        push!(mols, mol)
    end
    return unique(mols)
end

function synthesize_reactions(
        molecules::Vector{Molecule}, settings::SynthesizerSettings;
        known_molecules::Vector{Molecule} = Molecule[],
        partial_reaction::Union{Nothing, Reaction} = nothing
)::Vector{Reaction}
    metric_name = get(settings.options, :similarity_metric, :none)

    # If a partial reaction is provided, do not add its known molecules to the set of possible molecules
    pool_molecules = molecules
    if !isnothing(partial_reaction)
        partial_mols = get_reaction_molecules(partial_reaction)
        pool_molecules = filter(m -> !(m in partial_mols), pool_molecules)
    end

    sorted_molecules = (isempty(known_molecules) || metric_name == :none) ? pool_molecules :
                       sort_molecules_by_similarity(pool_molecules, known_molecules; metric_name = metric_name)
    grammar = reaction_grammar(sorted_molecules; settings = settings, partial_reaction = partial_reaction)
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
        fragment_rules::OrderedDict{Int, Vector{Expr}} = OrderedDict{Int, Vector{Expr}}(),
        starting_fragments::Vector{Expr} = Expr[]
)
    start_time = time()

    molecules = OrderedSet{Molecule}()
    metric_name = get(reaction_settings.options, :similarity_metric, :none)

    # Maximum number of molecules to pass to the reaction grammar at once
    max_pool_size = get(reaction_settings.options, :max_reaction_pool_size, 250)

    for molecule in problem.known_molecules
        push!(molecules, molecule)
    end
    check_stop_condition(
        molecule_settings, start_time, molecules, nothing; check_all_candidates = true
    )

    molecule_grammar = SMILES_grammar(problem.atom_valences; settings = molecule_settings,
        fragment_rules = fragment_rules, starting_fragments = starting_fragments)
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

        pool_molecules = sort_molecules_by_similarity(
            collect(molecules), problem.known_molecules;
            metric_name = metric_name, cache = molecule_score_cache)
        if !isnothing(problem.partial_reaction)
            partial_mols = get_reaction_molecules(problem.partial_reaction)
            pool_molecules = filter(m -> !(m in partial_mols), pool_molecules)
        end

        # Cap the pool to prevent a combinatorial explosion in get_possibilities in balanced_reaction
        if length(pool_molecules) > max_pool_size
            pool_molecules = pool_molecules[1:max_pool_size]
        end

        reactions_grammar = reaction_grammar(
            pool_molecules; settings = reaction_settings, partial_reaction = problem.partial_reaction
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
