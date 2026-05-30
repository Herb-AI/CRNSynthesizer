using CRNSynthesizer
using DataStructures

"""
    run_problem_synthesis(problem::ProblemDefinition; kwargs...) -> Bool

Synthesize the molecular network for a problem and report success.
"""
function run_problem_synthesis(
        problem::ProblemDefinition;
        metric::Symbol = :none,
        combine_method::Symbol = :multiplicative,
        max_stage::Symbol = :networks,
        use_fragments::Bool = true
)
    fragment_rules = OrderedDict{Int, Vector{Expr}}()
    starting_fragments = Expr[]

    if use_fragments
        for m in problem.known_molecules
            e, s = parse_molecule_to_fragment_rules(m.canonical_smiles)
            for (k, v) in e
                if haskey(fragment_rules, k)
                    append!(fragment_rules[k], v)
                    unique!(fragment_rules[k])
                else
                    fragment_rules[k] = v
                end
            end
            append!(starting_fragments, s)
        end
        starting_fragments = unique(starting_fragments)
    end

    # ---------------------------------------------------------
    # Atoms -> Molecules
    # ---------------------------------------------------------
    if max_stage != :reactions
        molecule_settings = SynthesizerSettings(; max_programs = 250, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
            options = Dict{Symbol, Any}(
                :unique_candidates => true, :similarity_metric => metric,
                :similarity_combine => combine_method)
        )

        molecules = Vector{Molecule}()
        try
            molecules = synthesize_molecules(
                problem.atom_valences; settings = molecule_settings, starting_fragments = starting_fragments,
                fragment_rules = fragment_rules
            )
        catch e
            println("Molecule synthesis failed with error: ", e)
            missing_goal_sizes = [length(m.atoms) for m in problem.goal_molecules]
            return (success = false, molecules_count = 0, reactions_count = 0, missing_goal_sizes = missing_goal_sizes)
        end

        println("[Atoms → Molecules] Found $(length(molecules)) molecules.")
        if issubset(problem.goal_molecules, molecules)
            println("      \033[32m✓ All goal molecules found.\033[0m")
        else
            println("      \033[31m✗ Missing goal molecules.\033[0m")
            missing_mols = setdiff(problem.goal_molecules, molecules)
            missing_goal_sizes = [length(m.atoms) for m in missing_mols]
            return (success = false, molecules_count = length(molecules), reactions_count = 0, missing_goal_sizes = missing_goal_sizes)
        end
        if max_stage == :molecules
            return (success = true, molecules_count = length(molecules), reactions_count = 0, missing_goal_sizes = Int[])
        end
    end

    # ---------------------------------------------------------
    # Atoms -> Molecules -> Reactions
    # ---------------------------------------------------------
    if max_stage == :reactions
        molecule_settings = SynthesizerSettings(; max_programs = 250, max_depth = 9,
            options = Dict{Symbol, Any}(
                :unique_candidates => true, :similarity_metric => metric,
                :similarity_combine => combine_method)
        )
        reaction_settings = SynthesizerSettings(; max_programs = 30000, max_depth = 5, goal = get_reactions(problem.goal_network),
            benchmark_type = UntilFound,
            options = Dict{Symbol, Any}(
                :unique_candidates => true, :similarity_metric => metric,
                :similarity_combine => combine_method)
        )

        candidates = Vector{CRNSynthesizer.Reaction}()
        found_molecules = OrderedSet{Molecule}()
        try
            candidates,
            found_molecules = synthesize_reactions(
                problem,
                molecule_settings,
                reaction_settings;
                initial_molecules_count = 250,
                fragment_rules = fragment_rules,
                starting_fragments = starting_fragments
            )
        catch e
            println("[Problem → Reactions] Reaction synthesis failed with error: ", e)
            return (success = false, molecules_count = 0, reactions_count = 0, missing_goal_sizes = Int[])
        end

        println("[Problem → Reactions] Found $(length(candidates)) reactions, $(length(found_molecules)) molecules.")
        if issubset(problem.goal_molecules, found_molecules)
            println("      \033[32m✓ All goal molecules found.\033[0m")
        else
            println("      \033[31m✗ Missing goal molecules.\033[0m")
            missing_mols = setdiff(problem.goal_molecules, found_molecules)
            missing_goal_sizes = [length(m.atoms) for m in missing_mols]
            return (success = false, molecules_count = length(found_molecules), reactions_count = length(candidates), missing_goal_sizes = missing_goal_sizes)
        end
        if issubset(get_reactions(problem.goal_network), candidates)
            println("      \033[32m✓ All goal reactions found.\033[0m")
        else
            println("      \033[31m✗ Missing goal reactions.\033[0m")
            return (success = false, molecules_count = length(found_molecules), reactions_count = length(candidates), missing_goal_sizes = Int[])
        end
        return (success = true, molecules_count = length(found_molecules), reactions_count = length(candidates), missing_goal_sizes = Int[])
    else
        # ---------------------------------------------------------
        # Molecules -> Reactions
        # ---------------------------------------------------------
        reaction_settings = SynthesizerSettings(; max_programs = 30000, max_depth = 5, goal = get_reactions(problem.goal_network),
            benchmark_type = UntilFound,
            options = Dict{Symbol, Any}(
                :unique_candidates => true, :similarity_metric => metric,
                :similarity_combine => combine_method)
        )
        molecules_for_reactions = unique(vcat(problem.known_molecules, molecules))
        candidates = Vector{CRNSynthesizer.Reaction}()
        try
            candidates = synthesize_reactions(
                molecules_for_reactions,
                reaction_settings;
                known_molecules = problem.known_molecules
            )
        catch e
            println("[Molecules → Reactions] Reaction synthesis failed with error: ", e)
            return (success = false, molecules_count = 0, reactions_count = 0, missing_goal_sizes = Int[])
        end

        println("[Molecules → Reactions] Found $(length(candidates)) reactions.")
        if issubset(get_reactions(problem.goal_network), candidates)
            println("      \033[32m✓ All goal reactions found.\033[0m")
        else
            println("      \033[31m✗ Missing goal reactions.\033[0m")
            return (success = false, molecules_count = length(molecules), reactions_count = length(candidates), missing_goal_sizes = Int[])
        end
    end

    # ---------------------------------------------------------
    # Full Pipeline (with oracle) (Atoms -> Molecules -> Reactions -> Networks)
    # ---------------------------------------------------------
    pipeline_molecule_settings = SynthesizerSettings(; max_programs = 250, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )
    pipeline_reaction_settings = SynthesizerSettings(; max_programs = 30000, max_depth = 5, goal = get_reactions(problem.goal_network),
        benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )
    network_settings = SynthesizerSettings(; max_programs = 5000, max_depth = 5, goal = [problem.goal_network], benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )

    networks = Vector{ReactionNetwork}()
    reactions_found = Vector{CRNSynthesizer.Reaction}()
    found_molecules = Vector{Molecule}()
    try
        (networks,
            reactions_found,
            found_molecules) = synthesize_networks(
            problem,
            pipeline_molecule_settings,
            pipeline_reaction_settings,
            network_settings;
            initial_molecules_count = 250,
            initial_reactions_count = 30000,
            fragment_rules = fragment_rules,
            starting_fragments = starting_fragments
        )
    catch e
        println("[Problem → Networks] Network synthesis failed with error: ", e)
        return (success = false, molecules_count = length(found_molecules), reactions_count = 0, missing_goal_sizes = Int[])
    end

    println("[Problem → Networks] Found $(length(networks)) networks, $(length(reactions_found)) reactions, $(length(found_molecules)) molecules.")
    if issubset([problem.goal_network], networks)
        println("      \033[32m✓ Target network found.\033[0m")
    else
        println("      \033[31m✗ Target network not found.\033[0m")
    end
    if issubset(get_reactions(problem.goal_network), reactions_found)
        println("      \033[32m✓ All goal reactions found in pipeline.\033[0m")
    else
        println("      \033[31m✗ Missing goal reactions in pipeline.\033[0m")
    end
    if issubset(problem.goal_molecules, found_molecules)
        println("      \033[32m✓ All goal molecules found in pipeline.\033[0m")
    else
        println("      \033[31m✗ Missing goal molecules in pipeline.\033[0m")
    end

    return (success = issubset([problem.goal_network], networks), molecules_count = length(found_molecules), reactions_count = length(reactions_found), missing_goal_sizes = Int[])
end
