using CRNSynthesizer

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")

max_time = 600

PROBLEMS = [
    (
        name = "Water problem with O2 missing",
        problem = water_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        molecules_goal = get_molecules(
            water_problem(; selected_known_indices = [2], selected_expected_indices = [2])
        ),
        reactions_goal = get_reactions(water_network()),
        network_goal = water_network()
    ),
    (
        name = "Methane Combustion problem with O2 and CO2 missing",
        problem = methane_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        molecules_goal = get_molecules(
            methane_problem(;
            selected_known_indices = [2, 4], selected_expected_indices = [2, 4]
        ),
        ),
        reactions_goal = get_reactions(methane_network()),
        network_goal = methane_network()
    ),
    (
        name = "Ethylene problem with C₂H₄O missing",
        problem = ethylene_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        molecules_goal = get_molecules(
            ethylene_problem(;
            selected_known_indices = [2], selected_expected_indices = [2])
        ),
        reactions_goal = get_reactions(ethylene_network()),
        network_goal = ethylene_network()
    ),
    (
        name = "Estherification problem with H2O, CH2O2 and CH4O missing",
        problem = estherification_problem(;
            selected_known_indices = [2, 3, 6], selected_expected_indices = [2, 3, 6]
        ),
        molecules_goal = get_molecules(
            estherification_problem(;
            selected_known_indices = [1, 4, 5], selected_expected_indices = [1, 4, 5]
        ),
        ),
        reactions_goal = get_reactions(estherification_network()),
        network_goal = estherification_network()
    )
]

for (name, problem, molecules_goal, reactions_goal, network_goal) in PROBLEMS
    # @profview for (name, problem, molecules_goal, reactions_goal, network_goal) in PROBLEMS

    println()
    println("-------------------------------------------------------")
    println("\033[1mBenchmarking problem: $name\033[0m")

    # -------------------------------------------------------
    # ----------------- Until Molecules Found ---------------
    # -------------------------------------------------------

    # Pipeline: Atoms -> Molecules
    settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound
    )
    elapsed_time = @elapsed molecules = synthesize_molecules(
        get_atoms(network_goal), settings
    )
    println(
        "[Atoms → Molecules] Found $(length(molecules)) molecules in $(elapsed_time) seconds.",
    )
    if issubset(molecules_goal, molecules)
        println("\033[32m  All goal molecules found.\033[0m")
    else
        println(
            "\033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(molecules)))\033[0m",
        )
    end

    fragment_rules = Dict{Int, Set{Expr}}()
    starting_fragments = Expr[]
    for m in molecules_goal
        e, s = parse_molecule_to_fragment_rules(m.canonical_smiles)
        for (k, v) in e
            if haskey(fragment_rules, k)
                append!(fragment_rules[k], v)
            else
                fragment_rules[k] = v
            end
        end
        append!(starting_fragments, s)
    end
    starting_fragments = unique(starting_fragments)

    for use_fragments in [false, true]
        println("\n  \033[1mUse Fragments: $use_fragments\033[0m")
        
        fr = use_fragments ? fragment_rules : Dict{Int, Set{Expr}}()
        sf = use_fragments ? starting_fragments : Expr[]

        for combine_method in [:multiplicative, :additive, :pooled, :sum]
            for metric in [:none, :simpson, :tanimoto, :tango]
                if metric == :none && combine_method != :multiplicative
                    continue # Only run :none once
                end
                println()
                println("    \033[1mSimilarity Metric: $metric, Combine: $combine_method\033[0m")

                # -------------------------------------------------------
                # ------------- Until Reactions Found -------------------
                # -------------------------------------------------------

                # Pipeline: Molecules -> Reactions
                reaction_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                molecules_for_reactions = unique(vcat(get_molecules(problem), molecules))
                elapsed_time = @elapsed candidates = synthesize_reactions(
                    unique(molecules_for_reactions), reaction_settings; known_molecules = problem.known_molecules
                )
                println(
                    "    [Molecules → Reactions] Found $(length(candidates)) reactions in $(elapsed_time) seconds.",
                )
                if issubset(reactions_goal, candidates)
                    println("    \033[32m  All goal reactions found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal reactions: $(setdiff(reactions_goal, candidates))\033[0m",
                    )
                end

                # -------------------------------------------------------
                # ----------------- Until Network Found -----------------
                # -------------------------------------------------------

                # Pipeline: Problem -> Molecules -> Reactions -> Networks
                molecule_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                reaction_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                network_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 10, goal = [network_goal], benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                elapsed_time = @elapsed (networks,
                    reactions,
                    found_molecules) = synthesize_networks(
                    problem,
                    molecule_settings,
                    reaction_settings,
                    network_settings;
                    initial_molecules_count = 1000,
                    initial_reactions_count = 10000,
                    fragment_rules = fr,
                    starting_fragments = sf
                )
                println(
                    "    [Problem → Molecules → Reactions → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
                )
                if issubset([network_goal], networks)
                    println("    \033[32m  All goal networks found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
                    )
                end
                if issubset(reactions_goal, reactions)
                    println("    \033[32m  All goal reactions found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal reactions: $(setdiff(reactions_goal, Set(reactions)))\033[0m",
                    )
                end
                if issubset(molecules_goal, found_molecules)
                    println("    \033[32m  All goal molecules found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(found_molecules)))\033[0m",
                    )
                end
            end
        end
    end
end
