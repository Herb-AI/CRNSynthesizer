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

    # -------------------------------------------------------
    # ------------- Until Reactions Found -------------------
    # -------------------------------------------------------

    # Pipelines without molecule step are not feasible based on Wijers conclusions
    #= Pipeline: Atoms -> Reactions
    settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound
    )
    elapsed_time = @elapsed candidates = synthesize_reactions(
        get_atoms(network_goal), settings
    )
    println(
        "[Atoms → Reactions] Found $(length(candidates)) reactions in $(elapsed_time) seconds.",
    )
    if issubset(reactions_goal, candidates)
        println("\033[32m  All goal reactions found.\033[0m")
    else
        println(
            "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, Set(candidates)))\033[0m",
        )
    end =#

    # Pipeline: Molecules -> Reactions
    reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound
    )
    # elapsed_time = @elapsed candidates = synthesize_reactions(get_molecules(network_goal), reaction_settings)
    molecules = unique(vcat(get_molecules(problem), molecules))
    elapsed_time = @elapsed candidates = synthesize_reactions(
        unique(molecules), reaction_settings; known_molecules = problem.known_molecules
    )
    println(
        "[Molecules → Reactions] Found $(length(candidates)) reactions in $(elapsed_time) seconds.",
    )
    if issubset(reactions_goal, candidates)
        println("\033[32m  All goal reactions found.\033[0m")
    else
        println(
            "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, candidates))\033[0m",
        )
    end

    # -------------------------------------------------------
    # ----------------- Until Network Found -----------------
    # -------------------------------------------------------

    #= DNF for every problem 
    # Pipeline: Problem -> Networks
    settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = [network_goal], benchmark_type = UntilFound
    )
    elapsed_time = @elapsed networks = synthesize_networks(problem, settings)
    println(
        "[Problem → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
    )
    if issubset([network_goal], networks)
        println("\033[32m  All goal networks found.\033[0m")
    else
        println(
            "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
        )
    end =#

    #= Pipeline: Problem -> Molecules -> Networks
    molecule_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound
    )
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = [network_goal], benchmark_type = UntilFound
    )
    elapsed_time = @elapsed (networks,
        molecules) = synthesize_networks(
        problem, molecule_settings, network_settings; initial_molecules_count = 100000
    )
    println(
        "[Problem → Molecules → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
    )
    if issubset([network_goal], networks)
        println("\033[32m  All goal networks found.\033[0m")
    else
        println(
            "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
        )
    end
    if issubset(molecules_goal, molecules)
        println("\033[32m  All goal molecules found.\033[0m")
    else
        println(
            "\033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(molecules)))\033[0m",
        )
    end =#

    #= DNF for every problem except water (still takes 500 seconds for water)
    # Pipeline: Problem -> Reactions -> Networks
    reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound
    )
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = [network_goal], benchmark_type = UntilFound
    )
    elapsed_time = @elapsed (networks,
        reactions) = synthesize_networks_2(
        problem, reaction_settings, network_settings; initial_reactions_count = 100000
    )
    println(
        "[Problem → Reactions → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
    )
    if issubset([network_goal], networks)
        println("\033[32m  All goal networks found.\033[0m")
    else
        println(
            "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
        )
    end
    if issubset(reactions_goal, reactions)
        println("\033[32m  All goal reactions found.\033[0m")
    else
        println(
            "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, Set(reactions)))\033[0m",
        )
    end =#

    # Pipeline: Problem -> Molecules -> Reactions -> Networks
    molecules_goal = unique(vcat(get_molecules(problem), molecules_goal))
    molecule_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound
    )
    reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = reactions_goal, benchmark_type = UntilFound
    )
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 10, goal = [network_goal], benchmark_type = UntilFound
    )
    elapsed_time = @elapsed (networks,
        reactions,
        molecules) = synthesize_networks(
        problem,
        molecule_settings,
        reaction_settings,
        network_settings;
        initial_molecules_count = 1000,
        initial_reactions_count = 100000
    )
    println(
        "[Problem → Molecules → Reactions → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
    )
    if issubset([network_goal], networks)
        println("\033[32m  All goal networks found.\033[0m")
    else
        println(
            "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
        )
    end
    if issubset(reactions_goal, reactions)
        println("\033[32m  All goal reactions found.\033[0m")
    else
        println(
            "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, Set(reactions)))\033[0m",
        )
    end
    if issubset(molecules_goal, molecules)
        println("\033[32m  All goal molecules found.\033[0m")
    else
        println(
            "\033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(molecules)))\033[0m",
        )
    end
    # println("num_molecules: $(length(molecules))")
    # println("num_reactions: $(length(reactions))")
end
