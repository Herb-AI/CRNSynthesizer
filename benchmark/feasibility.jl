using CRNSynthesizer

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")
include("data/example_SYN.jl")

max_time = 120

PROBLEMS = [
    (
        name = "Water problem with O2 missing",
        problem = water_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        )
    ),
    (
        name = "Methane Combustion problem with O2 and CO2 missing",
        problem = methane_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        )
    ),
    (
        name = "Ethylene problem with C₂H₄O missing",
        problem = ethylene_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        )
    ),
    (
        name = "Estherification problem with H2O, CH2O2 and CH4O missing",
        problem = estherification_problem(;
            selected_known_indices = [2, 3, 6], selected_expected_indices = [2, 3, 6]
        )
    ),#=
    (
        name = "Example missing molecule problem from SynRXN",
        problem = syn_problem(),
    )=#
]

for (name, problem) in PROBLEMS
    # @profview for (name, problem) in PROBLEMS

    println()
    println("-------------------------------------------------------")
    println("\033[1mBenchmarking problem: $name\033[0m")

    # -------------------------------------------------------
    # ----------------- Until Molecules Found ---------------
    # -------------------------------------------------------

    fragment_rules = Dict{Int, Set{Expr}}()
    starting_fragments = Expr[]
    for m in problem.known_molecules
        e, s = parse_molecule_to_fragment_rules(m.canonical_smiles)
        for (k, v) in e
            if haskey(fragment_rules, k)
                union!(fragment_rules[k], v)
            else
                fragment_rules[k] = v
            end
        end
        append!(starting_fragments, s)
    end
    starting_fragments = unique(starting_fragments)

    for use_fragments in [true, false]
        println("\n  \033[1mUse Fragments: $use_fragments\033[0m")

        fr = use_fragments ? fragment_rules : Dict{Int, Set{Expr}}()
        sf = use_fragments ? starting_fragments : Expr[]

        println("Starting fragments: ", sf)
        println("Fragment rules: ", fr)

        # Pipeline: Atoms -> Molecules
        settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 10, goal = problem.goal_molecules, benchmark_type = UntilFound
        )
        elapsed_time = @elapsed molecules = synthesize_molecules(
            problem.atom_valences; settings, starting_fragments = sf, fragment_rules = fr
        )
        println(
            "[Atoms → Molecules] Found $(length(molecules)) molecules in $(elapsed_time) seconds.",
        )
        if issubset(problem.goal_molecules, molecules)
            println("\033[32m  All goal molecules found.\033[0m")
        else

            println(
                "\033[31m  Missing goal molecules: $(setdiff(problem.goal_molecules, Set(molecules)))\033[0m",
            )
        end

        for combine_method in [:multiplicative]#, :additive, :pooled, :sum]
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
                    max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network), benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                molecules_for_reactions = unique(vcat(problem.known_molecules, molecules))
                elapsed_time = @elapsed candidates = synthesize_reactions(
                    molecules_for_reactions, reaction_settings; known_molecules = problem.known_molecules
                )
                println(
                    "    [Molecules → Reactions] Found $(length(candidates)) reactions in $(elapsed_time) seconds.",
                )
                if issubset(get_reactions(problem.goal_network), candidates)
                    println("    \033[32m  All goal reactions found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal reactions: $(setdiff(get_reactions(problem.goal_network), candidates))\033[0m",
                    )
                end

                # -------------------------------------------------------
                # ----------------- Until Network Found -----------------
                # -------------------------------------------------------

                # Pipeline: Problem -> Molecules -> Reactions -> Networks
                molecule_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 5, goal = problem.goal_molecules, benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                reaction_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network), benchmark_type = UntilFound,
                    options = Dict{Symbol, Any}(:similarity_metric => metric, :similarity_combine => combine_method)
                )
                network_settings = SynthesizerSettings(;
                    max_time = max_time, max_depth = 5, goal = [problem.goal_network], benchmark_type = UntilFound,
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
                    initial_reactions_count = 100000,
                    fragment_rules = fr,
                    starting_fragments = sf
                )
                println(
                    "    [Problem → Molecules → Reactions → Networks] Found $(length(networks)) networks in $(elapsed_time) seconds.",
                )
                if issubset([problem.goal_network], networks)
                    println("    \033[32m  All goal networks found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal networks: $(setdiff([problem.goal_network], Set(networks)))\033[0m",
                    )
                end
                if issubset(get_reactions(problem.goal_network), reactions)
                    println("    \033[32m  All goal reactions found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal reactions: $(setdiff(get_reactions(problem.goal_network), Set(reactions)))\033[0m",
                    )
                end
                if issubset(problem.goal_molecules, found_molecules)
                    println("    \033[32m  All goal molecules found.\033[0m")
                else
                    println(
                        "    \033[31m  Missing goal molecules: $(setdiff(problem.goal_molecules, Set(found_molecules)))\033[0m",
                    )
                end
            end
        end
    end
end
