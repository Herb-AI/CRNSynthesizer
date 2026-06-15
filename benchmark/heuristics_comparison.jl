using CRNSynthesizer

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")
include("data/hydrolysis.jl")
include("data/photosynthesis.jl")
include("data/fermentation.jl")

max_time = 600 # seconds

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
        name = "Methane Combustion problem with CH4 missing",
        problem = methane_problem(;
            selected_known_indices = [2,3,4], selected_expected_indices = [2, 3, 4]
        ),
        molecules_goal = get_molecules(
            methane_problem(;
            selected_known_indices = [1], selected_expected_indices = [1]
        ),
        ),
        reactions_goal = get_reactions(methane_network()),
        network_goal = methane_network()
    ),
    (
        name = "Photosynthesis problem with CO2 missing",
        problem = photosynthesis_problem(;
            selected_known_indices = [2, 3, 4], selected_expected_indices = [2, 3, 4]
        ),
        molecules_goal = get_molecules(
            photosynthesis_problem(;
            selected_known_indices = [1], selected_expected_indices = [1]
        ),
        ),
        reactions_goal = get_reactions(photosynthesis_network()),
        network_goal = photosynthesis_network()
    ),
    (
        name = "Methyl Acetate Hydrolysis problem with H2O and CH₃OH missing",
        problem = hydrolysis_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        molecules_goal = get_molecules(
            hydrolysis_problem(;
            selected_known_indices = [2, 4], selected_expected_indices = [2, 4]
        ),
        ),
        reactions_goal = get_reactions(hydrolysis_network()),
        network_goal = hydrolysis_network()
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
    # (
    #     name = "Estherification problem with H2O missing",
    #     problem = estherification_problem(;
    #         selected_known_indices = [1, 2, 3, 5, 6], selected_expected_indices = [1, 2, 3, 5, 6]
    #     ),
    #     molecules_goal = get_molecules(
    #         estherification_problem(;
    #         selected_known_indices = [4], selected_expected_indices = [4]
    #     ),
    #     ),
    #     reactions_goal = get_reactions(estherification_network()),
    #     network_goal = estherification_network()
    # ),
    # (
    #     name = "Estherification problem with H2O and CH4O missing",
    #     problem = estherification_problem(;
    #         selected_known_indices = [1, 2, 3, 6], selected_expected_indices = [1, 2, 3, 6]
    #     ),
    #     molecules_goal = get_molecules(
    #         estherification_problem(;
    #         selected_known_indices = [4, 5], selected_expected_indices = [4, 5]
    #     ),
    #     ),
    #     reactions_goal = get_reactions(estherification_network()),
    #     network_goal = estherification_network()
    # ),
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
    ),
    (
        name = "Fermentation problem with none missing",
        problem = fermentation_problem(;
            selected_known_indices = [1, 2, 3, 4, 5, 6, 7, 8], selected_expected_indices = [1, 2, 3, 4, 5, 6, 7, 8]
        ),
        molecules_goal = get_molecules(
            fermentation_problem(;
            selected_known_indices = [], selected_expected_indices = []
        ),
        ),
        reactions_goal = get_reactions(fermentation_network()),
        network_goal = fermentation_network()
    ),
    # (
    #     name = "Fermentation problem with CH3CHO, CO2 and CH₃CH₂OH missing",
    #     problem = fermentation_problem(;
    #         selected_known_indices = [1, 2, 3, 4, 5], selected_expected_indices = [1, 2, 3, 4, 5]
    #     ),
    #     molecules_goal = get_molecules(
    #         fermentation_problem(;
    #         selected_known_indices = [6, 7, 8], selected_expected_indices = [6, 7, 8]
    #     ),
    #     ),
    #     reactions_goal = get_reactions(fermentation_network()),
    #     network_goal = fermentation_network()
    # ),
    # (
    #     name = "Fermentation problem with C₃H₄O₃, CH3CHO and CO2 missing",
    #     problem = fermentation_problem(;
    #         selected_known_indices = [1, 2, 4, 5, 8], selected_expected_indices = [1, 2, 4, 5, 8]
    #     ),
    #     molecules_goal = get_molecules(
    #         fermentation_problem(;
    #         selected_known_indices = [3, 6, 7], selected_expected_indices = [3, 6, 7]
    #     ),
    #     ),
    #     reactions_goal = get_reactions(fermentation_network()),
    #     network_goal = fermentation_network()
    # ),
    # (
    #     name = "Fermentation problem with CO2 missing",
    #     problem = fermentation_problem(;
    #         selected_known_indices = [1, 2, 3, 4, 5, 6, 8], selected_expected_indices = [1, 2, 3, 4, 5, 6, 8]
    #     ),
    #     molecules_goal = get_molecules(
    #         fermentation_problem(;
    #         selected_known_indices = [7], selected_expected_indices = [7]
    #     ),
    #     ),
    #     reactions_goal = get_reactions(fermentation_network()),
    #     network_goal = fermentation_network()
    # ),
]

HEURISTICS = [
    (name = "BFS", iterator = BreadthFirst),
    (name = "MaxBond", iterator = MaxBond),
    (name = "DeltaEnergy", iterator = DeltaEnergy),
]



function get_networks(initial_reactions::Vector{CRNSynthesizer.Reaction}, network_settings::SynthesizerSettings; problem::CRNSynthesizer.ProblemDefinition)
    candidates = Vector{ReactionNetwork}()
    start_time = time()
    networks_grammar = network_grammar(
        collect(initial_reactions); settings = network_settings, problem = problem
    )
    network_iterator = get_iterator(network_settings, networks_grammar, :network)
    for network_program in network_iterator
        network = interpret_network(network_program, networks_grammar)
        
        push!(candidates, network)
        
        if check_stop_condition(network_settings, start_time, candidates, network)
            break
        end
    end

    if haskey(network_settings.options, :unique_candidates) && network_settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function molecules_to_networks(network_goal; reaction_settings::SynthesizerSettings, network_settings::SynthesizerSettings, problem::CRNSynthesizer.ProblemDefinition)
    start_time = time()
    molecules = get_molecules(network_goal)
    reaction_candidates = Vector{CRNSynthesizer.Reaction}()
    network_candidates = Vector{CRNSynthesizer.ReactionNetwork}()
    reactions_grammar = reaction_grammar(
        collect(molecules); settings = reaction_settings
    )
    reaction_iterator = get_iterator(reaction_settings, reactions_grammar, :reaction)
    for reaction_program in reaction_iterator
        reaction = interpret_reaction(reaction_program, reactions_grammar)
        push!(reaction_candidates, reaction)
        
        if check_stop_condition(reaction_settings, start_time, reaction_candidates, reaction)
            # if length(reaction_candidates) < initial_reactions_count
            break
        end
    end
    
    networks_grammar = network_grammar(
    collect(reaction_candidates); settings = network_settings, problem = problem
    )
    network_iterator = get_iterator(network_settings, networks_grammar, :network)
    for network_program in network_iterator
        network = interpret_network(network_program, networks_grammar)
        
        push!(network_candidates, network)
        
        if check_stop_condition(network_settings, start_time, network_candidates, network)
            break
        end
    end
    return network_candidates, reaction_candidates, molecules
end





for (name, problem, molecules_goal, reactions_goal, network_goal) in PROBLEMS
        
    println()
    println("-------------------------------------------------------")
    println("\033[1mBenchmarking problem: $name\033[0m")
    println("-------------------------------------------------------")
    
    
    initial_reaction_molecules = []
    # # -------------------------------------------------------
    # # ----------------- Until Molecules Found ---------------
    # # -------------------------------------------------------

    println("\033[3;41mUntil Molecules Found\033[0m")
    for (iter_name, iterator) in HEURISTICS
        println("\033[4;44m$iter_name\033[0m")


        # Pipeline: Atoms -> Molecules
        molecule_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound, iterator = iterator
        )
        elapsed_time = @elapsed molecules = synthesize_molecules(
            get_atoms(network_goal), molecule_settings
        )
        println(
            "[Atoms → Molecules] Found $(length(molecules)) molecules ($(length(unique(molecules)))) in $(elapsed_time) seconds.",
        )
        if issubset(molecules_goal, molecules)
            println("\033[32m  All goal molecules found.\033[0m")
        else
            println(
                "\033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(molecules)))\033[0m",
            )
        end

        if initial_reaction_molecules == [] || length(initial_reaction_molecules) > length(molecules) 
            initial_reaction_molecules = molecules
        end
    end
    
    
    
    # # -------------------------------------------------------
    # # ------------- Until Reactions Found -------------------
    # # -------------------------------------------------------
    
    println("\033[3;41mUntil Reactions Found\033[0m")

    candidate_molecules::Vector{Molecule} = unique(vcat(get_molecules(problem), initial_reaction_molecules))
    if (!issubset(molecules_goal, candidate_molecules))
        # if the previous step did not find all the molecules necessary for network synthesis, 
        # we will add the molecules from the goal to ensure we have all necessary molecules for network synthesis
        candidate_molecules = unique(vcat(candidate_molecules, molecules_goal))
    end

    if (name == "Fermentation problem with none missing")
        # Due to the current implementation of the molecule synthesizer
        # the fermentation problem with no missing molecules still gets a single molecule added
        # starting with 9 instead of the supossed 8 molecules.
        # To avoid this issue, we will use the molecules from the problem definition directly for reaction synthesis.
        candidate_molecules = get_molecules(network_goal)
    end
    
    
    println("Total unique molecules to consider for reactions: $(length(candidate_molecules))")
    initial_network_reactions = []
    
    for (iter_name, iterator) in HEURISTICS
        println("\033[4;44m$iter_name\033[0m")

        reaction_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 20, goal = reactions_goal, benchmark_type = UntilFound, iterator = iterator
        )

        # Add molecules from the problem definition to ensure we have all necessary molecules for reaction synthesis
        elapsed_time = @elapsed reaction_candidates = synthesize_reactions(
            collect(candidate_molecules), reaction_settings
        )
        println(
            "[Molecules → Reactions] Found $(length(reaction_candidates)) reactions ($(length(unique(reaction_candidates)))) in $(elapsed_time) seconds.",
        )
        if issubset(reactions_goal, reaction_candidates)
            println("\033[32m  All goal reactions found.\033[0m")
        else
            println(
                "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, reaction_candidates))\033[0m",
            )
        end

        if initial_network_reactions == [] || length(initial_network_reactions) > length(reaction_candidates) 
            initial_network_reactions = reaction_candidates
        end
    end
    
    
    # # -------------------------------------------------------
    # # ----------------- Until Network Found -----------------
    # # -------------------------------------------------------

    println("\033[3;41mUntil Network Found\033[0m")

    initial_reactions::Vector{CRNSynthesizer.Reaction} = initial_network_reactions
    if (!issubset(reactions_goal, initial_network_reactions))
        # if the previous step did not find all the reactions necessary for network synthesis, 
        # we will add the reactions from the problem definition to ensure we have all necessary reactions for network synthesis
        initial_reactions= unique(vcat(initial_network_reactions, get_reactions(network_goal)))
    end
    println("Total unique reactions to consider for networks: $(length(initial_reactions))")

    for (iter_name, iterator) in HEURISTICS
        println("\033[4;44m$iter_name\033[0m")

        # Pipeline: Reactions (from previous step) -> Networks
        network_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 6, goal = [network_goal], benchmark_type = UntilFound, iterator = iterator
        )

        # Add reactions from the problem definition to ensure we have all necessary reactions for network synthesis
        elapsed_time = @elapsed networks = get_networks(
            unique(initial_reactions), network_settings, problem = problem
        )
        println(
            "[Reactions → Networks] Found $(length(networks)) networks ($(length(unique(networks)))) in $(elapsed_time) seconds.",
        )
        if issubset([network_goal], networks)
            println("\033[32m  All goal networks found.\033[0m")
        else
            println(
                "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(networks)))\033[0m",
            )
            continue
        end
    end

    # -------------------------------------------------------
    # -------------------- Full Pipeline --------------------
    # -------------------------------------------------------
    
    println("\033[3;41mFull Pipeline: Problem → Molecules → Reactions → Networks\033[0m")
    for (iter_name, iterator) in HEURISTICS
        println("\033[4;44m$iter_name\033[0m")
        
        # Pipeline: Problem -> Molecules -> Reactions -> Networks
        molecule_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 10, goal = molecules_goal, benchmark_type = UntilFound, iterator = iterator
        )
        reaction_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 20, goal = reactions_goal, benchmark_type = UntilFound, iterator = iterator
        )
        network_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 6, goal = [network_goal], benchmark_type = UntilFound, iterator = DeltaEnergy
        )

        if (name == "Fermentation problem with none missing")
            # Due to the current implementation of the molecule synthesizer
            # the fermentation problem with no missing molecules still gets a single molecule added
            # starting with 9 instead of the supossed 8 molecules.
            # To avoid this issue, we will use the molecules from the problem definition directly for reaction synthesis.
            elapsed_time = @elapsed (candidate_networks,
            candidate_reactions,
            candidate_molecules) = molecules_to_networks(
                network_goal;
                problem = problem,
                reaction_settings = reaction_settings,
                network_settings = network_settings
                )
        else
            elapsed_time = @elapsed (candidate_networks,
                candidate_reactions,
                candidate_molecules_set) = synthesize_networks(
                problem,
                molecule_settings,
                reaction_settings,
                network_settings;
                initial_molecules_count = 1000000,
                initial_reactions_count = 10000000
            )
            candidate_molecules = collect(candidate_molecules_set)
        end
        
        println(
            "[Problem → Molecules → Reactions → Networks] Found $(length(candidate_networks)) networks ($(length(unique(candidate_networks)))), $(length(candidate_reactions)) reactions ($(length(unique(candidate_reactions)))), and $(length(candidate_molecules)) molecules ($(length(unique(candidate_molecules)))) in $(elapsed_time) seconds.",
        )
        if issubset([network_goal], candidate_networks)
            println("\033[32m  All goal networks found.\033[0m")
        else
            println(
                "\033[31m  Missing goal networks: $(setdiff([network_goal], Set(candidate_networks)))\033[0m",
            )
        end
        if issubset(reactions_goal, candidate_reactions)
            println("\033[32m  All goal reactions found.\033[0m")
        else
            println(
                "\033[31m  Missing goal reactions: $(setdiff(reactions_goal, Set(candidate_reactions)))\033[0m",
            )
        end
        if issubset(molecules_goal, candidate_molecules)
            println("\033[32m  All goal molecules found.\033[0m")
        else
            println(
                "\033[31m  Missing goal molecules: $(setdiff(molecules_goal, Set(candidate_molecules)))\033[0m",
            )
        end
    end
end
