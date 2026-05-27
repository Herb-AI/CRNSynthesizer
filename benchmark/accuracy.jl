using CRNSynthesizer
using Term.Progress
using Term: with, update!

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")

function safe_score(network, problem)
    try
        return score(network, problem)
    catch
        return Inf
    end
end

function rank_networks(pbar, networks, problem, network_goal)
    if !(network_goal in networks)
        println("\033[31mTarget network not found in the synthesized networks.\033[0m")
        # println(networks)
        return nothing
    end

    job = addjob!(pbar; N = length(networks))
    start!(pbar)
    network_scores = Dict()
    for network in networks
        network_scores[network] = safe_score(network, problem)
        update!(job)
        sleep(0.001)
        render(pbar)
    end
    stop!(pbar)

    sorted_network_scores = sort(collect(network_scores); by = x -> x[2], rev = false)
    println(
        "\033[32mTarget network is in position ",
        findfirst(x -> x[1] == network_goal, sorted_network_scores),
        " out of ",
        length(networks),
        " networks synthesized.\033[0m"
    )
end

max_time = 600

const PROBLEMS = [
    (
        name = "Water problem with O2 missing",
        problem = water_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        num_molecules = 5,
        num_reactions = 10,
        num_networks = 1000
    ),
    (
        name = "Methane Combustion problem with O2 and CO2 missing",
        problem = methane_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        num_molecules = 10,
        num_reactions = 100,
        num_networks = 1000
    ),
    (
        name = "Ethylene problem with C₂H₄O missing",
        problem = ethylene_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ),
        num_molecules = 300,
        num_reactions = 10000,
        num_networks = 1000
    ),
    (
        name = "Estherification problem with H2O, CH2O2 and CH4O missing",
        problem = estherification_problem(;
            selected_known_indices = [2, 3, 6], selected_expected_indices = [2, 3, 6]
        ),
        num_molecules = 70,
        num_reactions = 56000,
        num_networks = 3000
    )
]

# @profview for (name, problem, num_molecules, num_reactions, num_networks) in PROBLEMS
for (name, problem, num_molecules, num_reactions, num_networks) in PROBLEMS
    println()
    println("-------------------------------------------------------")
    println("\033[1mBenchmarking problem: $name\033[0m")

    network_goal = problem.goal_network

    # Setup
    pbar = ProgressBar()
    molecule_settings = SynthesizerSettings(; max_time = max_time, max_depth = 10)
    reaction_settings = SynthesizerSettings(; max_time = max_time, max_depth = 10)
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_programs = num_networks, max_depth = 10
    )

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

    #= DNF for every problem Pipelines without molecule step are not feasible based on Wijers conclusions
    # Pipeline: Problem -> Networks
    elapsed_time = @elapsed networks = synthesize_networks(problem, network_settings)
    println(
        "[Problem → Networks] Found ",
        length(networks),
        " networks in ",
        elapsed_time,
        " seconds."
    )
    rank_networks(pbar, networks, problem, network_goal) =#

    # Pipeline: Problem -> Molecules -> Networks
    elapsed_time = @elapsed (networks,
        molecules) = synthesize_networks(
        problem, molecule_settings, network_settings;
        initial_molecules_count = num_molecules,
        fragment_rules = fragment_rules,
        starting_fragments = starting_fragments
    )
    println(
        "[Problem → Molecules → Networks] Found ",
        length(networks),
        " networks in ",
        elapsed_time,
        " seconds."
    )
    rank_networks(pbar, networks, problem, network_goal)#= DNF for every problem Pipelines without molecule step are not feasible based on Wijers conclusions =#

    #= Pipeline: Problem -> Reactions -> Networks
    elapsed_time = @elapsed (networks,
        reactions) = synthesize_networks_2(
        problem, reaction_settings, network_settings; initial_reactions_count = num_reactions
    )
    println(
        "[Problem → Reactions → Networks] Found ",
        length(networks),
        " networks in ",
        elapsed_time,
        " seconds."
    )
    rank_networks(pbar, networks, problem, network_goal) =#

    # Pipeline: Problem -> Molecules -> Reactions -> Networks
    elapsed_time = @elapsed (networks,
        reactions,
        molecules) = synthesize_networks(
        problem,
        molecule_settings,
        reaction_settings,
        network_settings;
        initial_molecules_count = num_molecules,
        initial_reactions_count = num_reactions,
        fragment_rules = fragment_rules,
        starting_fragments = starting_fragments
    )
    println(
        "[Problem → Molecules → Reactions → Networks] Found ",
        length(networks),
        " networks in ",
        elapsed_time,
        " seconds."
    )
    rank_networks(pbar, networks, problem, network_goal)
end
