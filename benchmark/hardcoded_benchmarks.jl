using CRNSynthesizer
using CSV
using DataFrames
using Dates

const PROBLEMS = [
    #=(
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
    ),=#
    (
        name = "Estherification problem with H2O, CH2O2 and CH4O missing",
        problem = estherification_problem(;
            selected_known_indices = [2, 3, 6], selected_expected_indices = [2, 3, 6]
        )
    ),
    (
        name = "Ethylene problem with only C₂H₄O known",
        problem = ethylene_problem(;
            selected_known_indices = [3],
            selected_expected_indices = [3]
        )
    )
]

function run_hardcoded_benchmarks()
    println()
    println("=======================================================")
    println("Running Hardcoded Problems Tractability Benchmark")
    println("=======================================================")

    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(@__DIR__, "results", "hardcoded_" * timestamp)
    mkpath(results_dir)

    # Write metadata.txt with script commit
    metadata_path = joinpath(results_dir, "metadata.txt")
    open(metadata_path, "w") do io
        println(io, "task: hardcoded_benchmarks")
        println(io, "evaluation_script_commit: $(get_git_commit())")
    end

    ethylene_rows = []
    esterification_rows = []

    for (name, problem) in PROBLEMS
        println("\n-------------------------------------------------------")
        println("\033[1mBenchmarking problem: $name\033[0m")

        is_esterification = startswith(name, "Estherification")
        max_stage = is_esterification ? :networks : :reactions

        for metric in [:none, :simpson, :tanimoto, :both]
            println("\n    \033[1mSimilarity Metric: $metric\033[0m")
            elapsed = @elapsed res = run_problem_synthesis(
                problem; metric = metric,
                max_stage = max_stage, use_fragments = false
            )
            success = res.success
            if success
                println("    \033[32m✓ Target reaction network successfully synthesized in $(round(elapsed; digits=2))s!\033[0m")
            else
                println("    \033[31m✗ Synthesis failed or timed out in $(round(elapsed; digits=2))s.\033[0m")
            end

            if is_esterification
                # Esterification benchmark analysis
                target_rxns = get_reactions(problem.goal_network)
                rxn1_idx = findfirst(r -> r == target_rxns[1], res.reactions)
                rxn2_idx = findfirst(r -> r == target_rxns[2], res.reactions)
                num_rxns_until_target = if isnothing(rxn1_idx) || isnothing(rxn2_idx)
                    "N/A"
                else
                    max(rxn1_idx, rxn2_idx)
                end

                net_idx = findfirst(n -> n == problem.goal_network, res.networks)
                num_nets_until_target = isnothing(net_idx) ? "N/A" : net_idx

                push!(esterification_rows,
                    (
                        Metric = string(metric),
                        Reactions_Synthesised_Until_Target_Reactions = num_rxns_until_target,
                        Networks_Synthesised_Until_Target_Network = num_nets_until_target
                    ))
            else
                # Ethylene benchmark analysis
                target_rxn = get_reactions(problem.goal_network)[1]
                rxn_idx = findfirst(r -> r == target_rxn, res.reactions)
                num_rxns_before_target = isnothing(rxn_idx) ? "N/A" : rxn_idx

                push!(ethylene_rows,
                    (
                        Metric = string(metric),
                        Reactions_Synthesized_Before_Target_Count = num_rxns_before_target
                    ))
            end
        end
    end

    # Write the CSV files
    if !isempty(ethylene_rows)
        CSV.write(joinpath(results_dir, "ethylene_reactions_before_target.csv"), DataFrame(ethylene_rows))
    end
    if !isempty(esterification_rows)
        CSV.write(joinpath(results_dir, "esterification_benchmark.csv"), DataFrame(esterification_rows))
    end

    println("Saved hardcoded results & metadata.txt to: $results_dir")
    println("=======================================================")
end
