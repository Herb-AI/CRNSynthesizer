using CRNSynthesizer

const PROBLEMS = [
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
    )
]

function run_hardcoded_benchmarks()
    println()
    println("=======================================================")
    println("Running Hardcoded Problems Feasibility Benchmark (with oracle checkpoints)")
    println("=======================================================")

    for (name, problem) in PROBLEMS
        println("\n-------------------------------------------------------")
        println("\033[1mBenchmarking problem: $name\033[0m")

        is_feas, reason, _, _, _ = is_feasible_problem(problem)
        if is_feas
            println("  ✓ Feasibility Check: \033[32mFEASIBLE\033[0m ($reason)")
        else
            println("  ✗ Feasibility Check: \033[31mINFEASIBLE\033[0m ($reason)")
        end

        for metric in [:none, :simpson, :tanimoto, :both]
            println("\n    \033[1mSimilarity Metric: $metric\033[0m")

            fragment_flags = startswith(name, "Estherification") ? [true, false] : [true]

            for use_fragments in fragment_flags
                frag_str = use_fragments ? "With fragments" : "Without fragments"
                println("    \033[1m[$frag_str]\033[0m")
                elapsed = @elapsed res = run_problem_synthesis(
                    problem; metric = metric,
                    max_stage = :networks, use_fragments = use_fragments
                )
                success = res.success
                if success
                    println("    \033[32m✓ Target reaction network successfully synthesized in $(round(elapsed; digits=2))s!\033[0m")
                else
                    println("    \033[31m✗ Synthesis failed or timed out in $(round(elapsed; digits=2))s.\033[0m")
                end
            end
        end
    end
    println("=======================================================")
end
