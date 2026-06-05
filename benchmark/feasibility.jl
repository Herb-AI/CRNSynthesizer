using CRNSynthesizer
using DataFrames
using Dates

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")
include("data/synrxn_loader.jl")
import .SynRXNLoader
include("syn_parsers.jl")
include("log_utils.jl")
include("synthesis_runner.jl")
include("hardcoded_benchmarks.jl")

# -------------------------------------------------------------
# Automated SynRXN rbl Benchmark
# -------------------------------------------------------------
function run_automated_rbl_benchmark(;
        dataset::String = "mos", max_scan::Int = 100, max_synthesis_runs::Int = 5)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    run_dir = joinpath(@__DIR__, "results", timestamp)
    mkpath(run_dir)

    log_file = open(joinpath(run_dir, "console.log"), "w")
    old_stdout = stdout
    old_stderr = stderr

    out_rd, out_wr = redirect_stdout()
    err_rd, err_wr = redirect_stderr()

    out_task = @async begin
        while !eof(out_rd)
            line = readline(out_rd; keep = true)
            write(old_stdout, line)
            write(log_file, line)
            flush(log_file)
        end
    end

    err_task = @async begin
        while !eof(err_rd)
            line = readline(err_rd; keep = true)
            write(old_stderr, line)
            write(log_file, line)
            flush(log_file)
        end
    end

    try
        println()
        println("=======================================================")
        println("Running Automated SynRXN rbl/$dataset Tractability Benchmark")
        println("=======================================================")

        df = SynRXNLoader.load_synrxn_dataset("rbl", dataset)
        total_records = nrow(df)
        scan_limit = min(total_records, max_scan)

        println("Successfully loaded SynRXN rbl/$dataset dataset with $total_records records.")
        println("Scanning the first $scan_limit records for parsability...")

        parsed_problems = Tuple{String, ProblemDefinition}[]
        unparsed_reasons = Dict{String, Int}()

        parsed_count = 0
        untractable_count = 0

        for idx in 1:scan_limit
            row = df[idx, :]
            r_id = row.r_id
            rxn_str = row.rxn
            gt_str = row.ground_truth

            if ismissing(rxn_str) || ismissing(gt_str)
                untractable_count += 1
                println("  \033[31m✗ Missing reaction data for record $r_id\033[0m")
                reason = "Missing reaction data"
                unparsed_reasons[reason] = get(unparsed_reasons, reason, 0) + 1
                continue
            end

            reaction_str = rxn_str * "," * gt_str
            problem = nothing
            try
                problem = parse_syn_problem(reaction_str)
            catch e
                untractable_count += 1
                println("  \033[31m✗ Failed to parse reaction for record $r_id: ", e, "\033[0m")
                reason = "Failed to parse reaction"
                unparsed_reasons[reason] = get(unparsed_reasons, reason, 0) + 1
                continue
            end
            parsed_count += 1
            push!(parsed_problems, (r_id, problem))
        end

        println("\nClassification Summary (first $scan_limit records):")
        println("  ✓ Parsed:   $parsed_count")
        println("  ✗ Failed parsing: $untractable_count")
        for (
            reason, count) in sort(collect(unparsed_reasons); by = x -> x[2], rev = true)
            println("    - $reason: $count")
        end

        # Data structures to accumulate detailed evaluation results
        results_df = DataFrame(
            r_id = String[],
            similarity_metric = String[],
            use_fragments = Bool[],
            success = Bool[],
            elapsed_time_seconds = Float64[],
            molecules_synthesized = Int[],
            reactions_synthesized = Int[]
        )

        # Dictionary to store missing molecule synthesis stats per use_fragments setting
        missing_molecule_synthesis_stats = Dict{Bool, Dict{Symbol, Any}}()

        # Run synthesis evaluation on a subset of parsed problems
        synthesis_eval_limit = min(length(parsed_problems), max_synthesis_runs)
        if synthesis_eval_limit > 0
            println("\nPerforming JIT warm-up on the first parsed problem...")
            # Warm-up run to compile necessary functions
            try
                run_problem_synthesis(
                    parsed_problems[1][2]; metric = :none,
                    max_stage = :molecules, use_fragments = true
                )
                run_problem_synthesis(
                    parsed_problems[1][2]; metric = :both,
                    max_stage = :reactions, use_fragments = true
                )
            catch e
                println("  \033[33m! Warm-up encountered an error (ignoring): ", e, "\033[0m")
            end
            println("Warm-up complete.\n")

            println("\nEvaluating synthesis (stages molecules -> reactions only) on the first $synthesis_eval_limit parsed problems...")

            for use_fragments in [true, false]
                frag_str = use_fragments ? "With fragments" : "Without fragments"
                println("\n    \033[1m[$frag_str]\033[0m")

                missing_molecule_synthesis_successes = Bool[]
                missing_molecule_synthesis_runtimes = Float64[]
                missing_molecule_synthesis_untractable_sizes = Int[]
                missing_molecule_synthesis_molecules_synthesized = Int[]

                println("    Running missing molecule synthesis subproblem (Atoms → Molecules)...")
                filter_results = Dict{String, Tuple{Bool, Float64, Vector{Int}, Int}}()
                for (i,
                    (r_id, problem)) in enumerate(parsed_problems[1:synthesis_eval_limit])
                    println("      Filtering reaction $r_id ($i/$synthesis_eval_limit)...")
                    elapsed = @elapsed res = run_problem_synthesis(
                        problem; metric = :none,
                        max_stage = :molecules, use_fragments = use_fragments
                    )
                    success = res.success
                    filter_results[r_id] = (
                        success, elapsed, res.missing_goal_sizes, res.molecules_count)

                    push!(missing_molecule_synthesis_successes, success)
                    push!(missing_molecule_synthesis_runtimes, elapsed)
                    push!(missing_molecule_synthesis_molecules_synthesized, res.molecules_count)
                    if !success
                        append!(missing_molecule_synthesis_untractable_sizes, res.missing_goal_sizes)
                    end

                    if success
                        println("        \033[32m✓ Missing molecule synthesis subproblem passed in $(round(elapsed; digits=2))s\033[0m")
                    else
                        println("        \033[31m✗ Missing molecule synthesis subproblem failed in $(round(elapsed; digits=2))s\033[0m")
                    end
                end

                missing_molecule_synthesis_stats[use_fragments] = Dict{Symbol, Any}(
                    :successes => missing_molecule_synthesis_successes,
                    :runtimes => missing_molecule_synthesis_runtimes,
                    :unfeasible_sizes => missing_molecule_synthesis_untractable_sizes,
                    :molecules_synthesized => missing_molecule_synthesis_molecules_synthesized
                )

                for metric in [:none, :simpson, :tanimoto, :both]
                    println("\n  \033[1mSimilarity Metric: $metric\033[0m")
                    successful_runs = 0
                    total_time = 0.0

                    for (i,
                        (r_id, problem)) in
                        enumerate(parsed_problems[1:synthesis_eval_limit])
                        success_filter, filter_elapsed, missing_sizes,
                        mol_count = filter_results[r_id]

                        if !success_filter
                            println("      Evaluating reaction $r_id ($i/$synthesis_eval_limit)...")
                            println("        \033[31m✗ Skipped reaction synthesis (missing molecule synthesis subproblem failed).\033[0m")

                            push!(results_df,
                                (
                                    r_id = r_id,
                                    similarity_metric = string(metric),
                                    use_fragments = use_fragments,
                                    success = false,
                                    elapsed_time_seconds = filter_elapsed,
                                    molecules_synthesized = mol_count,
                                    reactions_synthesized = 0
                                ))
                            continue
                        end

                        println("      Evaluating reaction $r_id ($i/$synthesis_eval_limit)...")
                        elapsed = @elapsed res = run_problem_synthesis(
                            problem; metric = metric,
                            max_stage = :reactions, use_fragments = use_fragments
                        )
                        success = res.success

                        if success
                            println("        \033[32m✓ Target reaction successfully synthesized in $(round(elapsed; digits=2))s!\033[0m")
                            successful_runs += 1
                            total_time += elapsed
                        else
                            println("        \033[31m✗ Synthesis failed or timed out in $(round(elapsed; digits=2))s.\033[0m")
                        end

                        # Record detailed run data
                        push!(results_df,
                            (
                                r_id = r_id,
                                similarity_metric = string(metric),
                                use_fragments = use_fragments,
                                success = success,
                                elapsed_time_seconds = elapsed,
                                molecules_synthesized = res.molecules_count,
                                reactions_synthesized = res.reactions_count
                            ))
                    end

                    success_rate = (successful_runs / synthesis_eval_limit) * 100
                    avg_time = successful_runs > 0 ? total_time / successful_runs : 0.0

                    println("\n      Synthesis Performance (Metric: $metric, $frag_str) on Tractable Subspace:")
                    println("        - Evaluated runs: $synthesis_eval_limit")
                    println("        - Success rate:   $(round(success_rate; digits=1))%")
                    println("        - Average time:   $(round(avg_time; digits=2))s")
                end
            end

            save_benchmark_results(run_dir, results_df, dataset, missing_molecule_synthesis_stats)
        else
            println("\nNo tractable problems found to evaluate.")
        end

        println("=======================================================")
    finally
        redirect_stdout(old_stdout)
        redirect_stderr(old_stderr)
        close(out_wr)
        close(err_wr)
        wait(out_task)
        wait(err_task)
        close(log_file)
    end
end

run_hardcoded_benchmarks()

#run_automated_rbl_benchmark(; dataset = "complex", max_scan = 20, max_synthesis_runs = 10)
#run_automated_rbl_benchmark(; dataset = "mbs", max_scan = 100, max_synthesis_runs = 10)
#run_automated_rbl_benchmark(; dataset = "mnc", max_scan = 100, max_synthesis_runs = 10)
#run_automated_rbl_benchmark(; dataset = "mos", max_scan = 100, max_synthesis_runs = 10)
