using CRNSynthesizer
using DataStructures
import TOML
using Dates
using DataFrames
using CSV
# Disable GUI plot pop-ups by running in headless mode
ENV["GKSwstype"] = "100"
using StatsPlots
using Statistics
using Serialization
# -------------------------------------------------------------
# Metadata Recovery Helpers
# -------------------------------------------------------------
function get_git_commit()
    try
        return strip(read(`git rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

function get_package_versions()
    synrxn_version = "unknown"
    try
        conda_path = joinpath(@__DIR__, "..", "CondaPkg.toml")
        if isfile(conda_path)
            conda_toml = TOML.parsefile(conda_path)
            if haskey(conda_toml, "pip") && haskey(conda_toml["pip"], "deps")
                deps = conda_toml["pip"]["deps"]
                if haskey(deps, "synrxn")
                    synrxn_version = replace(deps["synrxn"], "==" => "")
                end
            end
        end
    catch e
    end
    return "synrxn==$synrxn_version"
end


function plot_success_rate_helper(
        run_dir::String, summary_df::DataFrame, metrics::Vector{String},
        filename::String, dataset::String, num_problems::Int)

    all_groups = [
        (false, "none", "base grammar, no similarity"),
        (false, "tanimoto", "base grammar, Tanimoto + Morgan2"),
        (true, "none", "with fragments, no similarity"),
        (true, "tanimoto", "with fragments, Tanimoto + Morgan2")
    ]
    
    active_groups = filter(g -> g[2] in metrics, all_groups)
    N = length(active_groups)
    
    vals = Float64[]
    counts = Int[]
    for g in active_groups
        use_frag, metric, label = g
        row = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, summary_df)
        push!(vals, isempty(row) ? 0.0 : row[1, :success_rate])
        push!(counts, isempty(row) ? 0 : row[1, :successful_runs])
    end

    color_map = Dict(
        "base grammar, no similarity" => "#E69F00",
        "base grammar, Tanimoto + Morgan2" => "#D55E00",
        "with fragments, no similarity" => "#56B4E9",
        "with fragments, Tanimoto + Morgan2" => "#0072B2"
    )
    colors = reshape([color_map[g[3]] for g in active_groups], 1, N)
    labels = reshape([g[3] for g in active_groups], 1, N)
    
    val_matrix = reshape(vals, 1, N)
    
    annots = reshape([Plots.text("$(round(vals[i]; digits=1))%\n($(counts[i]))", 7, :center, :bottom) for i in 1:N], 1, N)
    
    p = groupedbar(
        ["$dataset\n($num_problems)"],
        val_matrix,
        xlabel = "Dataset",
        ylabel = "Success Rate (%)",
        title = "Reaction Rebalancing (SynRXN rbl)\nSuccess Rate",
        label = labels,
        legend = :outerbottom,
        ylimits = (0, 105),
        series_annotations = annots,
        color = colors,
        linewidth = 1.2,
        edgecolor = :black,
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8,
        size = (600, 450)
    )
    savefig(p, joinpath(run_dir, filename))
end

function plot_comparison_boxplot_helper(run_dir::String, succ_rxn_df::DataFrame,
        metrics::Vector{String}, col_name::Symbol, ylabel_str::String, 
        filename::String, title::String, dataset::String; is_log_scale::Bool=false, val_digits::Int=1, val_suffix::String="")

    all_groups = [
        (false, "none", "base grammar, no similarity", "#E69F00"),
        (false, "tanimoto", "base grammar, Tanimoto + Morgan2", "#D55E00"),
        (true, "none", "with fragments, no similarity", "#56B4E9"),
        (true, "tanimoto", "with fragments, Tanimoto + Morgan2", "#0072B2")
    ]
    
    active_groups = filter(g -> g[2] in metrics, all_groups)
    N = length(active_groups)
    
    p = boxplot(
        xticks = (1:1, [dataset]),
        xlims = (0.5, 1.5),
        xlabel = "Dataset",
        ylabel = ylabel_str,
        title = title,
        legend = :outerbottom,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8,
        yscale = is_log_scale ? :log10 : :identity
    )
    
    # Pre-populate legend entries to guarantee the correct legend order
    for g in active_groups
        boxplot!(p, [0], [0], label = g[3], color = g[4],
            seriestype = :shape, fillalpha = 0.7, linecolor = :black)
    end
    
    max_val = 0.0
    xs = range(1 - 0.22, 1 + 0.22, length=N)
    bar_w = 0.44 / (N - 1)
    for (idx, g) in enumerate(active_groups)
        use_frag, metric, label, col = g
        sub_df = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, succ_rxn_df)
        vals = sub_df[!, col_name]
        if !isempty(vals)
            boxplot!(p, fill(xs[idx], length(vals)), vals,
                label = "", color = col, fillalpha = 0.7,
                bar_width = bar_w, linecolor = :black)
            max_val = max(max_val, maximum(vals))
        end
    end
    
    if max_val == 0.0
        # If no data was plotted, show a message
        annotate!(p, 1, is_log_scale ? 1.0 : 0.5, text("No successful runs", 10, :center))
        savefig(p, joinpath(run_dir, filename))
        return
    end

    y_max_multiplier = is_log_scale ? 2.5 : 1.3
    y_min = is_log_scale ? 0.5 : 0
    boxplot!(p, ylimits = (y_min, max_val * y_max_multiplier))
    
    annot_y = is_log_scale ? max_val * 1.5 : max_val * 1.15
    for (idx, g) in enumerate(active_groups)
        use_frag, metric, label, col = g
        sub_df = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, succ_rxn_df)
        vals = sub_df[!, col_name]
        if !isempty(vals)
            annotate!(p, xs[idx], annot_y,
                text("Avg: $(round(mean(vals); digits=val_digits))$(val_suffix)", 7, :center, :bottom))
        end
    end

    savefig(p, joinpath(run_dir, filename))
end

# -------------------------------------------------------------
# Save Results & Generate Plots
# -------------------------------------------------------------

"""
    save_benchmark_results(
        run_dir::String,
        results_df::DataFrame,
        dataset::String,
        missing_molecule_synthesis_stats::Dict{Bool, Dict{Symbol, Any}}
    )

Save tractability benchmark results (metadata, raw csv, summary csv, and SVG plots/CSVs).
"""
function save_benchmark_results(
        run_dir::String,
        results_df::DataFrame,
        dataset::String,
        missing_molecule_synthesis_stats::Dict{Bool, Dict{Symbol, Any}}
)
    # Write metadata.txt
    metadata_path = joinpath(run_dir, "metadata.txt")
    open(metadata_path, "w") do io
        println(io, "package: $(get_package_versions())")
        println(io, "task: rbl")
        println(io, "dataset: $dataset")
        println(io, "source: github")
        println(io, "version: v1.0.0")
        println(io, "cache: $(SynRXNLoader._SYNRXN_CACHE_DIR)")
        println(io, "evaluation_script_commit: $(get_git_commit())")
    end

    # Write detailed results CSV
    CSV.write(joinpath(run_dir, "results.csv"), results_df)

    # Serialize raw data for standalone plot generation
    Serialization.serialize(joinpath(run_dir, "benchmark_data.jls"),
        (results_df, dataset, missing_molecule_synthesis_stats))

    generate_plots(run_dir, results_df, dataset, missing_molecule_synthesis_stats)
end

function generate_plots(run_dir::String)
    data = Serialization.deserialize(joinpath(run_dir, "benchmark_data.jls"))
    generate_plots(run_dir, data...)
end

function generate_plots(
        run_dir::String,
        results_df::DataFrame,
        dataset::String,
        missing_molecule_synthesis_stats::Dict{Bool, Dict{Symbol, Any}}
)
    # Compute and write summary CSV
    summary_df = combine(groupby(results_df, [
        :similarity_metric, :use_fragments])) do sdf
        eval_runs = nrow(sdf)
        success_runs = sum(sdf.success)
        success_rate = eval_runs > 0 ? (success_runs / eval_runs) * 100 : 0.0

        # Calculate average time only for successful runs
        success_times = sdf.elapsed_time_seconds[sdf.success]
        avg_time = isempty(success_times) ? 0.0 : sum(success_times) / length(success_times)

        total_time = sum(sdf.elapsed_time_seconds)

        return (
            evaluated_runs = eval_runs,
            successful_runs = success_runs,
            success_rate = round(success_rate; digits = 1),
            avg_synthesis_time_seconds = round(avg_time; digits = 2),
            total_synthesis_time_seconds = round(total_time; digits = 2)
        )
    end
    CSV.write(joinpath(run_dir, "summary.csv"), summary_df)

    # 3. Bar plot & CSV: Success rate of missing molecule synthesis subproblem
    success_with = sum(missing_molecule_synthesis_stats[true][:successes]) /
                   max(1, length(missing_molecule_synthesis_stats[true][:successes])) * 100
    success_without = sum(missing_molecule_synthesis_stats[false][:successes]) /
                      max(1, length(missing_molecule_synthesis_stats[false][:successes])) *
                      100
    synthesis_eval_limit = length(missing_molecule_synthesis_stats[true][:successes])

    succ_count_with = sum(missing_molecule_synthesis_stats[true][:successes])
    succ_count_without = sum(missing_molecule_synthesis_stats[false][:successes])

    success_rate_df = DataFrame(
        Use_Fragments = ["base grammar", "with BRICS fragments"],
        Success_Rate_Percentage = [success_without, success_with],
        Successful_Problems = [succ_count_without, succ_count_with]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_success_rate.csv"), success_rate_df)

    p_success_rate = groupedbar(
        ["$dataset\n($synthesis_eval_limit)"],
        [success_without success_with],
        xlabel = "Dataset",
        ylabel = "Success Rate (%)",
        title = "Missing Molecule Synthesis (SynRXN rbl)\nSuccess Rate",
        label = ["base grammar" "with BRICS fragments"],
        legend = :outerbottom,
        ylimits = (0, 105),
        series_annotations = [Plots.text("$(round(success_without; digits=1))% ($succ_count_without)", 8, :center, :bottom) Plots.text("$(round(success_with; digits=1))% ($succ_count_with)", 8, :center, :bottom)],
        color = ["#E69F00" "#0072B2"],
        edgecolor = :black,
        linewidth = 1.2,
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8,
        size = (600, 450)
    )
    savefig(p_success_rate, joinpath(run_dir, "missing_molecule_synthesis_success_rate.svg"))

    # 4. Sizes of untractable molecules summary CSV
    sizes_with = missing_molecule_synthesis_stats[true][:unfeasible_sizes]
    sizes_without = missing_molecule_synthesis_stats[false][:unfeasible_sizes]

    total_with = length(sizes_with)
    total_without = length(sizes_without)

    pct_le_3_with = total_with > 0 ? count(s -> s <= 3, sizes_with) / total_with * 100 : 0.0
    pct_4_6_with = total_with > 0 ? count(s -> 4 <= s <= 6, sizes_with) / total_with * 100 : 0.0
    pct_7_9_with = total_with > 0 ? count(s -> 7 <= s <= 9, sizes_with) / total_with * 100 : 0.0
    pct_10_13_with = total_with > 0 ? count(s -> 10 <= s <= 13, sizes_with) / total_with * 100 : 0.0
    pct_ge_14_with = total_with > 0 ? count(s -> s >= 14, sizes_with) / total_with * 100 : 0.0

    pct_le_3_without = total_without > 0 ? count(s -> s <= 3, sizes_without) / total_without * 100 : 0.0
    pct_4_6_without = total_without > 0 ? count(s -> 4 <= s <= 6, sizes_without) / total_without * 100 : 0.0
    pct_7_9_without = total_without > 0 ? count(s -> 7 <= s <= 9, sizes_without) / total_without * 100 : 0.0
    pct_10_13_without = total_without > 0 ? count(s -> 10 <= s <= 13, sizes_without) / total_without * 100 : 0.0
    pct_ge_14_without = total_without > 0 ? count(s -> s >= 14, sizes_without) / total_without * 100 : 0.0

    unsynthesised_sizes_df = DataFrame(
        Symbol("dataset") => [dataset, dataset],
        Symbol("with fragments or base grammar") => ["base grammar", "with fragments"],
        Symbol("total number of unsynthesised molecules") => [total_without, total_with],
        Symbol("percentage of missing molecules with sizes <= 3") => [round(pct_le_3_without; digits=1), round(pct_le_3_with; digits=1)],
        Symbol("percentage of missing molecules with sizes 4-6") => [round(pct_4_6_without; digits=1), round(pct_4_6_with; digits=1)],
        Symbol("percentage of missing molecules with sizes 7-9") => [round(pct_7_9_without; digits=1), round(pct_7_9_with; digits=1)],
        Symbol("percentage of missing molecules with sizes 10-13") => [round(pct_10_13_without; digits=1), round(pct_10_13_with; digits=1)],
        Symbol("percentage of missing molecules with sizes >= 14") => [round(pct_ge_14_without; digits=1), round(pct_ge_14_with; digits=1)]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_unsynthesised_sizes.csv"), unsynthesised_sizes_df)

    # 5. Box plot & CSV: Runtime of missing molecule synthesis subproblem
    runtimes_succ_with = missing_molecule_synthesis_stats[true][:runtimes][missing_molecule_synthesis_stats[true][:successes]]
    runtimes_succ_without = missing_molecule_synthesis_stats[false][:runtimes][missing_molecule_synthesis_stats[false][:successes]]

    runtime_with = isempty(runtimes_succ_with) ? 0.0 : mean(runtimes_succ_with)
    runtime_without = isempty(runtimes_succ_without) ? 0.0 : mean(runtimes_succ_without)

    runtime_df = DataFrame(
        Use_Fragments = ["base grammar", "with BRICS fragments"],
        Average_Runtime_Seconds = [runtime_without, runtime_with]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_average_runtime.csv"), runtime_df)

    p_runtime = boxplot(
        xticks = (1:1, [dataset]),
        xlims = (0.5, 1.5),
        xlabel = "Dataset",
        ylabel = "Runtime (s)",
        title = "Missing Molecule Synthesis (SynRXN rbl)\nRuntime of Successful Molecule Syntheses",
        legend = :outerbottom,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8
    )
    # Legend
    boxplot!(p_runtime, [0], [0], label = "base grammar", color = "#E69F00",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black)
    boxplot!(p_runtime, [0], [0], label = "with BRICS fragments", color = "#0072B2",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black)

    if !isempty(runtimes_succ_without)
        boxplot!(p_runtime, fill(1 - 0.22, length(runtimes_succ_without)), runtimes_succ_without,
            label = "", color = "#E69F00", fillalpha = 0.7, bar_width = 0.44, linecolor = :black)
    end
    if !isempty(runtimes_succ_with)
        boxplot!(p_runtime, fill(1 + 0.22, length(runtimes_succ_with)), runtimes_succ_with,
            label = "", color = "#0072B2", fillalpha = 0.7, bar_width = 0.44, linecolor = :black)
    end

    max_rt_all = max(isempty(runtimes_succ_without) ? 0.0 : maximum(runtimes_succ_without),
        isempty(runtimes_succ_with) ? 0.0 : maximum(runtimes_succ_with))
    max_rt_all = max_rt_all == 0.0 ? 1.0 : max_rt_all
    boxplot!(p_runtime, ylims = (0, max_rt_all * 1.3))

    annotate!(p_runtime,
        1 - 0.22,
        max_rt_all * 1.15,
        isempty(runtimes_succ_without) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(runtime_without; digits=2))s", 8, :center, :bottom))
    annotate!(p_runtime,
        1 + 0.22,
        max_rt_all * 1.15,
        isempty(runtimes_succ_with) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(runtime_with; digits=2))s", 8, :center, :bottom))

    savefig(p_runtime, joinpath(run_dir, "missing_molecule_synthesis_average_runtime.svg"))

    # 6. Box plot & CSV: Number of molecules synthesized
    mols_succ_with = missing_molecule_synthesis_stats[true][:molecules_synthesized][missing_molecule_synthesis_stats[true][:successes]]
    mols_succ_without = missing_molecule_synthesis_stats[false][:molecules_synthesized][missing_molecule_synthesis_stats[false][:successes]]

    avg_mols_with = isempty(mols_succ_with) ? 0.0 : mean(mols_succ_with)
    avg_mols_without = isempty(mols_succ_without) ? 0.0 : mean(mols_succ_without)

    mols_df = DataFrame(
        Use_Fragments = ["base grammar", "with BRICS fragments"],
        Average_Molecules = [avg_mols_without, avg_mols_with]
    )
    CSV.write(joinpath(run_dir, "average_molecules_synthesized.csv"), mols_df)

    p_mols = boxplot(
        xticks = (1:1, [dataset]),
        xlims = (0.5, 1.5),
        xlabel = "Dataset",
        ylabel = "Molecules Synthesised",
        title = "Missing Molecule Synthesis (SynRXN rbl)\nNumber of Molecules Synthesised until All Targets Found for Successful Problems",
        legend = :outerbottom,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8
    )
    # Legend
    boxplot!(p_mols, [0], [0], label = "base grammar", color = "#E69F00",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black)
    boxplot!(p_mols, [0], [0], label = "with BRICS fragments", color = "#0072B2",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black)

    if !isempty(mols_succ_without)
        boxplot!(p_mols, fill(1 - 0.22, length(mols_succ_without)), mols_succ_without,
            label = "", color = "#E69F00", fillalpha = 0.7, bar_width = 0.44, linecolor = :black)
    end
    if !isempty(mols_succ_with)
        boxplot!(p_mols, fill(1 + 0.22, length(mols_succ_with)), mols_succ_with,
            label = "", color = "#0072B2", fillalpha = 0.7, bar_width = 0.44, linecolor = :black)
    end

    max_mols_all = max(isempty(mols_succ_without) ? 0.0 : maximum(mols_succ_without),
        isempty(mols_succ_with) ? 0.0 : maximum(mols_succ_with))
    max_mols_all = max_mols_all == 0.0 ? 1.0 : max_mols_all
    boxplot!(p_mols, ylims = (0, max_mols_all * 1.3))

    annotate!(p_mols,
        1 - 0.22,
        max_mols_all * 1.15,
        isempty(mols_succ_without) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(avg_mols_without; digits=1))", 8, :center, :bottom))
    annotate!(p_mols,
        1 + 0.22,
        max_mols_all * 1.15,
        isempty(mols_succ_with) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(avg_mols_with; digits=1))", 8, :center, :bottom))

    savefig(p_mols, joinpath(run_dir, "average_molecules_synthesized.svg"))

    # 7. Box plots & CSV: Reactions synthesized per metric
    metrics_order = ["none", "tanimoto"]
    succ_rxn_df = filter(r -> r.success == true, results_df)

    avg_rxns_with = Float64[]
    avg_rxns_without = Float64[]
    for m in metrics_order
        sub_df_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        push!(avg_rxns_with, isempty(sub_df_with) ? 0.0 : mean(sub_df_with.reactions_synthesized))
        
        sub_df_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        push!(avg_rxns_without, isempty(sub_df_without) ? 0.0 : mean(sub_df_without.reactions_synthesized))
    end
    df_rxns_comp = DataFrame(
        Similarity_Metric = metrics_order,
        Average_Reactions_Without_Fragments = avg_rxns_without,
        Average_Reactions_With_Fragments = avg_rxns_with
    )
    CSV.write(joinpath(run_dir, "average_reactions_synthesized.csv"), df_rxns_comp)

    plot_comparison_boxplot_helper(
        run_dir, succ_rxn_df, metrics_order, :reactions_synthesized, "Reactions Synthesised",
        "average_reactions_synthesized.svg",
        "Reaction Rebalancing (SynRXN rbl)\nNumber of Reactions Synthesised Until Target Found for Successful Problems",
        dataset;
        is_log_scale=true, val_digits=1, val_suffix=""
    )

    # 8. Success Rate Comparison Plots
    # All metrics combined
    plot_success_rate_helper(
        run_dir, summary_df, metrics_order, "synthesis_success_rate.svg",
        dataset, synthesis_eval_limit
    )
    # 9. Runtime Comparison Plots
    # All metrics combined
    plot_comparison_boxplot_helper(
        run_dir, succ_rxn_df, metrics_order, :elapsed_time_seconds, "Synthesis Time (s)",
        "synthesis_average_time.svg",
        "Reaction Rebalancing (SynRXN rbl)\nRuntime of Successful Reaction Syntheses",
        dataset;
        is_log_scale=false, val_digits=2, val_suffix="s"
    )

    # 11. Save corresponding CSV comparison tables
    # For success rate and average synthesis time (all metrics combined)
    success_with_metric = Float64[]
    success_without_metric = Float64[]
    success_counts_with = Int[]
    success_counts_without = Int[]
    time_with_metric = Float64[]
    time_without_metric = Float64[]

    for m in metrics_order
        row_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, summary_df)
        push!(success_with_metric, isempty(row_with) ? 0.0 : row_with[1, :success_rate])
        push!(success_counts_with, isempty(row_with) ? 0 : row_with[1, :successful_runs])
        push!(time_with_metric, isempty(row_with) ? 0.0 :
                                row_with[1, :avg_synthesis_time_seconds])

        row_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, summary_df)
        push!(success_without_metric, isempty(row_without) ? 0.0 :
                                      row_without[1, :success_rate])
        push!(success_counts_without, isempty(row_without) ? 0 :
                                      row_without[1, :successful_runs])
        push!(time_without_metric, isempty(row_without) ? 0.0 :
                                   row_without[1, :avg_synthesis_time_seconds])
    end

    success_comparison_df = DataFrame(
        Similarity_Metric = metrics_order,
        Success_Rate_Without_Fragments_Percentage = success_without_metric,
        Successful_Runs_Without_Fragments = success_counts_without,
        Success_Rate_With_Fragments_Percentage = success_with_metric,
        Successful_Runs_With_Fragments = success_counts_with
    )
    CSV.write(joinpath(run_dir, "synthesis_success_rate.csv"), success_comparison_df)

    time_comparison_df = DataFrame(
        Similarity_Metric = metrics_order,
        Average_Synthesis_Time_Without_Fragments_Seconds = time_without_metric,
        Average_Synthesis_Time_With_Fragments_Seconds = time_with_metric
    )
    CSV.write(joinpath(run_dir, "synthesis_average_time.csv"), time_comparison_df)

    println("\n=======================================================")
    println("✓ Benchmark results successfully recorded!")
    println("  - Run Folder:  $run_dir")
    println("=======================================================")
end
