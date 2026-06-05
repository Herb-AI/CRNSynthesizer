using CRNSynthesizer
using DataStructures
import TOML
using Dates
using DataFrames
using CSV
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

# -------------------------------------------------------------
# Plotting Helpers to Maximize Code Reuse
# -------------------------------------------------------------

function create_pie_plot(labels_raw::Vector{String}, counts::Vector{Int}, title_str::String)
    total = sum(counts)
    if total == 0
        return pie(title = title_str, legend = false,
            grid = false, xaxis = false, yaxis = false,
            annotation = (0.5, 0.5, "No data available"))
    end

    # Create pie chart with no legend
    p = pie(labels_raw, counts,
        title = title_str,
        legend = false,
        titlefontsize = 9,
        size = (600, 450)
    )

    # Calculate slice angles to annotate directly on the slices
    θ = (cumsum(counts) - counts/2) .* 360 / total
    for (i, θ_i) in enumerate(θ)
        s, c = sincosd(θ_i)
        lbl = labels_raw[i]
        val = counts[i]
        pct = round(val / total * 100; digits = 1)
        annotate!(p, 0.6 * c, 0.6 * s, text("$lbl\n$val ($pct%)", 8, :black, :center))
    end
    return p
end

function make_size_dist_df(sizes::Vector{Int})
    counts_dict = Dict{Int, Int}()
    for s in sizes
        counts_dict[s] = get(counts_dict, s, 0) + 1
    end
    sorted_pairs = sort(collect(counts_dict); by = x -> x[1])
    lbls = [string(x[1]) for x in sorted_pairs]
    vals = [x[2] for x in sorted_pairs]
    total = sum(vals)
    pcts = total > 0 ? [v / total * 100 for v in vals] : Float64[]
    return lbls, vals,
    DataFrame(Molecule_Size = [x[1] for x in sorted_pairs], Count = vals, Percentage = pcts)
end

"""
    clean_plot_title(title::String) -> String

Remove :none/none references and trailing hyphens from plot titles.
"""
function clean_plot_title(title::String)::String
    t = replace(title, " (:none)" => "")
    t = replace(t, ":none" => "")
    t = replace(t, " none" => "")
    t = replace(t, "none" => "")
    t = replace(t, r"\s*-\s*\n" => "\n")
    t = replace(t, r"\s*-\s*$" => "")
    return t
end

function plot_success_rate_helper(
        run_dir::String, summary_df::DataFrame, metrics::Vector{String},
        filename::String, title::String, eval_limit::Int)
    clean_title = clean_plot_title(title)

    if length(metrics) == 1
        m = metrics[1]
        row_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, summary_df)
        val_with = isempty(row_with) ? 0.0 : row_with[1, :success_rate]
        count_with = isempty(row_with) ? 0 : row_with[1, :successful_runs]

        row_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, summary_df)
        val_without = isempty(row_without) ? 0.0 : row_without[1, :success_rate]
        count_without = isempty(row_without) ? 0 : row_without[1, :successful_runs]

        p = bar(
            ["Without Fragments", "With Fragments"],
            [val_without, val_with],
            xlabel = "With and Without Fragments",
            ylabel = "Success Rate (%)",
            title = clean_title,
            legend = false,
            ylimits = (0, 105),
            series_annotations = map(
                (v, c) -> Plots.text("$(round(v; digits=1))%\n($c)", 8, :center, :bottom),
                [val_without, val_with], [count_without, count_with]),
            color = [:lightgrey, :dodgerblue],
            linewidth = 1.2,
            edgecolor = :black,
            titlefontsize = 9,
            guidefontsize = 9,
            tickfontsize = 8,
            size = (600, 450)
        )
        savefig(p, joinpath(run_dir, filename))
    else
        success_with_metric = Float64[]
        success_without_metric = Float64[]
        success_counts_with = Int[]
        success_counts_without = Int[]

        for m in metrics
            row_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, summary_df)
            push!(success_with_metric, isempty(row_with) ? 0.0 : row_with[1, :success_rate])
            push!(success_counts_with, isempty(row_with) ? 0 :
                                       row_with[1, :successful_runs])

            row_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, summary_df)
            push!(success_without_metric, isempty(row_without) ? 0.0 :
                                          row_without[1, :success_rate])
            push!(success_counts_without, isempty(row_without) ? 0 :
                                          row_without[1, :successful_runs])
        end

        success_matrix = hcat(success_without_metric, success_with_metric)
        success_counts_matrix = hcat(success_counts_without, success_counts_with)
        p = groupedbar(metrics, success_matrix,
            xlabel = "Similarity Guidance Metric",
            ylabel = "Success Rate (%)",
            label = ["Without Fragments" "With Fragments"],
            title = clean_title,
            legend = :right,
            ylimits = (0, 105),
            series_annotations = map(
                (v, c) -> Plots.text("$(round(v; digits=1))%\n($c)", 7, :center, :bottom), success_matrix, success_counts_matrix),
            color = [:lightgrey :dodgerblue],
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
end

function plot_synthesis_time_helper(run_dir::String, succ_rxn_df::DataFrame,
        metrics::Vector{String}, filename::String, title::String)
    plot_df = filter(r -> r.similarity_metric in metrics, succ_rxn_df)
    num_success = length(unique(plot_df.r_id))

    clean_title = clean_plot_title(title)
    clean_title = replace(clean_title, "\n" => " - $num_success Successful Problems\n")

    if length(metrics) == 1
        p = boxplot(
            xticks = (1:2, ["Without Fragments", "With Fragments"]),
            xlims = (0.5, 2.5),
            xlabel = "With and Without Fragments",
            ylabel = "Synthesis Time (s)",
            title = clean_title,
            legend = false,
            size = (600, 450),
            titlefontsize = 9,
            guidefontsize = 9,
            tickfontsize = 8
        )

        vals_without = filter(
            r -> r.similarity_metric == metrics[1] &&
                 r.use_fragments == false, succ_rxn_df).elapsed_time_seconds
        if !isempty(vals_without)
            boxplot!(p, fill(1, length(vals_without)), vals_without, color = :lightgrey,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
        end

        vals_with = filter(
            r -> r.similarity_metric == metrics[1] &&
                 r.use_fragments == true, succ_rxn_df).elapsed_time_seconds
        if !isempty(vals_with)
            boxplot!(p, fill(2, length(vals_with)), vals_with, color = :dodgerblue,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
        end

        max_val = max(isempty(vals_without) ? 0.0 : maximum(vals_without),
            isempty(vals_with) ? 0.0 : maximum(vals_with))
        max_val = max_val == 0.0 ? 1.0 : max_val
        boxplot!(p, ylimits = (0, max_val * 1.3))

        annotate!(p,
            1,
            max_val * 1.15,
            isempty(vals_without) ? text("N/A", 8, :center, :bottom) :
            text("Avg: $(round(mean(vals_without); digits=2))s", 8, :center, :bottom))
        annotate!(p,
            2,
            max_val * 1.15,
            isempty(vals_with) ? text("N/A", 8, :center, :bottom) :
            text("Avg: $(round(mean(vals_with); digits=2))s", 8, :center, :bottom))
        savefig(p, joinpath(run_dir, filename))
    else
        p = boxplot(
            xticks = (1:length(metrics), metrics),
            xlims = (0.5, length(metrics) + 0.5),
            xlabel = "Similarity Guidance Metric",
            ylabel = "Synthesis Time (s)",
            title = clean_title,
            legend = :topright,
            size = (600, 450),
            titlefontsize = 9,
            guidefontsize = 9,
            tickfontsize = 8,
            legendfontsize = 8
        )
        boxplot!(p, [0], [0], label = "Without Fragments", color = :lightgrey,
            seriestype = :shape, fillalpha = 0.7, linecolor = :black)
        boxplot!(p, [0], [0], label = "With Fragments", color = :dodgerblue,
            seriestype = :shape, fillalpha = 0.7, linecolor = :black)

        max_time_val = 0.0
        for (i, m) in enumerate(metrics)
            sub_df_without = filter(
                r -> r.similarity_metric == m &&
                     r.use_fragments == false, succ_rxn_df)
            vals_without = sub_df_without.elapsed_time_seconds
            if !isempty(vals_without)
                boxplot!(p, fill(i - 0.22, length(vals_without)),
                    vals_without, label = "", color = :lightgrey,
                    fillalpha = 0.7, bar_width = 0.35, linecolor = :black)
                max_time_val = max(max_time_val, maximum(vals_without))
            end

            sub_df_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
            vals_with = sub_df_with.elapsed_time_seconds
            if !isempty(vals_with)
                boxplot!(p, fill(i + 0.22, length(vals_with)), vals_with,
                    label = "", color = :dodgerblue, fillalpha = 0.7,
                    bar_width = 0.35, linecolor = :black)
                max_time_val = max(max_time_val, maximum(vals_with))
            end
        end

        max_time_val = max_time_val == 0.0 ? 1.0 : max_time_val
        boxplot!(p, ylimits = (0, max_time_val * 1.3))
        savefig(p, joinpath(run_dir, filename))
    end
end

function plot_reactions_synthesized_comparison_helper(
        run_dir::String, succ_rxn_df::DataFrame,
        metrics::Vector{String}, filename::String, title::String)
    plot_df = filter(r -> r.similarity_metric in metrics, succ_rxn_df)
    num_success = length(unique(plot_df.r_id))

    clean_title = clean_plot_title(title)
    clean_title = replace(clean_title, "\n" => " - $num_success Successful Problems\n")

    if length(metrics) == 1
        p = boxplot(
            xticks = (1:2, ["Without Fragments", "With Fragments"]),
            xlims = (0.5, 2.5),
            xlabel = "With and Without Fragments",
            ylabel = "Reactions Synthesized",
            title = clean_title,
            legend = false,
            size = (600, 450),
            titlefontsize = 9,
            guidefontsize = 9,
            tickfontsize = 8,
            yscale = :log10
        )

        vals_without = filter(
            r -> r.similarity_metric == metrics[1] &&
                 r.use_fragments == false, succ_rxn_df).reactions_synthesized
        if !isempty(vals_without)
            boxplot!(p, fill(1, length(vals_without)), vals_without, color = :lightgrey,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
        end

        vals_with = filter(
            r -> r.similarity_metric == metrics[1] &&
                 r.use_fragments == true, succ_rxn_df).reactions_synthesized
        if !isempty(vals_with)
            boxplot!(p, fill(2, length(vals_with)), vals_with, color = :dodgerblue,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
        end

        max_val = max(isempty(vals_without) ? 0.0 : maximum(vals_without),
            isempty(vals_with) ? 0.0 : maximum(vals_with))
        max_val = max_val == 0.0 ? 1.0 : max_val
        boxplot!(p, ylimits = (0.5, max_val * 2.5))

        annotate!(p,
            1,
            max_val * 1.5,
            isempty(vals_without) ? text("N/A", 8, :center, :bottom) :
            text("Avg: $(round(mean(vals_without); digits=1))", 8, :center, :bottom))
        annotate!(p,
            2,
            max_val * 1.5,
            isempty(vals_with) ? text("N/A", 8, :center, :bottom) :
            text("Avg: $(round(mean(vals_with); digits=1))", 8, :center, :bottom))
        savefig(p, joinpath(run_dir, filename))
    else
        p = boxplot(
            xticks = (1:length(metrics), metrics),
            xlims = (0.5, length(metrics) + 0.5),
            xlabel = "Similarity Guidance Metric",
            ylabel = "Reactions Synthesized",
            title = clean_title,
            legend = :topright,
            size = (600, 450),
            titlefontsize = 9,
            guidefontsize = 9,
            tickfontsize = 8,
            legendfontsize = 8,
            yscale = :log10
        )
        boxplot!(p, [0], [0], label = "Without Fragments", color = :lightgrey,
            seriestype = :shape, fillalpha = 0.7, linecolor = :black)
        boxplot!(p, [0], [0], label = "With Fragments", color = :dodgerblue,
            seriestype = :shape, fillalpha = 0.7, linecolor = :black)

        max_rxn_val = 0.0
        for (i, m) in enumerate(metrics)
            sub_df_without = filter(
                r -> r.similarity_metric == m &&
                     r.use_fragments == false, succ_rxn_df)
            vals_without = sub_df_without.reactions_synthesized
            if !isempty(vals_without)
                boxplot!(p, fill(i - 0.22, length(vals_without)),
                    vals_without, label = "", color = :lightgrey,
                    fillalpha = 0.7, bar_width = 0.35, linecolor = :black)
                max_rxn_val = max(max_rxn_val, maximum(vals_without))
            end

            sub_df_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
            vals_with = sub_df_with.reactions_synthesized
            if !isempty(vals_with)
                boxplot!(p, fill(i + 0.22, length(vals_with)), vals_with,
                    label = "", color = :dodgerblue, fillalpha = 0.7,
                    bar_width = 0.35, linecolor = :black)
                max_rxn_val = max(max_rxn_val, maximum(vals_with))
            end
        end

        max_rxn_val = max_rxn_val == 0.0 ? 1.0 : max_rxn_val
        boxplot!(p, ylimits = (0.5, max_rxn_val * 2.5))

        for (i, m) in enumerate(metrics)
            sub_df_without = filter(
                r -> r.similarity_metric == m &&
                     r.use_fragments == false, succ_rxn_df)
            vals_without = sub_df_without.reactions_synthesized
            if !isempty(vals_without)
                annotate!(p,
                    i - 0.22,
                    max_rxn_val * 1.5,
                    text("Avg: $(round(mean(vals_without); digits=1))", 7, :center, :bottom))
            end

            sub_df_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
            vals_with = sub_df_with.reactions_synthesized
            if !isempty(vals_with)
                annotate!(p, i + 0.22, max_rxn_val * 1.5,
                    text("Avg: $(round(mean(vals_with); digits=1))", 7, :center, :bottom))
            end
        end

        savefig(p, joinpath(run_dir, filename))
    end
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
        Use_Fragments = ["Without Fragments", "With Fragments"],
        Success_Rate_Percentage = [success_without, success_with],
        Successful_Problems = [succ_count_without, succ_count_with]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_success_rate.csv"), success_rate_df)

    p_success_rate = bar(
        ["Without Fragments", "With Fragments"],
        [success_without, success_with],
        xlabel = "With and Without Fragments",
        ylabel = "Success Rate (%)",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset)\nSuccess Rate ($synthesis_eval_limit Problems Analyzed)",
        legend = false,
        ylimits = (0, 105),
        series_annotations = map(
            (v, c) -> Plots.text("$(round(v; digits=1))% ($c)", 8, :center, :bottom),
            [success_without, success_with], [succ_count_without, succ_count_with]),
        color = [:lightgrey, :dodgerblue],
        edgecolor = :black,
        linewidth = 1.2,
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        size = (600, 450)
    )
    savefig(p_success_rate, joinpath(run_dir, "missing_molecule_synthesis_success_rate.svg"))

    # 4. Pie plots & CSV: Sizes of untractable molecules based on missing molecule synthesis subproblem
    # With fragments
    sizes_with_lbls, sizes_with_vals,
    sizes_with_df = make_size_dist_df(missing_molecule_synthesis_stats[true][:unfeasible_sizes])
    CSV.write(
        joinpath(run_dir, "missing_molecule_synthesis_untractable_sizes_with_fragments.csv"),
        sizes_with_df)

    failed_with_count = sum(.!missing_molecule_synthesis_stats[true][:successes])
    p_sizes_with = create_pie_plot(
        sizes_with_lbls,
        sizes_with_vals,
        "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $failed_with_count Problems Failed\nUnsynthesized Molecule Sizes (With Fragments)"
    )
    savefig(p_sizes_with,
        joinpath(run_dir, "missing_molecule_synthesis_untractable_sizes_with_fragments.svg"))

    # Without fragments
    sizes_without_lbls, sizes_without_vals,
    sizes_without_df = make_size_dist_df(missing_molecule_synthesis_stats[false][:unfeasible_sizes])
    CSV.write(
        joinpath(run_dir, "missing_molecule_synthesis_untractable_sizes_without_fragments.csv"),
        sizes_without_df)

    failed_without_count = sum(.!missing_molecule_synthesis_stats[false][:successes])
    p_sizes_without = create_pie_plot(
        sizes_without_lbls,
        sizes_without_vals,
        "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $failed_without_count Problems Failed\nUnsynthesized Molecule Sizes (Without Fragments)"
    )
    savefig(p_sizes_without,
        joinpath(run_dir, "missing_molecule_synthesis_untractable_sizes_without_fragments.svg"))

    # 4b. Combined CSV: unfeasible molecule sizes comparison
    sizes_with = missing_molecule_synthesis_stats[true][:unfeasible_sizes]
    sizes_without = missing_molecule_synthesis_stats[false][:unfeasible_sizes]
    unique_sizes = sort(unique(vcat(sizes_with, sizes_without)))

    total_with = length(sizes_with)
    total_without = length(sizes_without)

    comp_rows = []
    for sz in unique_sizes
        cnt_with = count(==(sz), sizes_with)
        pct_with = total_with > 0 ? (cnt_with / total_with) * 100 : 0.0

        cnt_without = count(==(sz), sizes_without)
        pct_without = total_without > 0 ? (cnt_without / total_without) * 100 : 0.0

        push!(comp_rows,
            (
                Molecule_Size = sz,
                Count_With_Fragments = cnt_with,
                Percentage_With_Fragments = round(pct_with; digits = 1),
                Count_Without_Fragments = cnt_without,
                Percentage_Without_Fragments = round(pct_without; digits = 1)
            ))
    end
    comp_df = DataFrame(comp_rows)
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_untractable_sizes_comparison.csv"), comp_df)

    # Calculate count of successful problems for molecules/reactions
    num_success_mols = sum(missing_molecule_synthesis_stats[true][:successes] .|
                           missing_molecule_synthesis_stats[false][:successes])

    # 5. Box plot & CSV: Runtime of missing molecule synthesis subproblem
    runtimes_succ_with = missing_molecule_synthesis_stats[true][:runtimes][missing_molecule_synthesis_stats[true][:successes]]
    runtimes_succ_without = missing_molecule_synthesis_stats[false][:runtimes][missing_molecule_synthesis_stats[false][:successes]]

    runtime_with = isempty(runtimes_succ_with) ? 0.0 : mean(runtimes_succ_with)
    runtime_without = isempty(runtimes_succ_without) ? 0.0 : mean(runtimes_succ_without)

    runtime_df = DataFrame(
        Use_Fragments = ["Without Fragments", "With Fragments"],
        Average_Runtime_Seconds = [runtime_without, runtime_with]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_average_runtime.csv"), runtime_df)

    p_runtime = boxplot(
        xticks = (1:2, ["Without Fragments", "With Fragments"]),
        xlims = (0.5, 2.5),
        xlabel = "With and Without Fragments",
        ylabel = "Runtime (s)",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $num_success_mols Successful Problems\nRuntime of Successful Syntheses",
        legend = false,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8
    )
    if !isempty(runtimes_succ_without)
        boxplot!(p_runtime, fill(1, length(runtimes_succ_without)), runtimes_succ_without,
            color = :lightgrey, fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
    end
    if !isempty(runtimes_succ_with)
        boxplot!(p_runtime, fill(2, length(runtimes_succ_with)), runtimes_succ_with,
            color = :dodgerblue, fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
    end

    max_rt_all = max(isempty(runtimes_succ_without) ? 0.0 : maximum(runtimes_succ_without),
        isempty(runtimes_succ_with) ? 0.0 : maximum(runtimes_succ_with))
    max_rt_all = max_rt_all == 0.0 ? 1.0 : max_rt_all
    boxplot!(p_runtime, ylims = (0, max_rt_all * 1.3))

    annotate!(p_runtime,
        1,
        max_rt_all * 1.15,
        isempty(runtimes_succ_without) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(runtime_without; digits=2))s", 8, :center, :bottom))
    annotate!(p_runtime,
        2,
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
        Use_Fragments = ["Without Fragments", "With Fragments"],
        Average_Molecules = [avg_mols_without, avg_mols_with]
    )
    CSV.write(joinpath(run_dir, "average_molecules_synthesized.csv"), mols_df)

    p_mols = boxplot(
        xticks = (1:2, ["Without Fragments", "With Fragments"]),
        xlims = (0.5, 2.5),
        xlabel = "With and Without Fragments",
        ylabel = "Molecules Synthesized",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $num_success_mols Successful Problems\nNumber of Molecules Synthesized for Successful Problems",
        legend = false,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8
    )
    if !isempty(mols_succ_without)
        boxplot!(p_mols, fill(1, length(mols_succ_without)), mols_succ_without,
            color = :lightgrey, fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
    end
    if !isempty(mols_succ_with)
        boxplot!(p_mols, fill(2, length(mols_succ_with)), mols_succ_with,
            color = :dodgerblue, fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
    end

    max_mols_all = max(isempty(mols_succ_without) ? 0.0 : maximum(mols_succ_without),
        isempty(mols_succ_with) ? 0.0 : maximum(mols_succ_with))
    max_mols_all = max_mols_all == 0.0 ? 1.0 : max_mols_all
    boxplot!(p_mols, ylims = (0, max_mols_all * 1.3))

    annotate!(p_mols,
        1,
        max_mols_all * 1.15,
        isempty(mols_succ_without) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(avg_mols_without; digits=1))", 8, :center, :bottom))
    annotate!(p_mols,
        2,
        max_mols_all * 1.15,
        isempty(mols_succ_with) ? text("N/A", 8, :center, :bottom) :
        text("Avg: $(round(avg_mols_with; digits=1))", 8, :center, :bottom))

    savefig(p_mols, joinpath(run_dir, "average_molecules_synthesized.svg"))

    # 7. Box plots & CSV: Reactions synthesized per metric
    metrics_order = ["none", "simpson", "tanimoto", "both"]
    succ_rxn_df = filter(r -> r.success == true, results_df)
    num_success_rxns = length(unique(succ_rxn_df.r_id))

    # With fragments
    avg_rxns_with = Float64[]
    for m in metrics_order
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        push!(avg_rxns_with, isempty(sub_df) ? 0.0 : mean(sub_df.reactions_synthesized))
    end
    df_rxns_with = DataFrame(
        Similarity_Metric = metrics_order,
        Average_Reactions = avg_rxns_with
    )
    CSV.write(joinpath(run_dir, "average_reactions_synthesized_with_fragments.csv"), df_rxns_with)

    p_rxns_with = boxplot(
        xticks = (1:4, metrics_order),
        xlims = (0.5, 4.5),
        xlabel = "Similarity Guidance Metric",
        ylabel = "Reactions Synthesized",
        title = "Reaction Rebalancing (SynRXN rbl/$dataset) - $num_success_rxns Successful Problems\nNumber of Reactions Synthesized (With Fragments)",
        legend = false,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        yscale = :log10
    )
    max_rxn_with_val = 0.0
    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        if !isempty(vals)
            boxplot!(p_rxns_with, fill(i, length(vals)), vals, color = :dodgerblue,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
            max_rxn_with_val = max(max_rxn_with_val, maximum(vals))
        end
    end
    max_rxn_with_val = max_rxn_with_val == 0.0 ? 1.0 : max_rxn_with_val
    boxplot!(p_rxns_with, ylims = (0.5, max_rxn_with_val * 2.5))

    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        annot = isempty(vals) ? text("N/A", 8, :center, :bottom) :
                text("Avg: $(round(mean(vals); digits=1))", 8, :center, :bottom)
        annotate!(p_rxns_with, i, max_rxn_with_val * 1.5, annot)
    end
    savefig(p_rxns_with, joinpath(run_dir, "average_reactions_synthesized_with_fragments.svg"))

    # Without fragments
    avg_rxns_without = Float64[]
    for m in metrics_order
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        push!(avg_rxns_without, isempty(sub_df) ? 0.0 : mean(sub_df.reactions_synthesized))
    end
    df_rxns_without = DataFrame(
        Similarity_Metric = metrics_order,
        Average_Reactions = avg_rxns_without
    )
    CSV.write(joinpath(run_dir, "average_reactions_synthesized_without_fragments.csv"), df_rxns_without)

    p_rxns_without = boxplot(
        xticks = (1:4, metrics_order),
        xlims = (0.5, 4.5),
        xlabel = "Similarity Guidance Metric",
        ylabel = "Reactions Synthesized",
        title = "Reaction Rebalancing (SynRXN rbl/$dataset) - $num_success_rxns Successful Problems\nNumber of Reactions Synthesized (Without Fragments)",
        legend = false,
        size = (600, 450),
        titlefontsize = 9,
        guidefontsize = 9,
        tickfontsize = 8,
        yscale = :log10
    )
    max_rxn_without_val = 0.0
    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        if !isempty(vals)
            boxplot!(p_rxns_without, fill(i, length(vals)), vals, color = :lightgrey,
                fillalpha = 0.7, bar_width = 0.6, linecolor = :black)
            max_rxn_without_val = max(max_rxn_without_val, maximum(vals))
        end
    end
    max_rxn_without_val = max_rxn_without_val == 0.0 ? 1.0 : max_rxn_without_val
    boxplot!(p_rxns_without, ylims = (0.5, max_rxn_without_val * 2.5))

    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        annot = isempty(vals) ? text("N/A", 8, :center, :bottom) :
                text("Avg: $(round(mean(vals); digits=1))", 8, :center, :bottom)
        annotate!(p_rxns_without, i, max_rxn_without_val * 1.5, annot)
    end
    savefig(p_rxns_without, joinpath(run_dir, "average_reactions_synthesized_without_fragments.svg"))

    # 8. Success Rate Comparison Plots
    # All metrics combined
    plot_success_rate_helper(
        run_dir, summary_df, metrics_order, "synthesis_success_rate.svg",
        "Reaction Rebalancing (SynRXN rbl/$dataset)\nSuccess Rate ($synthesis_eval_limit Problems Analyzed)",
        synthesis_eval_limit
    )
    # :none metric only
    plot_success_rate_helper(
        run_dir, summary_df, ["none"], "synthesis_success_rate_none.svg",
        "Reaction Rebalancing (SynRXN rbl/$dataset)\nSuccess Rate with No Similarity Guidance",
        synthesis_eval_limit
    )

    # 9. Runtime Comparison Plots
    # All metrics combined
    plot_synthesis_time_helper(
        run_dir, succ_rxn_df, metrics_order, "synthesis_average_time.svg",
        "Reaction Rebalancing (SynRXN rbl/$dataset)\nRuntime of Successful Reaction Syntheses"
    )
    # :none metric only
    plot_synthesis_time_helper(
        run_dir, succ_rxn_df, ["none"], "synthesis_average_time_none.svg",
        "Reaction Rebalancing (SynRXN rbl/$dataset)\nRuntime with No Similarity Guidance"
    )

    # 10. Reactions Synthesized Comparison Plots
    # :none metric only
    plot_reactions_synthesized_comparison_helper(
        run_dir, succ_rxn_df, ["none"], "average_reactions_synthesized_none.svg",
        "Reaction Rebalancing (SynRXN rbl/$dataset)\nReactions Synthesized with No Similarity Guidance"
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

    # For :none metric specifically
    none_row_with = filter(r -> r.similarity_metric == "none" && r.use_fragments == true, summary_df)
    none_row_without = filter(
        r -> r.similarity_metric == "none" &&
             r.use_fragments == false, summary_df)

    succ_rxn_none = filter(r -> r.similarity_metric == "none", succ_rxn_df)
    none_rxns_with = filter(r -> r.use_fragments == true, succ_rxn_none)
    none_rxns_without = filter(r -> r.use_fragments == false, succ_rxn_none)

    none_comp_df = DataFrame(
        Fragments_Setting = ["Without Fragments", "With Fragments"],
        Success_Rate_Percentage = [
            isempty(none_row_without) ? 0.0 : none_row_without[1, :success_rate],
            isempty(none_row_with) ? 0.0 : none_row_with[1, :success_rate]
        ],
        Average_Time_Seconds = [
            isempty(none_row_without) ? 0.0 :
            none_row_without[1, :avg_synthesis_time_seconds],
            isempty(none_row_with) ? 0.0 : none_row_with[1, :avg_synthesis_time_seconds]
        ],
        Average_Reactions_Synthesized = [
            isempty(none_rxns_without) ? 0.0 :
            mean(none_rxns_without.reactions_synthesized),
            isempty(none_rxns_with) ? 0.0 : mean(none_rxns_with.reactions_synthesized)
        ]
    )
    CSV.write(joinpath(run_dir, "reaction_synthesis_none_comparison.csv"), none_comp_df)

    println("\n=======================================================")
    println("✓ Benchmark results successfully recorded!")
    println("  - Run Folder:  $run_dir")
    println("=======================================================")
end
