using CRNSynthesizer
import TOML
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


function get_boxplot_positions_and_width(i::Int, N::Int; total_span::Float64 = 0.44)
    W = total_span / N
    xs = [i - total_span / 2 + (k - 0.5) * W for k in 1:N]
    return xs, W
end

function plot_success_rate_helper(
        output_dir::String, datasets::Vector{String}, num_problems::Vector{Int},
        summary_dfs::Vector{DataFrame}, metrics::Vector{String},
        filename::String;
        is_paper::Bool = false,
        title_str::String = "Reaction Rebalancing (SynRXN rbl)\nSuccess Rate"
)
    all_groups = [
        (false, "none", "Base, no similarity"),
        (false, "tanimoto", "Base, Tanimoto"),
        (true, "none", "w/ BRICS, no similarity"),
        (true, "tanimoto", "w/ BRICS, Tanimoto")
    ]
    
    active_groups = filter(g -> g[2] in metrics, all_groups)
    N = length(active_groups)
    M = length(datasets)
    
    val_matrix = zeros(M, N)
    count_matrix = zeros(Int, M, N)
    
    for j in 1:M
        summary_df = summary_dfs[j]
        for (i, g) in enumerate(active_groups)
            use_frag, metric, label = g
            row = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, summary_df)
            val_matrix[j, i] = isempty(row) ? 0.0 : row[1, :success_rate]
            count_matrix[j, i] = isempty(row) ? 0 : row[1, :successful_runs]
        end
    end

    color_map = Dict(
        "Base, no similarity" => "#E69F00",
        "Base, Tanimoto" => "#D55E00",
        "w/ BRICS, no similarity" => "#56B4E9",
        "w/ BRICS, Tanimoto" => "#0072B2"
    )
    colors = reshape([color_map[g[3]] for g in active_groups], 1, N)
    labels = reshape([g[3] for g in active_groups], 1, N)
    
    if is_paper
        legend_pos = :outertop
        p_kwargs = (
            size = (1050, 950),
            legend_columns = N > 2 ? 2 : -1,
            title = "",
            tickfontsize = 20,
            guidefontsize = 22,
            legendfontsize = 20
        )
        annot_fontsize = 20
    else
        legend_pos = :outerbottom
        p_kwargs = (
            size = (600, 450),
            title = title_str,
            tickfontsize = 8,
            guidefontsize = 9,
            legendfontsize = 8
        )
        annot_fontsize = 10
    end
    
    annots = Matrix{Plots.PlotText}(undef, M, N)
    for j in 1:M
        for i in 1:N
            val = val_matrix[j, i]
            count = count_matrix[j, i]
            annots[j, i] = Plots.text("$(round(Int, val))%", annot_fontsize, :left, :bottom, rotation=60)
        end
    end
    
    x_ticks = ["$(datasets[j])\n($(num_problems[j]))" for j in 1:M]
    
    p = groupedbar(
        x_ticks,
        val_matrix,
        xlabel = "SynRXN Reaction Rebalancing Dataset",
        ylabel = "Reaction Completion Success Rate (%)",
        label = labels,
        legend = legend_pos,
        ylimits = (0, 100),
        xlims = (0.5 - 0.85 / 2, (M - 0.5) + 0.85 / 2),
        widen = false,
        series_annotations = annots,
        color = colors,
        linewidth = is_paper ? 3.5 : 1.2,
        bar_width = 0.85,
        edgecolor = :black;
        p_kwargs...
    )
    savefig(p, joinpath(output_dir, filename))
end

function plot_missing_molecule_success_rate_helper(
        output_dir::String, datasets::Vector{String}, num_problems::Vector{Int},
        success_rates_without::Vector{Float64}, success_rates_with::Vector{Float64},
        counts_without::Vector{Int}, counts_with::Vector{Int}, filename::String;
        is_paper::Bool = false,
        title_str::String = "Missing Molecule Synthesis (SynRXN rbl)\nSuccess Rate"
)
    M = length(datasets)
    val_matrix = [success_rates_without success_rates_with] # M x 2
    
    if is_paper
        p_kwargs = (
            size = (1050, 950),
            legend_columns = -1,
            title = "",
            tickfontsize = 20,
            guidefontsize = 22,
            legendfontsize = 20
        )
        annot_fontsize = 20
    else
        p_kwargs = (
            size = (600, 450),
            title = title_str,
            tickfontsize = 8,
            guidefontsize = 9,
            legendfontsize = 8
        )
        annot_fontsize = 10
    end
    
    annots = Matrix{Plots.PlotText}(undef, M, 2)
    for j in 1:M
        annots[j, 1] = Plots.text("$(round(success_rates_without[j]; digits=1))%", annot_fontsize, :center, :bottom)
        annots[j, 2] = Plots.text("$(round(success_rates_with[j]; digits=1))%", annot_fontsize, :center, :bottom)
    end
    
    x_ticks = ["$(datasets[j])\n($(num_problems[j]))" for j in 1:M]
    
    p = groupedbar(
        x_ticks,
        val_matrix,
        xlabel = "SynRXN Reaction Rebalancing Dataset",
        ylabel = "Missing Molecules Found Success Rate (%)",
        label = ["base grammar" "with BRICS fragments"],
        legend = :outertop,
        ylimits = (0, 100),
        xlims = (0.5 - 0.85 / 2, (M - 0.5) + 0.85 / 2),
        widen = false,
        series_annotations = annots,
        color = ["#E69F00" "#0072B2"],
        edgecolor = :black,
        bar_width = 0.85,
        linewidth = is_paper ? 3.5 : 1.2;
        p_kwargs...
    )
    savefig(p, joinpath(output_dir, filename))
end

function plot_comparison_boxplot_helper(
        output_dir::String, succ_rxn_dfs::Vector{DataFrame},
        metrics::Vector{String}, col_name::Symbol, ylabel_str::String, 
        filename::String, title_str::String, datasets::Vector{String};
        is_log_scale::Bool=false, val_digits::Int=1, val_suffix::String="",
        is_paper::Bool = false, y_cap::Union{Real, Nothing} = nothing,
        legend_position = nothing, legend_columns = nothing
)
    all_groups = [
        (false, "none", "Base, no similarity", "#E69F00"),
        (false, "tanimoto", "Base, Tanimoto", "#D55E00"),
        (true, "none", "w/ BRICS, no similarity", "#56B4E9"),
        (true, "tanimoto", "w/ BRICS, Tanimoto", "#0072B2")
    ]
    
    active_groups = filter(g -> g[2] in metrics, all_groups)
    N = length(active_groups)
    M = length(datasets)
    
    total_span = is_paper ? 0.95 : 0.44
    
    if is_paper
        legend_pos = legend_position !== nothing ? legend_position : :bottom
        p_kwargs = (
            size = (1050, 950),
            legend_columns = legend_columns !== nothing ? legend_columns : (N > 2 ? 2 : 1),
            title = "",
            tickfontsize = 20,
            guidefontsize = 22,
            legendfontsize = 20
        )
        annot_fontsize = 15
    else
        legend_pos = legend_position !== nothing ? legend_position : :outerbottom
        p_kwargs = (
            size = (600, 450),
            title = title_str,
            tickfontsize = 8,
            guidefontsize = 9,
            legendfontsize = 8
        )
        annot_fontsize = 7
    end
    
    p = boxplot(
        xticks = (1:M, datasets),
        xlims = (1 - total_span / 2, M + total_span / 2),
        widen = false,
        xlabel = "SynRXN Reaction Rebalancing Dataset",
        ylabel = ylabel_str,
        legend = legend_pos,
        yscale = is_log_scale ? :log10 : :identity;
        p_kwargs...
    )
    
    # Pre-populate legend entries to guarantee the correct legend order
    for g in active_groups
        boxplot!(p, Float64[], Float64[], label = g[3], color = g[4],
            seriestype = :shape, fillalpha = 0.7, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
    end
    
    max_val = 0.0
    
    for j in 1:M
        succ_rxn_df = succ_rxn_dfs[j]
        xs, bar_w = get_boxplot_positions_and_width(j, N; total_span = total_span)
        for (idx, g) in enumerate(active_groups)
            use_frag, metric, label, col = g
            sub_df = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, succ_rxn_df)
            vals = sub_df[!, col_name]
            if !isempty(vals)
                boxplot!(p, fill(xs[idx], length(vals)), vals,
                    label = "", color = col, fillalpha = 0.7,
                    bar_width = bar_w, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
                max_val = max(max_val, maximum(vals))
            end
        end
    end
    
    if max_val == 0.0
        annotate!(p, (M + 1) / 2, is_log_scale ? 1.0 : 0.5, text("No successful runs", 10, :center))
        savefig(p, joinpath(output_dir, filename))
        return
    end

    y_max_multiplier = is_paper ? (is_log_scale ? 10.0 : 1.45) : (is_log_scale ? 2.5 : 1.3)
    y_min = is_log_scale ? 0.5 : 0
    final_y_max = y_cap !== nothing ? y_cap : max_val * y_max_multiplier
    boxplot!(p, ylimits = (y_min, final_y_max))
    
    if !is_paper
        annot_y = is_log_scale ? max_val * 1.5 : max_val * 1.15
        for j in 1:M
            succ_rxn_df = succ_rxn_dfs[j]
            xs, _ = get_boxplot_positions_and_width(j, N; total_span = total_span)
            for (idx, g) in enumerate(active_groups)
                use_frag, metric, label, col = g
                sub_df = filter(r -> r.similarity_metric == metric && r.use_fragments == use_frag, succ_rxn_df)
                vals = sub_df[!, col_name]
                if !isempty(vals)
                    annotate!(p, xs[idx], annot_y,
                        text("Avg: $(round(mean(vals); digits=val_digits))$(val_suffix)", annot_fontsize, :center, :bottom))
                end
            end
        end
    end

    savefig(p, joinpath(output_dir, filename))
end

function plot_missing_molecule_boxplot_helper(
        output_dir::String,
        runtimes_or_mols_by_dataset::Vector,
        ylabel_str::String, filename::String, title_str::String,
        datasets::Vector{String}; val_digits::Int=1, val_suffix::String="",
        is_paper::Bool = false, y_cap::Union{Real, Nothing} = nothing
)
    M = length(datasets)
    N = 2
    
    total_span = is_paper ? 0.95 : 0.44
    
    if is_paper
        p_kwargs = (
            size = (1050, 950),
            legend_columns = -1,
            title = "",
            tickfontsize = 20,
            guidefontsize = 22,
            legendfontsize = 20
        )
        annot_fontsize = 18
    else
        p_kwargs = (
            size = (600, 450),
            title = title_str,
            tickfontsize = 8,
            guidefontsize = 9,
            legendfontsize = 8
        )
        annot_fontsize = 8
    end
    
    p = boxplot(
        xticks = (1:M, datasets),
        xlims = (1 - total_span / 2, M + total_span / 2),
        widen = false,
        xlabel = "SynRXN Reaction Rebalancing Dataset",
        ylabel = ylabel_str,
        legend = :outertop;
        p_kwargs...
    )
    boxplot!(p, Float64[], Float64[], label = "base grammar", color = "#E69F00",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
    boxplot!(p, Float64[], Float64[], label = "with BRICS fragments", color = "#0072B2",
        seriestype = :shape, fillalpha = 0.7, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
        
    max_val = 0.0
    
    for j in 1:M
        data_dict = runtimes_or_mols_by_dataset[j]
        xs, bar_w = get_boxplot_positions_and_width(j, N; total_span = total_span)
        
        vals_without = data_dict[false]
        vals_with = data_dict[true]
        
        if !isempty(vals_without)
            boxplot!(p, fill(xs[1], length(vals_without)), vals_without,
                label = "", color = "#E69F00", fillalpha = 0.7, bar_width = bar_w, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
            max_val = max(max_val, maximum(vals_without))
        end
        if !isempty(vals_with)
            boxplot!(p, fill(xs[2], length(vals_with)), vals_with,
                label = "", color = "#0072B2", fillalpha = 0.7, bar_width = bar_w, linecolor = :black, linewidth = is_paper ? 3.5 : 1.2)
            max_val = max(max_val, maximum(vals_with))
        end
    end
    
    max_val = max_val == 0.0 ? 1.0 : max_val
    max_val_multiplier = is_paper ? 1.4 : 1.3
    final_y_max = y_cap !== nothing ? y_cap : max_val * max_val_multiplier
    boxplot!(p, ylimits = (0, final_y_max))
    
    if !is_paper
        for j in 1:M
            data_dict = runtimes_or_mols_by_dataset[j]
            xs, _ = get_boxplot_positions_and_width(j, N; total_span = total_span)
            
            vals_without = data_dict[false]
            vals_with = data_dict[true]
            
            if !isempty(vals_without)
                annotate!(p, xs[1], max_val * 1.15,
                    text("Avg: $(round(mean(vals_without); digits=val_digits))$(val_suffix)", annot_fontsize, :center, :bottom))
            else
                annotate!(p, xs[1], max_val * 1.15, text("N/A", annot_fontsize, :center, :bottom))
            end
            
            if !isempty(vals_with)
                annotate!(p, xs[2], max_val * 1.15,
                    text("Avg: $(round(mean(vals_with); digits=val_digits))$(val_suffix)", annot_fontsize, :center, :bottom))
            else
                annotate!(p, xs[2], max_val * 1.15, text("N/A", annot_fontsize, :center, :bottom))
            end
        end
    end
    
    savefig(p, joinpath(output_dir, filename))
end

# -------------------------------------------------------------
# Data Processing Helpers
# -------------------------------------------------------------

function compute_summary_df(results_df::DataFrame)
    return combine(groupby(results_df, [:similarity_metric, :use_fragments])) do sdf
        eval_runs = nrow(sdf)
        success_runs = sum(sdf.success)
        success_rate = eval_runs > 0 ? (success_runs / eval_runs) * 100 : 0.0

        # Calculate average time only for successful runs
        success_times = sdf.elapsed_time_seconds[sdf.success]
        avg_time = isempty(success_times) ? 0.0 : mean(success_times)

        total_time = sum(sdf.elapsed_time_seconds)

        return (
            evaluated_runs = eval_runs,
            successful_runs = success_runs,
            success_rate = round(success_rate; digits = 1),
            avg_synthesis_time_seconds = round(avg_time; digits = 2),
            total_synthesis_time_seconds = round(total_time; digits = 2)
        )
    end
end

function compute_missing_molecule_success_stats(missing_molecule_synthesis_stats::Dict{Bool, Dict{Symbol, Any}})
    synthesis_eval_limit = length(missing_molecule_synthesis_stats[true][:successes])
    
    succ_count_without = sum(missing_molecule_synthesis_stats[false][:successes])
    succ_count_with = sum(missing_molecule_synthesis_stats[true][:successes])
    
    success_without = succ_count_without / max(1, length(missing_molecule_synthesis_stats[false][:successes])) * 100
    success_with = succ_count_with / max(1, length(missing_molecule_synthesis_stats[true][:successes])) * 100
    
    return synthesis_eval_limit, succ_count_without, succ_count_with, success_without, success_with
end

function compute_unsynthesised_sizes_df(dataset::String, missing_molecule_synthesis_stats::Dict{Bool, Dict{Symbol, Any}})
    sizes_with = missing_molecule_synthesis_stats[true][:unfeasible_sizes]
    sizes_without = missing_molecule_synthesis_stats[false][:unfeasible_sizes]

    total_with = length(sizes_with)
    total_without = length(sizes_without)

    pct_le_3_with = total_with > 0 ? count(s -> s <= 3, sizes_with) / total_with * 100 : 0.0
    pct_4_6_with = total_with > 0 ? count(s -> 4 <= s <= 6, sizes_with) / total_with * 100 : 0.0
    pct_7_9_with = total_with > 0 ? count(s -> 7 <= s <= 9, sizes_with) / total_with * 100 : 0.0
    pct_ge_10_with = total_with > 0 ? count(s -> s >= 10, sizes_with) / total_with * 100 : 0.0

    pct_le_3_without = total_without > 0 ? count(s -> s <= 3, sizes_without) / total_without * 100 : 0.0
    pct_4_6_without = total_without > 0 ? count(s -> 4 <= s <= 6, sizes_without) / total_without * 100 : 0.0
    pct_7_9_without = total_without > 0 ? count(s -> 7 <= s <= 9, sizes_without) / total_without * 100 : 0.0
    pct_ge_10_without = total_without > 0 ? count(s -> s >= 10, sizes_without) / total_without * 100 : 0.0

    return DataFrame(
        Symbol("Dataset") => [dataset, dataset],
        Symbol("Methodology") => ["Base Molecule Grammar", "with BRICS Fragments"],
        Symbol("TotalUnsynthesised") => [total_without, total_with],
        Symbol("SizeLE3") => [round(pct_le_3_without; digits=1), round(pct_le_3_with; digits=1)],
        Symbol("Size4to6") => [round(pct_4_6_without; digits=1), round(pct_4_6_with; digits=1)],
        Symbol("Size7to9") => [round(pct_7_9_without; digits=1), round(pct_7_9_with; digits=1)],
        Symbol("SizeGE10") => [round(pct_ge_10_without; digits=1), round(pct_ge_10_with; digits=1)]
    )
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
    summary_df = compute_summary_df(results_df)
    CSV.write(joinpath(run_dir, "summary.csv"), summary_df)

    # 3. Bar plot & CSV: Success rate of missing molecule synthesis subproblem
    synthesis_eval_limit, succ_count_without, succ_count_with, success_without, success_with = 
        compute_missing_molecule_success_stats(missing_molecule_synthesis_stats)

    success_rate_df = DataFrame(
        Use_Fragments = ["base grammar", "with BRICS fragments"],
        Success_Rate_Percentage = [success_without, success_with],
        Successful_Problems = [succ_count_without, succ_count_with]
    )
    CSV.write(joinpath(run_dir, "missing_molecule_synthesis_success_rate.csv"), success_rate_df)

    plot_missing_molecule_success_rate_helper(
        run_dir, [dataset], [synthesis_eval_limit],
        [success_without], [success_with],
        [succ_count_without], [succ_count_with],
        "missing_molecule_synthesis_success_rate.svg";
        is_paper = false
    )

    # 4. Sizes of untractable molecules summary CSV
    unsynthesised_sizes_df = compute_unsynthesised_sizes_df(dataset, missing_molecule_synthesis_stats)
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

    runtimes_dict = Dict(false => runtimes_succ_without, true => runtimes_succ_with)
    plot_missing_molecule_boxplot_helper(
        run_dir, [runtimes_dict], "Runtime (s)", "missing_molecule_synthesis_average_runtime.svg",
        "Missing Molecule Synthesis (SynRXN rbl)\nRuntime of Successful Molecule Syntheses",
        [dataset]; val_digits=2, val_suffix="s", is_paper=false
    )

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

    mols_dict = Dict(false => mols_succ_without, true => mols_succ_with)
    plot_missing_molecule_boxplot_helper(
        run_dir, [mols_dict], "Molecules Synthesised", "average_molecules_synthesized.svg",
        "Missing Molecule Synthesis (SynRXN rbl)\nNumber of Molecules Synthesised until All Targets Found for Successful Problems",
        [dataset]; val_digits=1, val_suffix="", is_paper=false
    )

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
        run_dir, [succ_rxn_df], metrics_order, :reactions_synthesized, "Reactions Synthesised",
        "average_reactions_synthesized.svg",
        "Reaction Rebalancing (SynRXN rbl)\nNumber of Reactions Synthesised Until Target Found for Successful Problems",
        [dataset];
        is_log_scale=true, val_digits=1, val_suffix="", is_paper=false
    )

    # 8. Success Rate Comparison Plots
    # All metrics combined
    plot_success_rate_helper(
        run_dir, [dataset], [synthesis_eval_limit], [summary_df], metrics_order, "synthesis_success_rate.svg";
        is_paper=false
    )
    
    # 9. Runtime Comparison Plots
    # All metrics combined
    plot_comparison_boxplot_helper(
        run_dir, [succ_rxn_df], metrics_order, :elapsed_time_seconds, "Synthesis Time (s)",
        "synthesis_average_time.svg",
        "Reaction Rebalancing (SynRXN rbl)\nRuntime of Successful Reaction Syntheses",
        [dataset];
        is_log_scale=false, val_digits=2, val_suffix="s", is_paper=false
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
