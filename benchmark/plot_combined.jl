using DataFrames
using Serialization

if length(ARGS) < 4
    println("Usage: julia --project=. benchmark/plot_combined.jl <dir1> <dir2> <dir3> <dir4> [output_dir]")
    println("Example: julia --project=. benchmark/plot_combined.jl benchmark/results/dir1 benchmark/results/dir2 benchmark/results/dir3 benchmark/results/dir4 benchmark/results/combined")
    exit(1)
end

input_dirs = ARGS[1:4]
output_dir = length(ARGS) >= 5 ? ARGS[5] : joinpath(@__DIR__, "results", "combined")
mkpath(output_dir)

include("log_utils.jl")

loaded_data = Dict{String, Tuple{DataFrame, String, Dict{Bool, Dict{Symbol, Any}}}}()

for dir in input_dirs
    if !isdir(dir)
        println("Error: Directory '$dir' does not exist.")
        exit(1)
    end
    jls_path = joinpath(dir, "benchmark_data.jls")
    if !isfile(jls_path)
        println("Error: File '$jls_path' not found.")
        exit(1)
    end
    data = Serialization.deserialize(jls_path)
    # data is (results_df, dataset, missing_molecule_synthesis_stats)
    dataset_name = data[2]
    loaded_data[dataset_name] = data
end

# Check for the expected datasets
expected_datasets = ["complex", "mbs", "mnc", "mos"]
all_loaded = keys(loaded_data)
datasets_to_plot = [d for d in expected_datasets if d in all_loaded]
other_datasets = sort([d for d in all_loaded if !(d in datasets_to_plot)])
append!(datasets_to_plot, other_datasets)

if length(datasets_to_plot) < 4
    println("Warning: Only found datasets: $(join(datasets_to_plot, ", ")). Expected at least 4.")
end

# Organize the data in the ordered datasets_to_plot
results_dfs = DataFrame[]
missing_stats = Dict{Bool, Dict{Symbol, Any}}[]
num_problems = Int[]

success_rates_without = Float64[]
success_rates_with = Float64[]
counts_without = Int[]
counts_with = Int[]

runtimes_by_dataset = Dict{Bool, Vector{Float64}}[]
mols_by_dataset = Dict{Bool, Vector{Int}}[]

for d in datasets_to_plot
    rdf, _, m_stats = loaded_data[d]
    push!(results_dfs, rdf)
    push!(missing_stats, m_stats)
    
    # Missing molecule calculations
    synthesis_eval_limit, succ_count_without, succ_count_with, success_without, success_with = 
        compute_missing_molecule_success_stats(m_stats)
    
    push!(num_problems, synthesis_eval_limit)
    push!(success_rates_without, success_without)
    push!(success_rates_with, success_with)
    push!(counts_without, succ_count_without)
    push!(counts_with, succ_count_with)
    
    # Runtimes dict
    runtimes_succ_with = m_stats[true][:runtimes][m_stats[true][:successes]]
    runtimes_succ_without = m_stats[false][:runtimes][m_stats[false][:successes]]
    push!(runtimes_by_dataset, Dict(false => runtimes_succ_without, true => runtimes_succ_with))
    
    # Molecules dict
    mols_succ_with = m_stats[true][:molecules_synthesized][m_stats[true][:successes]]
    mols_succ_without = m_stats[false][:molecules_synthesized][m_stats[false][:successes]]
    push!(mols_by_dataset, Dict(false => mols_succ_without, true => mols_succ_with))
end

# Compute summary DataFrames for reaction rebalancing plots
summary_dfs = DataFrame[]
succ_rxn_dfs = DataFrame[]
for results_df in results_dfs
    push!(summary_dfs, compute_summary_df(results_df))
    push!(succ_rxn_dfs, filter(r -> r.success == true, results_df))
end

metrics_order = ["none", "tanimoto"]

println("Generating combined plots in: $output_dir")

# 1. Missing molecule success rate
plot_missing_molecule_success_rate_helper(
    output_dir, datasets_to_plot, num_problems,
    success_rates_without, success_rates_with,
    counts_without, counts_with,
    "molecules_success_rate.pdf";
    is_paper = true
)

# 2. Missing molecule average runtime
plot_missing_molecule_boxplot_helper(
    output_dir, runtimes_by_dataset, "Runtime Until Missing Molecules Found (s)",
    "molecules_runtime.pdf",
    "",
    datasets_to_plot; val_digits=2, val_suffix="s", is_paper=true, y_cap=4.0
)

# 3. Average molecules synthesized
plot_missing_molecule_boxplot_helper(
    output_dir, mols_by_dataset, "Molecules Synthesised Until Targets Found",
    "boxplot_molecules_synthesized.pdf",
    "",
    datasets_to_plot; val_digits=1, val_suffix="", is_paper=true, y_cap=250.0
)

# 4. Average reactions synthesized
plot_comparison_boxplot_helper(
    output_dir, succ_rxn_dfs, metrics_order, :reactions_synthesized, "Reactions Synthesised Until Target Found",
    "boxplot_reactions_synthesized.pdf",
    "",
    datasets_to_plot; is_log_scale=true, val_digits=1, val_suffix="", is_paper=true, y_cap=30000.0
)

# 5. Synthesis success rate
plot_success_rate_helper(
    output_dir, datasets_to_plot, num_problems, summary_dfs, metrics_order,
    "reaction_success_rate.pdf";
    is_paper=true
)

# 6. Synthesis average time
plot_comparison_boxplot_helper(
    output_dir, succ_rxn_dfs, metrics_order, :elapsed_time_seconds, "Runtime Until Target Reaction Found (s)",
    "reaction_average_time.pdf",
    "",
    datasets_to_plot; is_log_scale=false, val_digits=2, val_suffix="s", is_paper=true, y_cap=50.0,
    legend_position = :outertop
)

# 7. Combined missing molecule synthesis unsynthesised sizes
combined_unsynthesised_sizes_df = DataFrame()
for (idx, d) in enumerate(datasets_to_plot)
    m_stats = missing_stats[idx]
    df = compute_unsynthesised_sizes_df(d, m_stats)
    append!(combined_unsynthesised_sizes_df, df)
end
CSV.write(joinpath(output_dir, "molecules_unsynthesised_sizes.csv"), combined_unsynthesised_sizes_df)

println("✓ Combined plots generated successfully!")
