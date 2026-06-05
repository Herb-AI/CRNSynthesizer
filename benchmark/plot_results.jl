using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

if length(ARGS) != 1
    println("Usage: julia --project=. benchmark/plot_results.jl <run_dir>")
    println("Example: julia --project=. benchmark/plot_results.jl benchmark/results/2026-05-30_12-34-56")
    exit(1)
end

run_dir = ARGS[1]
if !isdir(run_dir)
    println("Error: Directory '$run_dir' does not exist.")
    exit(1)
end

if !isfile(joinpath(run_dir, "benchmark_data.jls"))
    println("Error: No 'benchmark_data.jls' found in '$run_dir'.")
    exit(1)
end

include("log_utils.jl")

println("Re-generating plots and CSVs for run: $run_dir...")
generate_plots(run_dir)
