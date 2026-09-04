using CRNSynthesizer
using DataStructures
using BenchmarkTools
using StatsPlots

ENV["GKSwstype"] = "100"

function run_benchmarks()
    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2, "[C]" => 4, "[N]" => 3)

    depths = 6:9
    runtimes = Float64[]

    for max_depth in depths
        settings = SynthesizerSettings(
            max_depth = max_depth, options = Dict{Symbol, Any}(:unique_candidates => true))
        println("Benchmarking max_depth = $max_depth...")
        
        # Benchmark with a limit to avoid excessive runtimes on large depths
        b = @benchmark synthesize_molecules($atom_valences; settings = $settings) samples=3 seconds=10 evals=1
        t = minimum(b).time / 1e9 # convert nanoseconds to seconds
        push!(runtimes, t)
        println("Runtime for max_depth = $max_depth: ", t, " seconds")
    end

    p = bar(depths, runtimes,
        yscale = :log10,
        xlabel = "Maximum Depth",
        ylabel = "Runtime (seconds)",
        title = "Time Required to Exhaust a Molecule Synthesiser at Increasing Depth",
        legend = false,
        color = :dodgerblue,
        linewidth = 1,
        edgecolor = :black,
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        size = (600, 450),
        series_annotations = map(t -> Plots.text("$(round(t, sigdigits=3))s", 8, :center, :bottom), runtimes)
    )

    output_path = joinpath(@__DIR__, "molecule_tractability_runtime.svg")
    savefig(p, output_path)
    println("Plot saved to ", output_path)
end

run_benchmarks()