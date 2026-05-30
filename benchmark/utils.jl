using CRNSynthesizer
using DataStructures
using LRUCache
import MoleculeFlow
import TOML
using Dates
using DataFrames
using CSV
using StatsPlots
using Statistics
using Serialization

# -------------------------------------------------------------
# Caching wrappers to avoid expensive RDKit/MoleculeFlow calls
# -------------------------------------------------------------
const CACHE_SIZE = 1024

const MOL_CACHE = LRU{String, Union{MoleculeFlow.Molecule, Missing}}(maxsize = CACHE_SIZE)
const BRICS_CACHE = LRU{String, Union{Vector{String}, Missing}}(maxsize = CACHE_SIZE)
const CLEANED_FRAG_CACHE = LRU{String, Union{MoleculeFlow.Molecule, Missing}}(maxsize = CACHE_SIZE)
const PARSED_MOLECULE_CACHE = LRU{String, Molecule}(maxsize = CACHE_SIZE)
const FRAG_ATOM_COUNT_CACHE = LRU{String, Int}(maxsize = CACHE_SIZE)
const MOL_ATOM_COUNT_CACHE = LRU{String, Int}(maxsize = CACHE_SIZE)

function get_mol_cached(smiles::String)::Union{MoleculeFlow.Molecule, Missing}
    get!(MOL_CACHE, smiles) do
        MoleculeFlow.mol_from_smiles(smiles)
    end
end

function get_brics_cached(smiles::String)::Union{Vector{String}, Missing}
    get!(BRICS_CACHE, smiles) do
        mol = get_mol_cached(smiles)
        ismissing(mol) ? missing : MoleculeFlow.brics_decompose(mol; min_fragment_size = 2)
    end
end

function get_cleaned_frag_cached(frag_smiles::String)::Union{MoleculeFlow.Molecule, Missing}
    get!(CLEANED_FRAG_CACHE, frag_smiles) do
        frag_m = get_mol_cached(frag_smiles)
        ismissing(frag_m) ? missing : MoleculeFlow.delete_substructs(frag_m, "[#0]")
    end
end

function get_frag_atom_count_cached(frag_smiles::String)::Int
    get!(FRAG_ATOM_COUNT_CACHE, frag_smiles) do
        cleaned_m = get_cleaned_frag_cached(frag_smiles)
        ismissing(cleaned_m) ? 0 : length(collect(MoleculeFlow.get_atoms(cleaned_m)))
    end
end

function get_mol_atom_count_cached(smiles::String)::Int
    get!(MOL_ATOM_COUNT_CACHE, smiles) do
        mol = get_mol_cached(smiles)
        ismissing(mol) ? 0 : length(collect(MoleculeFlow.get_atoms(mol)))
    end
end

function get_parsed_molecule_cached(smiles::String)::Molecule
    get!(PARSED_MOLECULE_CACHE, smiles) do
        from_SMILES(make_smiles_custom_explicit(make_smiles_rdkit_explicit(smiles)))
    end
end

function parse_rxn_cached(rxn_str::SubString{String})::Tuple{
        Vector{Molecule}, Vector{Molecule}}
    # Pre-process to convert free hydrogen/oxygen representations to stable molecules
    # We use lookarounds to ensure we only replace isolated atoms, making an exception
    # if the bond is already explicit (e.g. [H][H] or [H]-[H] or [O]=[O]).
    processed_rxn = replace(rxn_str, r"(?<=^|\.|>>)(\[H\])\.(\[H\])(?=$|\.|>>)" => s"\1-\2")
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)(\[O\])\.(\[O\])(?=$|\.|>>)" =>
        s"\1=\2")

    # Replace any leftover standalone [H] and [O]
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)\[H\](?=$|\.|>>)" => "[H]-[H]")
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)\[O\](?=$|\.|>>)" => "[O]=[O]")

    rxn_parts = split(processed_rxn, ">>")
    if length(rxn_parts) != 2
        error("Invalid reaction string: $rxn_str")
    end
    reactants_part = strip(rxn_parts[1])
    products_part = strip(rxn_parts[2])

    reactants_smiles = filter(!isempty, map(strip, split(reactants_part, ".")))
    products_smiles = filter(!isempty, map(strip, split(products_part, ".")))

    reactants = [get_parsed_molecule_cached(String(s)) for s in reactants_smiles]
    products = [get_parsed_molecule_cached(String(s)) for s in products_smiles]

    return reactants, products
end

function parse_syn_problem(reaction_str::AbstractString)::ProblemDefinition
    parts = split(reaction_str, ",")
    if length(parts) != 2
        error("Invalid input format: expected two reaction strings separated by a comma")
    end
    incomplete_rxn_str = strip(parts[1])
    target_rxn_str = strip(parts[2])

    incomplete_reactants, incomplete_products = parse_rxn_cached(incomplete_rxn_str)
    target_reactants, target_products = parse_rxn_cached(target_rxn_str)

    # Collect unique molecules from the target reaction
    all_molecules = Molecule[]
    for m in vcat(target_reactants, target_products)
        if !(m in all_molecules)
            push!(all_molecules, m)
        end
    end

    # Construct the target reaction with proper stoichiometry
    reactant_counts = OrderedDict{Molecule, Int}()
    for m in target_reactants
        reactant_counts[m] = get(reactant_counts, m, 0) + 1
    end
    inputs = [(reactant_counts[m], m) for m in unique(target_reactants)]

    product_counts = OrderedDict{Molecule, Int}()
    for m in target_products
        product_counts[m] = get(product_counts, m, 0) + 1
    end
    outputs = [(product_counts[m], m) for m in unique(target_products)]

    reaction = CRNSynthesizer.Reaction(nothing, inputs, outputs, false)
    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Calculate known indices
    incomplete_molecules = OrderedSet{Molecule}(vcat(incomplete_reactants, incomplete_products))
    selected_known_indices = Int[]
    for (i, m) in enumerate(all_molecules)
        if m in incomplete_molecules
            push!(selected_known_indices, i)
        end
    end

    known_molecules = all_molecules[selected_known_indices]
    goal_molecules = Molecule[]
    for i in eachindex(all_molecules)
        if !(i in selected_known_indices)
            push!(goal_molecules, all_molecules[i])
        end
    end

    atom_valences = get_valences_from_molecules(all_molecules)

    inc_reactant_counts = OrderedDict{Molecule, Int}()
    for m in incomplete_reactants
        inc_reactant_counts[m] = get(inc_reactant_counts, m, 0) + 1
    end
    partial_inputs = [(inc_reactant_counts[m], m) for m in unique(incomplete_reactants)]

    inc_product_counts = OrderedDict{Molecule, Int}()
    for m in incomplete_products
        inc_product_counts[m] = get(inc_product_counts, m, 0) + 1
    end
    partial_outputs = [(inc_product_counts[m], m) for m in unique(incomplete_products)]

    partial_reaction = if isempty(partial_inputs) && isempty(partial_outputs)
        nothing
    else
        CRNSynthesizer.Reaction(nothing, partial_inputs, partial_outputs, false)
    end

    return ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network;
        partial_reaction = partial_reaction
    )
end

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

"""
    save_benchmark_results(
        results_df::DataFrame,
        dataset::String,
        first_filter_stats::Dict{Symbol, Any},
        second_filter_stats::Dict{Bool, Dict{Symbol, Any}}
    )

Save feasibility benchmark results (metadata, raw csv, summary csv, and advanced comparison plots/CSVs).
"""
function save_benchmark_results(
        run_dir::String,
        results_df::DataFrame,
        dataset::String,
        first_filter_stats::Dict{Symbol, Any},
        second_filter_stats::Dict{Bool, Dict{Symbol, Any}}
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
        println(io,
            "preprocessing: Feasibility classification (molecule size <= 6 or sharing BRICS fragments with known molecules; min_fragment_size = 2)")
        println(io, "evaluation_script_commit: $(get_git_commit())")
    end

    # Write detailed results CSV
    CSV.write(joinpath(run_dir, "results.csv"), results_df)

    # Serialize raw data for standalone plot generation
    Serialization.serialize(joinpath(run_dir, "benchmark_data.jls"), (results_df, dataset, first_filter_stats, second_filter_stats))

    generate_plots(run_dir, results_df, dataset, first_filter_stats, second_filter_stats)
end

function generate_plots(run_dir::String)
    data = Serialization.deserialize(joinpath(run_dir, "benchmark_data.jls"))
    generate_plots(run_dir, data...)
end

function generate_plots(
        run_dir::String,
        results_df::DataFrame,
        dataset::String,
        first_filter_stats::Dict{Symbol, Any},
        second_filter_stats::Dict{Bool, Dict{Symbol, Any}}
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

    # Helpers for plotting and label formatting
    function create_pie_plot(labels_raw::Vector{String}, counts::Vector{Int}, title_str::String)
        total = sum(counts)
        if total == 0
            return pie(title = title_str, legend = false,
                grid = false, xaxis = false, yaxis = false,
                annotation = (0.5, 0.5, "No data available"))
        end
        labels_with_pct = ["$lbl: $val ($(round(val / total * 100; digits=1))%)"
                           for (lbl, val) in zip(labels_raw, counts)]
        return pie(labels_with_pct, counts,
            title = title_str,
            legend = :topright,
            titlefontsize = 10,
            legendfontsize = 8,
            size = (600, 450)
        )
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

    # 1. Pie plot & CSV: Feasible vs Infeasible (first filter)
    feasible_count = first_filter_stats[:feasible_count]
    infeasible_count = first_filter_stats[:infeasible_count]
    total_problems = feasible_count + infeasible_count

    first_filter_feas_df = DataFrame(
        Status = ["Feasible", "Infeasible"],
        Count = [feasible_count, infeasible_count],
        Percentage = total_problems > 0 ?
                     [
            feasible_count / total_problems * 100, infeasible_count / total_problems * 100] :
                     [0.0, 0.0]
    )
    CSV.write(joinpath(run_dir, "first_filter_feasibility.csv"), first_filter_feas_df)

    p_first_feas = create_pie_plot(
        ["Feasible", "Infeasible"],
        [feasible_count, infeasible_count],
        "BRICS and Size <= 6 Filter (SynRXN rbl/$dataset) - $total_problems Problems Analyzed\nFilter Outcome Distribution"
    )
    savefig(p_first_feas, joinpath(run_dir, "first_filter_feasibility.pdf"))

    # 2. Pie plot & CSV: Distribution of unfeasible molecule sizes (first filter)
    first_sizes_lbls, first_sizes_vals,
    first_sizes_df = make_size_dist_df(first_filter_stats[:unfeasible_molecule_sizes])
    CSV.write(joinpath(run_dir, "first_filter_unfeasible_sizes.csv"), first_sizes_df)

    p_first_sizes = create_pie_plot(
        first_sizes_lbls,
        first_sizes_vals,
        "BRICS and Size <= 6 Filter (SynRXN rbl/$dataset) - $infeasible_count Problems Analyzed\nDistribution of Unfeasible Molecule Sizes"
    )
    savefig(p_first_sizes, joinpath(run_dir, "first_filter_unfeasible_sizes.pdf"))

    # 3. Bar plot & CSV: Success rate of Atoms -> Molecules filter
    success_with = sum(second_filter_stats[true][:successes]) /
                   max(1, length(second_filter_stats[true][:successes])) * 100
    success_without = sum(second_filter_stats[false][:successes]) /
                      max(1, length(second_filter_stats[false][:successes])) * 100
    synthesis_eval_limit = length(second_filter_stats[true][:successes])

    succ_count_with = sum(second_filter_stats[true][:successes])
    succ_count_without = sum(second_filter_stats[false][:successes])

    success_rate_df = DataFrame(
        Use_Fragments = ["Without Fragments", "With Fragments"],
        Success_Rate_Percentage = [success_without, success_with],
        Successful_Problems = [succ_count_without, succ_count_with]
    )
    CSV.write(joinpath(run_dir, "second_filter_success_rate.csv"), success_rate_df)

    p_success_rate = bar(
        ["Without Fragments", "With Fragments"],
        [success_without, success_with],
        xlabel = "Fragments Setting",
        ylabel = "Success Rate (%)",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset)\nSuccess Rate of Missing Molecule Synthesis ($synthesis_eval_limit Problems Analyzed)",
        legend = false,
        ylimits = (0, 100),
        series_annotations = map(
            (v, c) -> Plots.text("$(round(v; digits=1))% ($c)", 8, :center, :bottom),
            [success_without, success_with], [succ_count_without, succ_count_with]),
        color = [:lightgrey, :dodgerblue],
        edgecolor = :black,
        linewidth = 1.2,
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        size = (600, 450)
    )
    savefig(p_success_rate, joinpath(run_dir, "second_filter_success_rate.pdf"))

    # 4. Pie plots & CSV: Sizes of unfeasible molecules based on second filter (with and without fragments)
    # With fragments
    sizes_with_lbls, sizes_with_vals,
    sizes_with_df = make_size_dist_df(second_filter_stats[true][:unfeasible_sizes])
    CSV.write(joinpath(run_dir, "second_filter_unfeasible_sizes_with_fragments.csv"), sizes_with_df)

    failed_with_count = sum(.!second_filter_stats[true][:successes])
    p_sizes_with = create_pie_plot(
        sizes_with_lbls,
        sizes_with_vals,
        "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $failed_with_count Problems Analyzed\nDistribution of Unfeasible Molecule Sizes (With Fragments)"
    )
    savefig(p_sizes_with, joinpath(run_dir, "second_filter_unfeasible_sizes_with_fragments.pdf"))

    # Without fragments
    sizes_without_lbls, sizes_without_vals,
    sizes_without_df = make_size_dist_df(second_filter_stats[false][:unfeasible_sizes])
    CSV.write(joinpath(run_dir, "second_filter_unfeasible_sizes_without_fragments.csv"), sizes_without_df)

    failed_without_count = sum(.!second_filter_stats[false][:successes])
    p_sizes_without = create_pie_plot(
        sizes_without_lbls,
        sizes_without_vals,
        "Missing Molecule Synthesis (SynRXN rbl/$dataset) - $failed_without_count Problems Analyzed\nDistribution of Unfeasible Molecule Sizes (Without Fragments)"
    )
    savefig(p_sizes_without, joinpath(run_dir, "second_filter_unfeasible_sizes_without_fragments.pdf"))

    # 5. Box plot & CSV: Runtime of Atoms -> Molecules synthesis
    runtimes_succ_with = second_filter_stats[true][:runtimes][second_filter_stats[true][:successes]]
    runtimes_succ_without = second_filter_stats[false][:runtimes][second_filter_stats[false][:successes]]

    runtime_with = isempty(runtimes_succ_with) ? 0.0 : mean(runtimes_succ_with)
    runtime_without = isempty(runtimes_succ_without) ? 0.0 : mean(runtimes_succ_without)

    runtime_df = DataFrame(
        Use_Fragments = ["Without Fragments", "With Fragments"],
        Average_Runtime_Seconds = [runtime_without, runtime_with]
    )
    CSV.write(joinpath(run_dir, "second_filter_average_runtime.csv"), runtime_df)

    p_runtime = boxplot(
        xticks = (1:2, ["Without Fragments", "With Fragments"]),
        xlims = (0.5, 2.5),
        xlabel = "Fragments Setting",
        ylabel = "Runtime (s)",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset)\nRuntime of Successful Missing Molecule Syntheses",
        legend = false,
        size = (600, 450),
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8
    )
    if !isempty(runtimes_succ_without)
        boxplot!(p_runtime, fill(1, length(runtimes_succ_without)), runtimes_succ_without, color=:lightgrey, fillalpha=0.7, bar_width=0.6, linecolor=:black)
    end
    if !isempty(runtimes_succ_with)
        boxplot!(p_runtime, fill(2, length(runtimes_succ_with)), runtimes_succ_with, color=:dodgerblue, fillalpha=0.7, bar_width=0.6, linecolor=:black)
    end

    max_rt_all = max(isempty(runtimes_succ_without) ? 0.0 : maximum(runtimes_succ_without), isempty(runtimes_succ_with) ? 0.0 : maximum(runtimes_succ_with))
    max_rt_all = max_rt_all == 0.0 ? 1.0 : max_rt_all
    boxplot!(p_runtime, ylims=(0, max_rt_all * 1.3))

    annotate!(p_runtime, 1, max_rt_all * 1.15, isempty(runtimes_succ_without) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(runtime_without; digits=2))s", 8, :center, :bottom))
    annotate!(p_runtime, 2, max_rt_all * 1.15, isempty(runtimes_succ_with) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(runtime_with; digits=2))s", 8, :center, :bottom))
    
    savefig(p_runtime, joinpath(run_dir, "second_filter_average_runtime.pdf"))

    # 6. Box plot & CSV: Number of molecules synthesized (with and without fragments)
    # Using the results from the second filter:
    mols_succ_with = second_filter_stats[true][:molecules_synthesized][second_filter_stats[true][:successes]]
    mols_succ_without = second_filter_stats[false][:molecules_synthesized][second_filter_stats[false][:successes]]

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
        xlabel = "Fragments Setting",
        ylabel = "Molecules Synthesized",
        title = "Missing Molecule Synthesis (SynRXN rbl/$dataset)\nNumber of Molecules Synthesized for Successful Problems",
        legend = false,
        size = (600, 450),
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8
    )
    if !isempty(mols_succ_without)
        boxplot!(p_mols, fill(1, length(mols_succ_without)), mols_succ_without, color=:lightgrey, fillalpha=0.7, bar_width=0.6, linecolor=:black)
    end
    if !isempty(mols_succ_with)
        boxplot!(p_mols, fill(2, length(mols_succ_with)), mols_succ_with, color=:dodgerblue, fillalpha=0.7, bar_width=0.6, linecolor=:black)
    end

    max_mols_all = max(isempty(mols_succ_without) ? 0.0 : maximum(mols_succ_without), isempty(mols_succ_with) ? 0.0 : maximum(mols_succ_with))
    max_mols_all = max_mols_all == 0.0 ? 1.0 : max_mols_all
    boxplot!(p_mols, ylims=(0, max_mols_all * 1.3))

    annotate!(p_mols, 1, max_mols_all * 1.15, isempty(mols_succ_without) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(avg_mols_without; digits=1))", 8, :center, :bottom))
    annotate!(p_mols, 2, max_mols_all * 1.15, isempty(mols_succ_with) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(avg_mols_with; digits=1))", 8, :center, :bottom))

    savefig(p_mols, joinpath(run_dir, "average_molecules_synthesized.pdf"))

    # 7. Box plots & CSV: Reactions synthesized per metric
    metrics_order = ["none", "simpson", "tanimoto", "both"]
    succ_rxn_df = filter(r -> r.success == true, results_df)

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
        title = "Reaction Rebalancing (SynRXN rbl/$dataset)\nNumber of Reactions Synthesized for Successful Problems\n(With Fragments)",
        legend = false,
        size = (600, 450),
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        yscale = :log10
    )
    max_rxn_with_val = 0.0
    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        if !isempty(vals)
            boxplot!(p_rxns_with, fill(i, length(vals)), vals, color=:dodgerblue, fillalpha=0.7, bar_width=0.6, linecolor=:black)
            max_rxn_with_val = max(max_rxn_with_val, maximum(vals))
        end
    end
    max_rxn_with_val = max_rxn_with_val == 0.0 ? 1.0 : max_rxn_with_val
    boxplot!(p_rxns_with, ylims=(0.5, max_rxn_with_val * 2.5))
    
    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        annot = isempty(vals) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(mean(vals); digits=1))", 8, :center, :bottom)
        annotate!(p_rxns_with, i, max_rxn_with_val * 1.5, annot)
    end
    savefig(p_rxns_with, joinpath(run_dir, "average_reactions_synthesized_with_fragments.pdf"))

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
        title = "Reaction Rebalancing (SynRXN rbl/$dataset)\nNumber of Reactions Synthesized for Successful Problems\n(Without Fragments)",
        legend = false,
        size = (600, 450),
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        yscale = :log10
    )
    max_rxn_without_val = 0.0
    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        if !isempty(vals)
            boxplot!(p_rxns_without, fill(i, length(vals)), vals, color=:lightgrey, fillalpha=0.7, bar_width=0.6, linecolor=:black)
            max_rxn_without_val = max(max_rxn_without_val, maximum(vals))
        end
    end
    max_rxn_without_val = max_rxn_without_val == 0.0 ? 1.0 : max_rxn_without_val
    boxplot!(p_rxns_without, ylims=(0.5, max_rxn_without_val * 2.5))

    for (i, m) in enumerate(metrics_order)
        sub_df = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        vals = sub_df.reactions_synthesized
        annot = isempty(vals) ? text("N/A", 8, :center, :bottom) : text("Avg: $(round(mean(vals); digits=1))", 8, :center, :bottom)
        annotate!(p_rxns_without, i, max_rxn_without_val * 1.5, annot)
    end
    savefig(p_rxns_without, joinpath(run_dir, "average_reactions_synthesized_without_fragments.pdf"))

    # Existing success and time plots
    success_with_metric = Float64[]
    success_without_metric = Float64[]
    success_counts_with = Int[]
    success_counts_without = Int[]
    time_with_metric = Float64[]
    time_without_metric = Float64[]
    time_with_empty = Bool[]
    time_without_empty = Bool[]

    for m in metrics_order
        row_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, summary_df)
        push!(success_with_metric, isempty(row_with) ? 0.0 : row_with[1, :success_rate])
        push!(success_counts_with, isempty(row_with) ? 0 : row_with[1, :successful_runs])
        push!(time_with_metric, isempty(row_with) ? 0.0 : row_with[1, :avg_synthesis_time_seconds])
        push!(time_with_empty, isempty(row_with) || row_with[1, :successful_runs] == 0)

        row_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, summary_df)
        push!(success_without_metric, isempty(row_without) ? 0.0 : row_without[1, :success_rate])
        push!(success_counts_without, isempty(row_without) ? 0 : row_without[1, :successful_runs])
        push!(time_without_metric, isempty(row_without) ? 0.0 : row_without[1, :avg_synthesis_time_seconds])
        push!(time_without_empty, isempty(row_without) || row_without[1, :successful_runs] == 0)
    end

    success_matrix = hcat(success_without_metric, success_with_metric)
    success_counts_matrix = hcat(success_counts_without, success_counts_with)
    p1 = groupedbar(metrics_order, success_matrix,
        xlabel = "Similarity Guidance Metric",
        ylabel = "Success Rate (%)",
        label = ["Without Fragments" "With Fragments"],
        title = "Reaction Rebalancing (SynRXN rbl/$dataset)\nSuccess Rate of Reaction Rebalancing ($synthesis_eval_limit Problems Analyzed)",
        legend = :right,
        ylimits = (0, 100),
        series_annotations = map(
            (v, c) -> Plots.text("$(round(v; digits=1))%\n($c)", 7, :center, :bottom), success_matrix, success_counts_matrix),
        color = [:lightgrey :dodgerblue],
        linewidth = 1.2,
        edgecolor = :black,
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8,
        size = (600, 450)
    )

    p2 = boxplot(
        xticks = (1:4, metrics_order),
        xlims = (0.5, 4.5),
        xlabel = "Similarity Guidance Metric",
        ylabel = "Synthesis Time (s)",
        title = "Reaction Rebalancing (SynRXN rbl/$dataset)\nRuntime of Successful Reaction Syntheses",
        legend = :topright,
        size = (600, 450),
        titlefontsize = 10,
        guidefontsize = 9,
        tickfontsize = 8,
        legendfontsize = 8
    )
    boxplot!(p2, [0], [0], label="Without Fragments", color=:lightgrey, seriestype=:shape, fillalpha=0.7, linecolor=:black)
    boxplot!(p2, [0], [0], label="With Fragments", color=:dodgerblue, seriestype=:shape, fillalpha=0.7, linecolor=:black)

    max_time_val = 0.0
    for (i, m) in enumerate(metrics_order)
        sub_df_without = filter(r -> r.similarity_metric == m && r.use_fragments == false, succ_rxn_df)
        vals_without = sub_df_without.elapsed_time_seconds
        if !isempty(vals_without)
            boxplot!(p2, fill(i - 0.22, length(vals_without)), vals_without, label="", color=:lightgrey, fillalpha=0.7, bar_width=0.35, linecolor=:black)
            max_time_val = max(max_time_val, maximum(vals_without))
        end

        sub_df_with = filter(r -> r.similarity_metric == m && r.use_fragments == true, succ_rxn_df)
        vals_with = sub_df_with.elapsed_time_seconds
        if !isempty(vals_with)
            boxplot!(p2, fill(i + 0.22, length(vals_with)), vals_with, label="", color=:dodgerblue, fillalpha=0.7, bar_width=0.35, linecolor=:black)
            max_time_val = max(max_time_val, maximum(vals_with))
        end
    end

    max_time_val = max_time_val == 0.0 ? 1.0 : max_time_val
    boxplot!(p2, ylims=(0, max_time_val * 1.3))

    # Save plots as PDF
    pdf_path1 = joinpath(run_dir, "synthesis_success_rate.pdf")
    savefig(p1, pdf_path1)

    pdf_path2 = joinpath(run_dir, "synthesis_average_time.pdf")
    savefig(p2, pdf_path2)

    # Save corresponding CSVs for success rate and average synthesis time
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
    println("  - Metadata:    $(joinpath(run_dir, "metadata.txt"))")
    println("  - Detailed CSV:$(joinpath(run_dir, "results.csv"))")
    println("  - Summary CSV: $(joinpath(run_dir, "summary.csv"))")
    println("  - Success CSV/Plot: $(joinpath(run_dir, "synthesis_success_rate.csv")) / .pdf")
    println("  - Time CSV/Plot:    $(joinpath(run_dir, "synthesis_average_time.csv")) / .pdf")
    println("  - 1st Feasibility CSV/Plot: $(joinpath(run_dir, "first_filter_feasibility.csv")) / .pdf")
    println("  - 1st Sizes CSV/Plot:       $(joinpath(run_dir, "first_filter_unfeasible_sizes.csv")) / .pdf")
    println("  - 2nd Success CSV/Plot:     $(joinpath(run_dir, "second_filter_success_rate.csv")) / .pdf")
    println("  - 2nd Sizes With CSV/Plot:   $(joinpath(run_dir, "second_filter_unfeasible_sizes_with_fragments.csv")) / .pdf")
    println("  - 2nd Sizes Without CSV/Plot:$(joinpath(run_dir, "second_filter_unfeasible_sizes_without_fragments.csv")) / .pdf")
    println("  - 2nd Runtime CSV/Plot:     $(joinpath(run_dir, "second_filter_average_runtime.csv")) / .pdf")
    println("  - Mols CSV/Plot:            $(joinpath(run_dir, "average_molecules_synthesized.csv")) / .pdf")
    println("  - Rxns With CSV/Plot:       $(joinpath(run_dir, "average_reactions_synthesized_with_fragments.csv")) / .pdf")
    println("  - Rxns Without CSV/Plot:    $(joinpath(run_dir, "average_reactions_synthesized_without_fragments.csv")) / .pdf")
    println("=======================================================")
end
