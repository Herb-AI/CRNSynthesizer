using CRNSynthesizer
using DataFrames
import MoleculeFlow
using DataStructures

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")
include("data/SYN_problems.jl")

include("data/synrxn_loader.jl")
import .SynRXNLoader

# -------------------------------------------------------------
# Caching wrappers to avoid expensive RDKit/MoleculeFlow calls
# -------------------------------------------------------------
const MOL_CACHE = Dict{String, MoleculeFlow.Molecule}()
const BRICS_CACHE = Dict{String, Vector{String}}()
const CLEANED_FRAG_CACHE = Dict{String, MoleculeFlow.Molecule}()
const PARSED_MOLECULE_CACHE = Dict{String, Molecule}()
const FRAG_ATOM_COUNT_CACHE = Dict{String, Int}()
const MOL_ATOM_COUNT_CACHE = Dict{String, Int}()

function get_mol_cached(smiles::String)::MoleculeFlow.Molecule
    get!(MOL_CACHE, smiles) do
        MoleculeFlow.mol_from_smiles(smiles)
    end
end

function get_brics_cached(smiles::String)::Vector{String}
    get!(BRICS_CACHE, smiles) do
        mol = get_mol_cached(smiles)
        if ismissing(mol)
            missing
        else
            MoleculeFlow.brics_decompose(mol; min_fragment_size = 2)
        end
    end
end

function get_cleaned_frag_cached(frag_smiles::String)::MoleculeFlow.Molecule
    get!(CLEANED_FRAG_CACHE, frag_smiles) do
        frag_m = get_mol_cached(frag_smiles)
        if ismissing(frag_m)
            missing
        else
            MoleculeFlow.delete_substructs(frag_m, "[#0]")
        end
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
    rxn_parts = split(rxn_str, ">>")
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

function parse_syn_problem(reaction_str::String)::ProblemDefinition
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
    reactant_counts = Dict{Molecule, Int}()
    for m in target_reactants
        reactant_counts[m] = get(reactant_counts, m, 0) + 1
    end
    inputs = [(reactant_counts[m], m) for m in unique(target_reactants)]

    product_counts = Dict{Molecule, Int}()
    for m in target_products
        product_counts[m] = get(product_counts, m, 0) + 1
    end
    outputs = [(product_counts[m], m) for m in unique(target_products)]

    reaction = CRNSynthesizer.Reaction(nothing, inputs, outputs, false)
    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Calculate known indices
    incomplete_molecules = Set{Molecule}(vcat(incomplete_reactants, incomplete_products))
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

    return ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network
    )
end

function is_feasible_problem(problem::ProblemDefinition)::Tuple{Bool, String}
    # Decompose all known molecules into BRICS fragments
    known_frags = Set{String}()
    for m in problem.known_molecules
        frags = get_brics_cached(m.canonical_smiles)
        if !ismissing(frags)
            for f in frags
                push!(known_frags, f)
            end
        end
    end

    # Check each goal molecule
    for goal_mol in problem.goal_molecules
        total_atoms = length(goal_mol.atoms)

        # A missing molecule is feasible if it has size <= 6 
        # small size due to increased branching factor after introduction of fragments
        if total_atoms <= 6
            continue
        end

        # Or if it shares a BRICS fragment (i.e. has a substructure match with any known BRICS fragment)
        goal_m = get_mol_cached(goal_mol.canonical_smiles)
        if ismissing(goal_m)
            return false,
            "Goal molecule $(goal_mol.canonical_smiles) could not be parsed by MoleculeFlow"
        end

        goal_atom_count = get_mol_atom_count_cached(goal_mol.canonical_smiles)

        shared_any = false
        for frag_smiles in known_frags
            # Pre-filter by atom count
            frag_atom_count = get_frag_atom_count_cached(frag_smiles)
            if frag_atom_count > goal_atom_count
                continue
            end

            cleaned_m = get_cleaned_frag_cached(frag_smiles)
            if !ismissing(cleaned_m)
                match_res = MoleculeFlow.has_substructure_match(goal_m, cleaned_m)
                if !ismissing(match_res) && match_res
                    shared_any = true
                    break
                end
            end
        end

        if !shared_any
            return false,
            "Goal molecule has > 6 atoms ($total_atoms) and does not share a BRICS fragment"
        end
    end

    return true, "Problem is feasible"
end

function run_problem_synthesis(
        problem::ProblemDefinition;
        max_time::Int = 60,
        metric::Symbol = :none,
        combine_method::Symbol = :multiplicative,
        max_stage::Symbol = :networks,
        use_fragments::Bool = true
)::Bool
    fragment_rules = Dict{Int, Set{Expr}}()
    starting_fragments = Expr[]

    if use_fragments
        for m in problem.known_molecules
            e, s = parse_molecule_to_fragment_rules(m.canonical_smiles)
            for (k, v) in e
                if haskey(fragment_rules, k)
                    union!(fragment_rules[k], v)
                else
                    fragment_rules[k] = v
                end
            end
            append!(starting_fragments, s)
        end
        starting_fragments = unique(starting_fragments)
    end

    if max_stage == :molecules
        # ---------------------------------------------------------
        # Atoms -> Molecules
        # ---------------------------------------------------------
        molecule_settings = SynthesizerSettings(;
            max_time = max_time, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
            options = Dict{Symbol, Any}(
                :unique_candidates => true, :similarity_metric => metric,
                :similarity_combine => combine_method)
        )

        molecules = Vector{Molecule}()
        try
            molecules = synthesize_molecules(
                problem.atom_valences; settings = molecule_settings, starting_fragments = starting_fragments,
                fragment_rules = fragment_rules
            )
        catch e
            println("Molecule synthesis failed with error: ", e)
            return false
        end

        println("[Atoms → Molecules] Found $(length(molecules)) molecules.")
        if issubset(problem.goal_molecules, molecules)
            println("      \033[32m✓ All goal molecules found.\033[0m")
        else
            println("      \033[31m✗ Missing goal molecules.\033[0m")
            return false
        end
        return true
    end

    # ---------------------------------------------------------
    # Molecules -> Reactions (using unified atom -> molecules -> reactions synthesizer)
    # ---------------------------------------------------------
    molecule_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )
    reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network),
        benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )

    candidates = Vector{CRNSynthesizer.Reaction}()
    found_molecules = OrderedSet{Molecule}()
    try
        candidates, found_molecules = synthesize_reactions(
            problem,
            molecule_settings,
            reaction_settings;
            initial_molecules_count = 1000,
            fragment_rules = fragment_rules,
            starting_fragments = starting_fragments
        )
    catch e
        println("[Problem → Reactions] Reaction synthesis failed with error: ", e)
        return false
    end

    println("[Problem → Reactions] Found $(length(candidates)) reactions, $(length(found_molecules)) molecules.")
    if issubset(problem.goal_molecules, found_molecules)
        println("      \033[32m✓ All goal molecules found.\033[0m")
    else
        println("      \033[31m✗ Missing goal molecules.\033[0m")
        return false
    end
    if issubset(get_reactions(problem.goal_network), candidates)
        println("      \033[32m✓ All goal reactions found.\033[0m")
    else
        println("      \033[31m✗ Missing goal reactions.\033[0m")
        return false
    end

    if max_stage == :reactions
        return true
    end

    # ---------------------------------------------------------
    # Full Pipeline (Atoms -> Molecules -> Reactions -> Networks)
    # ---------------------------------------------------------
    pipeline_molecule_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )
    pipeline_reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network),
        benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = [problem.goal_network], benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(
            :unique_candidates => true, :similarity_metric => metric,
            :similarity_combine => combine_method)
    )

    networks = Vector{ReactionNetwork}()
    reactions_found = Vector{CRNSynthesizer.Reaction}()
    found_molecules = Vector{Molecule}()
    try
        (networks,
            reactions_found,
            found_molecules) = synthesize_networks(
            problem,
            pipeline_molecule_settings,
            pipeline_reaction_settings,
            network_settings;
            initial_molecules_count = 1000,
            initial_reactions_count = 100000,
            fragment_rules = fragment_rules,
            starting_fragments = starting_fragments
        )
    catch e
        println("[Problem → Networks] Network synthesis failed with error: ", e)
        return false
    end

    println("[Problem → Networks] Found $(length(networks)) networks, $(length(reactions_found)) reactions, $(length(found_molecules)) molecules.")
    if issubset([problem.goal_network], networks)
        println("      \033[32m✓ Target network found.\033[0m")
    else
        println("      \033[31m✗ Target network not found.\033[0m")
    end
    if issubset(get_reactions(problem.goal_network), reactions_found)
        println("      \033[32m✓ All goal reactions found in pipeline.\033[0m")
    else
        println("      \033[31m✗ Missing goal reactions in pipeline.\033[0m")
    end
    if issubset(problem.goal_molecules, found_molecules)
        println("      \033[32m✓ All goal molecules found in pipeline.\033[0m")
    else
        println("      \033[31m✗ Missing goal molecules in pipeline.\033[0m")
    end

    return issubset([problem.goal_network], networks)
end

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
    ),
    (
        name = "Example missing molecule problem from SynRXN",
        problem = syn_problem()
    )
]

function run_hardcoded_benchmarks(max_time::Int = 60)
    println()
    println("=======================================================")
    println("Running Hardcoded Problems Feasibility Benchmark")
    println("=======================================================")

    for (name, problem) in PROBLEMS
        println("\n-------------------------------------------------------")
        println("\033[1mBenchmarking problem: $name\033[0m")

        is_feas, reason = is_feasible_problem(problem)
        if is_feas
            println("  ✓ Feasibility Check: \033[32mFEASIBLE\033[0m")
        else
            println("  ✗ Feasibility Check: \033[31mINFEASIBLE\033[0m ($reason)")
        end

        for metric in [:none, :simpson, :tanimoto, :both]
            println("\n    \033[1mSimilarity Metric: $metric\033[0m")

            fragment_flags = startswith(name, "Estherification") ? [true, false] : [true]

            for use_fragments in fragment_flags
                frag_str = use_fragments ? "With fragments" : "Without fragments"
                println("    \033[1m[$frag_str]\033[0m")
                println("    Running unified synthesis pipeline...")
                elapsed = @elapsed success = run_problem_synthesis(
                    problem; max_time = max_time, metric = metric,
                    max_stage = :networks, use_fragments = use_fragments
                )
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

# -------------------------------------------------------------
# Automated SynRXN mos Benchmark
# -------------------------------------------------------------
function run_automated_mos_benchmark(;
        max_time::Int = 120, max_scan::Int = 10, max_synthesis_runs::Int = 5)
    println()
    println("=======================================================")
    println("Running Automated SynRXN mos Feasibility Benchmark")
    println("=======================================================")

    df = SynRXNLoader.load_synrxn_dataset("rbl", "mos")
    total_records = nrow(df)
    scan_limit = min(total_records, max_scan)

    println("Successfully loaded SynRXN mos dataset with $total_records records.")
    println("Scanning the first $scan_limit records for feasibility...")

    feasible_problems = Tuple{String, ProblemDefinition}[]
    infeasible_reasons = Dict{String, Int}()

    feasible_count = 0
    infeasible_count = 0

    for idx in 1:scan_limit
        row = df[idx, :]
        r_id = row.r_id
        rxn_str = row.rxn
        gt_str = row.ground_truth

        reaction_str = rxn_str * "," * gt_str
        problem = nothing
        try
            problem = parse_syn_problem(reaction_str)
        catch e
            infeasible_count += 1
            println("  \033[31m✗ Failed to parse reaction for record $r_id: ", e, "\033[0m")
            reason = "Failed to parse reaction"
            infeasible_reasons[reason] = get(infeasible_reasons, reason, 0) + 1
            continue
        end
        try
            is_feas, reason = is_feasible_problem(problem)

            if is_feas
                push!(feasible_problems, (r_id, problem))
                feasible_count += 1
            else
                infeasible_count += 1
                infeasible_reasons[reason] = get(infeasible_reasons, reason, 0) + 1
            end
        catch e
            infeasible_count += 1
            println("  \033[31m✗ Failed to decide feasibility for record $r_id: ", e, "\033[0m")
            reason = "Failed to decide feasibility"
            infeasible_reasons[reason] = get(infeasible_reasons, reason, 0) + 1
        end
    end

    println("\nClassification Summary (first $scan_limit records):")
    println("  ✓ Feasible:   $feasible_count")
    println("  ✗ Infeasible: $infeasible_count")
    for (reason, count) in sort(collect(infeasible_reasons); by = x -> x[2], rev = true)
        println("    - $reason: $count")
    end

    # Run synthesis evaluation on a subset of feasible problems
    synthesis_eval_limit = min(length(feasible_problems), max_synthesis_runs)
    if synthesis_eval_limit > 0
        println("\nEvaluating synthesis (stages molecules -> reactions only) on the first $synthesis_eval_limit feasible problems...")

        successful_runs = 0
        total_time = 0.0

        for (i, (r_id, problem)) in enumerate(feasible_problems[1:synthesis_eval_limit])
            println("  Evaluating reaction $r_id ($i/$synthesis_eval_limit)...")
            elapsed = @elapsed success = run_problem_synthesis(
                problem; max_time = max_time, max_stage = :reactions
            )
            total_time += elapsed

            if success
                println("    \033[32m✓ Target reaction successfully synthesized in $(round(elapsed; digits=2))s!\033[0m")
                successful_runs += 1
            else
                println("    \033[31m✗ Synthesis failed or timed out in $(round(elapsed; digits=2))s.\033[0m")
            end
        end

        success_rate = (successful_runs / synthesis_eval_limit) * 100
        avg_time = total_time / synthesis_eval_limit

        println("\nSynthesis Performance on Feasible Subspace:")
        println("  - Evaluated runs: $synthesis_eval_limit")
        println("  - Success rate:   $(round(success_rate; digits=1))%")
        println("  - Average time:   $(round(avg_time; digits=2))s")
    else
        println("\nNo feasible problems found to evaluate.")
    end

    println("=======================================================")
end

run_hardcoded_benchmarks()

#run_automated_mos_benchmark(; max_scan = 10, max_synthesis_runs = 3)
