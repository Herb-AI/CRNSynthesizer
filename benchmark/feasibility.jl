using CRNSynthesizer
using DataFrames
import MoleculeFlow

include("data/estherification.jl")
include("data/water.jl")
include("data/methane.jl")
include("data/ethylene.jl")
include("data/SYN_problems.jl")

include("data/synrxn_loader.jl")
import .SynRXNLoader

function is_feasible_problem(problem::ProblemDefinition)::Tuple{Bool, String}
    # Decompose all known molecules into BRICS fragments
    known_frags = Set{String}()
    for m in problem.known_molecules
        frags = MoleculeFlow.brics_decompose(
            MoleculeFlow.mol_from_smiles(m.canonical_smiles); min_fragment_size = 2
        )
        if !ismissing(frags)
            for f in frags
                push!(known_frags, f)
            end
        end
    end
    
    # Check each goal molecule
    for goal_mol in problem.goal_molecules
        total_atoms = length(goal_mol.atoms)
        
        # A missing molecule is feasible if it has size <= 10
        if total_atoms <= 10
            continue
        end
        
        # Or if it shares a BRICS fragment (i.e. has a substructure match with any known BRICS fragment)
        goal_m = MoleculeFlow.mol_from_smiles(goal_mol.canonical_smiles)
        if ismissing(goal_m)
            return false, "Goal molecule $(goal_mol.canonical_smiles) could not be parsed by MoleculeFlow"
        end
        
        shared_any = false
        for frag_smiles in known_frags
            frag_m = MoleculeFlow.mol_from_smiles(frag_smiles)
            if !ismissing(frag_m)
                if MoleculeFlow.has_substructure_match(goal_m, frag_m)
                    shared_any = true
                    break
                end
            end
        end
        
        if !shared_any
            return false, "Goal molecule has > 10 atoms ($total_atoms) and does not share a BRICS fragment"
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

    # ---------------------------------------------------------
    # Atoms -> Molecules
    # ---------------------------------------------------------
    molecule_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 9, goal = problem.goal_molecules, benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(:unique_candidates => true, :similarity_metric => metric, :similarity_combine => combine_method)
    )
    
    molecules = Vector{Molecule}()
    try
        molecules = synthesize_molecules(
            problem.atom_valences; settings = molecule_settings, starting_fragments = starting_fragments, fragment_rules = fragment_rules
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
    
    if max_stage == :molecules
        return true
    end

    # ---------------------------------------------------------
    # Molecules -> Reactions
    # ---------------------------------------------------------
    reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network), benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(:unique_candidates => true, :similarity_metric => metric, :similarity_combine => combine_method)
    )
    molecules_for_reactions = unique(vcat(problem.known_molecules, molecules))
    
    candidates = Vector{CRNSynthesizer.Reaction}()
    try
        candidates = synthesize_reactions(
            molecules_for_reactions, reaction_settings; known_molecules = problem.known_molecules
        )
    catch e
        println("[Molecules → Reactions] Reaction synthesis failed with error: ", e)
        return false
    end
    
    println("[Molecules → Reactions] Found $(length(candidates)) reactions.")
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
        options = Dict{Symbol, Any}(:unique_candidates => true, :similarity_metric => metric, :similarity_combine => combine_method)
    )
    pipeline_reaction_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = get_reactions(problem.goal_network), benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(:unique_candidates => true, :similarity_metric => metric, :similarity_combine => combine_method)
    )
    network_settings = SynthesizerSettings(;
        max_time = max_time, max_depth = 5, goal = [problem.goal_network], benchmark_type = UntilFound,
        options = Dict{Symbol, Any}(:unique_candidates => true, :similarity_metric => metric, :similarity_combine => combine_method)
    )
    
    networks = Vector{ReactionNetwork}()
    reactions_found = Vector{CRNSynthesizer.Reaction}()
    found_molecules = Vector{Molecule}()
    try
        (networks, reactions_found, found_molecules) = synthesize_networks(
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
        problem = syn_problem(),
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
                    problem; max_time = max_time, metric = metric, max_stage = :networks, use_fragments = use_fragments
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
function run_automated_mos_benchmark(; max_time::Int = 60, max_scan::Int = typemax(Int), max_synthesis_runs::Int = 5)
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

#run_hardcoded_benchmarks()

run_automated_mos_benchmark(; max_scan = 8, max_synthesis_runs = 3)
