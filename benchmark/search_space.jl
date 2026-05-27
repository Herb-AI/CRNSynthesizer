using CRNSynthesizer, DataStructures
include("data/methane.jl")

# Problem Definition
atom_valences = OrderedDict("[H]" => 1, "[O]" => 2, "[C]" => 4)
max_time = 60

# ----------------------------------------------------------
# ----------------------- Molecules -----------------------
# ----------------------------------------------------------

# Synthesizer Settings
max_depth = 6

# Get a baseline for the amount of unique molecules that can be generated
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    molecules = synthesize_molecules(atom_valences; settings)
end
println(
    "Baseline: Generated $(length(molecules)) unique molecules in $(elapsed_time) seconds."
)

#= Gather the result of a synthesizer run with the ValidSMILES constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(max_time = max_time, max_depth = max_depth)
    molecules = synthesize_molecules(atoms, settings)
end
println(
    "With ValidSMILES constraint: Generated $(length(molecules)) molecules in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the ValidSMILES constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_valid_smiles => true)
    )
    molecules = synthesize_molecules(atoms, settings)
end
println(
    "Without ValidSMILES constraint: Generated $(length(molecules)) molecules in $(elapsed_time) seconds.",
) =#

# ----------------------------------------------------------
# ------------------- Reactions from Atoms -----------------
# ----------------------------------------------------------
# Reactions from Atoms is not supported with fragment rules
# Based on Wijers conclusions that pipelines without molecule step are not feasible
#= Problem Definition
max_depth = 8

# Get a baseline for the amount of unique reactions that can be generated
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    candidates = synthesize_reactions(atoms, settings)
end
println(
    "Baseline: Generated $(length(candidates)) unique reactions in $(elapsed_time) seconds."
)

# Gather the result of a synthesizer run with both the BalancedReaction constraint and the ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(max_time = max_time, max_depth = max_depth)
    candidates = synthesize_reactions(atoms, settings)
end
println(
    "With BalancedReaction and Ordered constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the BalancedReaction constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_balanced_reaction => true)
    )
    candidates = synthesize_reactions(atoms, settings)
end
println(
    "Without BalancedReaction constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the Ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_ordered_molecule_list => true)
    )
    candidates = synthesize_reactions(atoms, settings)
end
println(
    "Without Ordered constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without both the BalancedReaction and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(
            :disable_balanced_reaction => true, :disable_ordered_molecule_list => true
        )
    )
    candidates = synthesize_reactions(atoms, settings)
end
println(
    "Without BalancedReaction and Ordered constraints: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)
=#
# ----------------------------------------------------------
# ------------------- Reactions from Molecules -------------
# ----------------------------------------------------------

# Get a baseline for the amount of unique reactions that can be generated
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    candidates = synthesize_reactions(molecules, settings)
end
println(
    "Baseline: Generated $(length(candidates)) unique reactions in $(elapsed_time) seconds."
)

#= Gather the result of a synthesizer run with both the BalancedReaction constraint and the ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(max_time = max_time, max_depth = max_depth)
    candidates = synthesize_reactions(molecules, settings)
end
println(
    "With BalancedReaction and Ordered constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the BalancedReaction constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_balanced_reaction => true)
    )
    candidates = synthesize_reactions(molecules, settings)
end
println(
    "Without BalancedReaction constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the Ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_ordered_molecule_list => true)
    )
    candidates = synthesize_reactions(molecules, settings)
end
println(
    "Without Ordered constraint: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without both the BalancedReaction and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(
            :disable_balanced_reaction => true, :disable_ordered_molecule_list => true
        )
    )
    candidates = synthesize_reactions(molecules, settings)
end
println(
    "Without BalancedReaction and Ordered constraints: Generated $(length(candidates)) reactions in $(elapsed_time) seconds.",
) =#

# ----------------------------------------------------------
# ------------------- Networks from Atoms ------------------
# ----------------------------------------------------------

# # Problem Definition
# max_depth = 10

# # Get a baseline for the amount of unique networks that can be generated
# elapsed_time = @elapsed begin
#     settings = SynthesizerSettings(
#         max_time = max_time,
#         max_depth = max_depth,
#         options = Dict{Symbol, Any}(
#             :unique_candidates => true
#         ),
#     )
#     networks = synthesize_networks(atoms, settings; problem=problem)
# end
# println("Baseline: Generated $(length(networks)) unique networks in $(elapsed_time) seconds.")

# # Gather the result of a synthesizer run with both the ContainsMolecules and Ordered constraints
# elapsed_time = @elapsed begin
#     settings = SynthesizerSettings(
#         max_time = max_time,
#         max_depth = max_depth
#     )
#     networks = synthesize_networks(atoms, settings)
# end
# println("With ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.")

# # Gather the result of a synthesizer run without the ContainsMolecules constraint
# elapsed_time = @elapsed begin
#     settings = SynthesizerSettings(
#         max_time = max_time,
#         max_depth = max_depth,
#         options = Dict{Symbol, Any}(
#             :disable_contains_molecules => true,
#         ),
#     )
#     networks = synthesize_networks(atoms, settings)
# end
# println("Without ContainsMolecules constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.")

# # Gather the result of a synthesizer run without the Ordered constraint
# elapsed_time = @elapsed begin
#     settings = SynthesizerSettings(
#         max_time = max_time,
#         max_depth = max_depth,
#         options = Dict{Symbol, Any}(
#             :disable_ordered_reaction_list => true,
#         ),
#     )
#     networks = synthesize_networks(atoms, settings)
# end
# println("Without Ordered constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.")

# # Gather the result of a synthesizer run without both the ContainsMolecules and Ordered constraints
# elapsed_time = @elapsed begin
#     settings = SynthesizerSettings(
#         max_time = max_time,
#         max_depth = max_depth,
#         options = Dict{Symbol, Any}(
#             :disable_contains_molecules => true,
#             :disable_ordered_reaction_list => true,
#         ),
#     )
#     networks = synthesize_networks(atoms, settings)
# end
# println("Without ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.")

# ----------------------------------------------------------
# ------------------- Networks from Molecules -------------
# ----------------------------------------------------------

#= Problem Definition
max_depth = 7
molecules = synthesize_molecules(
    atoms, SynthesizerSettings(; max_programs = 100, max_depth = 9)
)
molecules = unique(molecules)
molecules = molecules[1:10]

# Get a baseline for the amount of unique networks that can be generated
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    networks = synthesize_networks(molecules, settings; problem = problem)
end
println(
    "Baseline: Generated $(length(networks)) unique networks in $(elapsed_time) seconds."
)

# Gather the result of a synthesizer run with both the ContainsMolecules and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(max_time = max_time, max_depth = max_depth)
    networks = synthesize_networks(molecules, settings; problem = problem)
end
println(
    "With ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the ContainsMolecules constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_contains_molecules => true)
    )
    networks = synthesize_networks(molecules, settings; problem = problem)
end
println(
    "Without ContainsMolecules constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the Ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_ordered_reaction_list => true)
    )
    networks = synthesize_networks(molecules, settings; problem = problem)
end
println(
    "Without Ordered constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without both the ContainsMolecules and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(
            :disable_contains_molecules => true, :disable_ordered_reaction_list => true
        )
    )
    networks = synthesize_networks(molecules, settings; problem = problem)
end
println(
    "Without ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
) =#

# ----------------------------------------------------------
# ------------------- Networks from Reactions --------------
# ----------------------------------------------------------

# Problem Definition
max_depth = 4
molecules = synthesize_molecules(
    atom_valences; settings = SynthesizerSettings(; max_programs = 100, max_depth = 9)
)
molecules = unique(molecules)
molecules = molecules[1:27]
reactions = synthesize_reactions(
    molecules, SynthesizerSettings(; max_programs = 200, max_depth = 8)
)
reactions = unique(reactions)

# Get a baseline for the amount of unique networks that can be generated
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    networks = synthesize_networks(reactions, settings, methane_problem(;
            selected_known_indices = [1, 3], selected_expected_indices = [1, 3]
        ))
end
println(
    "Baseline: Generated $(length(networks)) unique networks in $(elapsed_time) seconds."
)

#= Gather the result of a synthesizer run with both the ContainsMolecules and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(max_time = max_time, max_depth = max_depth)
    networks = synthesize_networks(reactions, settings; problem = problem)
end
println(
    "With ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the ContainsMolecules constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_contains_molecules => true)
    )
    networks = synthesize_networks(reactions, settings; problem = problem)
end
println(
    "Without ContainsMolecules constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without the Ordered constraint
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:disable_ordered_reaction_list => true)
    )
    networks = synthesize_networks(reactions, settings; problem = problem)
end
println(
    "Without Ordered constraint: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
)

# Gather the result of a synthesizer run without both the ContainsMolecules and Ordered constraints
elapsed_time = @elapsed begin
    settings = SynthesizerSettings(
        max_time = max_time,
        max_depth = max_depth,
        options = Dict{Symbol, Any}(
            :disable_contains_molecules => true, :disable_ordered_reaction_list => true
        )
    )
    networks = synthesize_networks(reactions, settings; problem = problem)
end
println(
    "Without ContainsMolecules and Ordered constraints: Generated $(length(networks)) networks in $(elapsed_time) seconds.",
) =#

# ----------------------------------------------------------
# ----------------- BRICS Fragments Evaluation -------------
# ----------------------------------------------------------
println("\n\n=======================================================")
println("Evaluating BRICS Fragments on Search Space Size")
println("=======================================================")

include("data/estherification.jl")
include("data/water.jl")
include("data/ethylene.jl")

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
    )
]

for (name, problem) in PROBLEMS
    println("\n-------------------------------------------------------")
    println("\033[1mBenchmarking Search Space: $name\033[0m")

    # Extract and merge fragments from the target molecules
    fragment_rules = Dict{Int, Set{Expr}}()
    starting_fragments = Expr[]
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
    println("Extracted $(length(starting_fragments)) unique starting fragments from target molecules.")

    for use_fragments in [false, true]
        println("\n  \033[1mUse Fragments: $use_fragments\033[0m")
        
        fr = use_fragments ? fragment_rules : Dict{Int, Set{Expr}}()
        sf = use_fragments ? starting_fragments : Expr[]

        # Molecules Search Space
        molecule_settings = SynthesizerSettings(
            max_depth = 6, options = Dict{Symbol, Any}(:unique_candidates => true)
        )
        elapsed_time_m = @elapsed mols = synthesize_molecules(
            problem.atom_valences; settings = molecule_settings, fragment_rules=fr, starting_fragments=sf
        )
        println("  [Atoms → Molecules] Exhaustively found $(length(mols)) valid molecules in $(elapsed_time_m) seconds.")

        # Reactions Search Space
        reaction_settings = SynthesizerSettings(
            max_depth = 6, options = Dict{Symbol, Any}(:unique_candidates => true)
        )
        elapsed_time_r = @elapsed reacts = synthesize_reactions(
            mols, reaction_settings; known_molecules = problem.known_molecules
        )
        println("  [Molecules → Reactions] Exhaustively found $(length(reacts)) valid reactions in $(elapsed_time_r) seconds.")
    end
end
