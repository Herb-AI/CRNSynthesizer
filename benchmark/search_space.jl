using CRNSynthesizer

# Problem Definition
include("data/methane.jl")
problem = methane_problem(;
    selected_known_indices = [1, 3], selected_expected_indices = [1, 3])
atoms = get_atoms(problem)
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
    molecules = synthesize_molecules(atoms, settings)
end
println(
    "Baseline: Generated $(length(molecules)) unique molecules in $(elapsed_time) seconds."
)

# Gather the result of a synthesizer run with the ValidSMILES constraint
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
)

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

# Problem Definition
max_depth = 5
molecules = synthesize_molecules(
    atoms, SynthesizerSettings(; max_programs = 100, max_depth = 9)
)
molecules = unique(molecules)
molecules = molecules[1:10]

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

# Gather the result of a synthesizer run with both the BalancedReaction constraint and the ordered constraint
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
)

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

# Problem Definition
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
)

# ----------------------------------------------------------
# ------------------- Networks from Reactions --------------
# ----------------------------------------------------------

# Problem Definition
max_depth = 4
molecules = synthesize_molecules(
    atoms, SynthesizerSettings(; max_programs = 100, max_depth = 9)
)
molecules = unique(molecules)
molecules = molecules[1:30]
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
    networks = synthesize_networks(reactions, settings; problem = problem)
end
println(
    "Baseline: Generated $(length(networks)) unique networks in $(elapsed_time) seconds."
)

# Gather the result of a synthesizer run with both the ContainsMolecules and Ordered constraints
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
)
