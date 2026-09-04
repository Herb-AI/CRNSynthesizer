@enum ReactionPosition begin
    INPUT
    OUTPUT
    UNKNOWN
end

struct RequiredMolecule
    molecule::Molecule
    position::ReactionPosition
end

struct ProblemDefinition
    atom_valences::OrderedDict{String, Int}
    known_molecules::Vector{Molecule}
    goal_molecules::Vector{Molecule}
    goal_network::ReactionNetwork
    expected_profiles::Dict{Molecule, Vector{Float64}}
    time_data::Vector{Float64}
    required_molecules::Vector{RequiredMolecule}
    partial_reaction::Union{Nothing, Reaction}

    function ProblemDefinition(
            atom_valences::OrderedDict{String, Int},
            known_molecules::Vector{Molecule},
            goal_molecules::Vector{Molecule},
            goal_network::ReactionNetwork;
            expected_profiles::Dict{Molecule, Vector{Float64}} = Dict{
                Molecule, Vector{Float64}}(),
            time_data::Vector{Float64} = Float64[],
            required_molecules::Vector{RequiredMolecule} = RequiredMolecule[],
            partial_reaction::Union{Nothing, Reaction} = nothing
    )
        if isempty(required_molecules)
            required_molecules = preprocess_problem(known_molecules, expected_profiles)
        end
        new(
            atom_valences,
            known_molecules,
            goal_molecules,
            goal_network,
            expected_profiles,
            time_data,
            required_molecules,
            partial_reaction
        )
    end
end

function preprocess_problem(
        known_molecules::Vector{Molecule}, expected_profiles::Dict{
            Molecule, Vector{Float64}}
)
    # Preprocess the problem to extract required molecules and their positions
    required_molecules = RequiredMolecule[]
    for mol in known_molecules
        if !(mol in keys(expected_profiles))
            push!(required_molecules, RequiredMolecule(mol, UNKNOWN))
        end
    end

    # TODO: look into the threshold and if that works with measurement noise
    for (molecule, profile) in expected_profiles
        # Calculate if there are significant increases and decreases
        increases = false
        decreases = false

        # Look for significant changes throughout the profile
        for i in eachindex(profile)[2:end]
            if profile[i] > profile[i - 1] * 1.2 # 20% increase threshold
                increases = true
            elseif profile[i] < profile[i - 1] * 0.8  # 20% decrease threshold
                decreases = true
            end
        end

        # Classify based on behavior
        if increases && decreases
            # Both increases and decreases - add as both input and output
            push!(required_molecules, RequiredMolecule(molecule, INPUT))
            push!(required_molecules, RequiredMolecule(molecule, OUTPUT))
        elseif increases
            push!(required_molecules, RequiredMolecule(molecule, OUTPUT))
        elseif decreases
            push!(required_molecules, RequiredMolecule(molecule, INPUT))
        else
            push!(required_molecules, RequiredMolecule(molecule, UNKNOWN))
        end
    end

    return required_molecules
end

#=function get_atoms(problem::ProblemDefinition)::Vector{String}
    atoms = Set{String}()
    for molecule in problem.known_molecules
        for atom in molecule.atoms
            push!(atoms, atom.name)
        end
    end
    for reaction in problem.known_reactions
        for input in reaction.inputs
            for atom in input[2].atoms
                push!(atoms, atom.name)
            end
        end
        for output in reaction.outputs
            for atom in output[2].atoms
                push!(atoms, atom.name)
            end
        end
    end

    return collect(atoms)
end=#

#=function get_molecules(problem::ProblemDefinition)
    molecules = Set{Molecule}()
    for molecule in problem.known_molecules
        push!(molecules, molecule)
    end
    for reaction in problem.known_reactions
        for input in reaction.inputs
            push!(molecules, input[2])
        end
        for output in reaction.outputs
            push!(molecules, output[2])
        end
    end

    return collect(molecules)
end=#

# function herb_cost_function(problem::ProblemDefinition)
#     sol -> begin
#         try
#             total_cost = 0.0
#             for (species_name, expected_profile) in problem.expected_profiles
#                 sol_species = sol(problem.time_data, idxs=[Symbol(species_name)])
#                 sol_species = collect(Iterators.flatten(sol_species))
#                 if isnothing(sol_species)
#                     total_cost += 1000.0
#                 else
#                     cost = sum((sol_species .- expected_profile).^2)
#                     total_cost += cost    end
#             end
#             return total_cost
#         catch e
#             println("Error in cost function: ", e)
#             return Inf
#         end
#     end
# end
