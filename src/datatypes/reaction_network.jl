struct Reaction
    rate::Union{Nothing, Float64}
    inputs::Vector{Tuple{Int, Molecule}}
    outputs::Vector{Tuple{Int, Molecule}}
    ignore_balanced::Bool
end
function Reaction(
        rate::Union{Nothing, Float64},
        inputs::Vector{Tuple{Int, Molecule}},
        outputs::Vector{Tuple{Int, Molecule}}
)
    Reaction(rate, inputs, outputs, false)
end
function Reaction(
        inputs::Vector{Tuple{Int, Molecule}}, outputs::Vector{Tuple{Int, Molecule}})
    Reaction(nothing, inputs, outputs, false)
end

import Base: hash, ==

function hash(reaction::Reaction, h::UInt)
    h = hash(reaction.rate, h)
    # Sort inputs and outputs to make hash order-insensitive
    sorted_inputs = sort(reaction.inputs; by = x -> x[2].canonical_smiles)
    sorted_outputs = sort(reaction.outputs; by = x -> x[2].canonical_smiles)
    for (count, molecule) in sorted_inputs
        h = hash((count, molecule), h)
    end
    for (count, molecule) in sorted_outputs
        h = hash((count, molecule), h)
    end
    return h
end

function ==(a::Reaction, b::Reaction)
    if a.rate != b.rate
        return false
    end
    if length(a.inputs) != length(b.inputs) || length(a.outputs) != length(b.outputs)
        return false
    end

    # Compare inputs (order-insensitive)
    inputs_a = sort(a.inputs; by = x -> x[2].canonical_smiles)
    inputs_b = sort(b.inputs; by = x -> x[2].canonical_smiles)
    if inputs_a != inputs_b
        return false
    end

    # Compare outputs (order-insensitive)
    outputs_a = sort(a.outputs; by = x -> x[2].canonical_smiles)
    outputs_b = sort(b.outputs; by = x -> x[2].canonical_smiles)
    if outputs_a != outputs_b
        return false
    end

    return true
end

function Reaction(inputs::Vector{Molecule}, outputs::Vector{Molecule})
    input_dict = Dict{Molecule, Int}()
    for (i, input) in enumerate(inputs)
        if haskey(input_dict, input)
            input_dict[input] += 1
        else
            input_dict[input] = 1
        end
    end
    inputs = [(input_dict[input], input) for input in keys(input_dict)]

    output_dict = Dict{Molecule, Int}()
    for (i, output) in enumerate(outputs)
        if haskey(output_dict, output)
            output_dict[output] += 1
        else
            output_dict[output] = 1
        end
    end
    outputs = [(output_dict[output], output) for output in keys(output_dict)]

    return Reaction(nothing, inputs, outputs, false)
end

struct ReactionNetwork
    reactions::Vector{Reaction}
end

function get_reactions(network::ReactionNetwork)
    return network.reactions
end

function get_molecules(network::ReactionNetwork)::Vector{Molecule}
    molecules = Set{Molecule}()
    for reaction in network.reactions
        for (count, molecule) in reaction.inputs
            push!(molecules, molecule)
        end
        for (count, molecule) in reaction.outputs
            push!(molecules, molecule)
        end
    end
    return collect(molecules)
end

function get_atoms(network::ReactionNetwork)::Vector{String}
    atoms = Set{String}()
    for reaction in network.reactions
        for (count, molecule) in reaction.inputs
            for atom in molecule.atoms
                push!(atoms, atom.name)
            end
        end
        for (count, molecule) in reaction.outputs
            for atom in molecule.atoms
                push!(atoms, atom.name)
            end
        end
    end
    return collect(atoms)
end

import Base.==
function ==(a::ReactionNetwork, b::ReactionNetwork)
    return Set(a.reactions) == Set(b.reactions)
end

import Base.hash
function hash(network::ReactionNetwork, h::UInt)
    # Order-insensitive hash: hash the set of reactions
    return hash(Set(network.reactions), h)
end

function add_reactions!(network::ReactionNetwork, new_reactions::Vector{Reaction})
    ReactionNetwork(vcat(network.reactions, new_reactions))
end

function is_valid(model)
    # Check if all the reactions are balanced
    for reaction in model.reactions
        # Skip if the reaction is not supposed to be balanced (e.g external model inputs and outputs)
        if reaction.ignore_balanced
            continue
        end

        input_species = Dict{String, Int}()
        output_species = Dict{String, Int}()

        for (count, input) in reaction.inputs
            for subspecies in input.atoms
                if haskey(input_species, subspecies.name)
                    input_species[subspecies.name] += count
                else
                    input_species[subspecies.name] = count
                end
            end
        end

        for (count, output) in reaction.outputs
            for subspecies in output.atoms
                if haskey(output_species, subspecies.name)
                    output_species[subspecies.name] += count
                else
                    output_species[subspecies.name] = count
                end
            end
        end

        if input_species != output_species
            return false
        end
    end

    # TODO: Check if a reaction is not the same on both sides

    return true
end

function compare_networks(
        network1::ReactionNetwork, network2::ReactionNetwork; check_rates::Bool = true
)
    # Check if the number of reactions is the same
    if length(network1.reactions) != length(network2.reactions)
        return false
    end

    molecules1 = get_molecules(network1)
    molecules2 = get_molecules(network2)

    if !issetequal(molecules1, molecules2)
        return false
    end

    # Check if the reactions are the same
    for reaction1 in network1.reactions
        found_reaction = false
        for reaction2 in network2.reactions
            if check_rates && reaction1.rate != reaction2.rate
                continue
            end

            if length(reaction1.inputs) != length(reaction2.inputs) ||
               length(reaction1.outputs) != length(reaction2.outputs)
                continue
            end

            found_reaction = true

            # println("inputs1: ", reaction1.inputs)
            # println("inputs2: ", reaction2.inputs)

            for input1 in reaction1.inputs
                found_input = false
                for input2 in reaction2.inputs
                    if input1[2] == input2[2] && input1[1] == input2[1]
                        println("found input: ", input1[2], " == ", input2[2])
                        found_input = true
                        break
                    end
                end
                if !found_input
                    found_reaction = false
                    break
                end
            end

            # for output1 in reaction1.outputs
            #     found_output = false
            #     for output2 in reaction2.outputs
            #         if output1[2] == output2[2] && output1[1] == output2[1]
            #             found_output = true
            #             break
            #         end
            #     end
            #     if !found_output
            #         found_reaction = false
            #         break
            #     end
            # end
        end
        if !found_reaction
            return false
        end
    end

    return true
end

# Function to print a single reaction
function Base.show(io::IO, reaction::Reaction)
    # Print complete reaction
    print(io, to_string(reaction))
end

function to_string(reaction::Reaction; compact::Bool = true)
    # Build reactants string
    reactants = if isempty(reaction.inputs)
        "∅"
    else
        join(
            ["$(input[1] > 1 ? "$(input[1])" : "")$(compact ? to_compact(input[2]) : to_SMILES(input[2]))"
             for input in reaction.inputs],
            " + "
        )
    end

    # Build arrow with rate
    arrow = isnothing(reaction.rate) ? " → " : " -($(reaction.rate))→ "

    # Build products string
    products = if isempty(reaction.outputs)
        "∅"
    else
        join(
            ["$(output[1] > 1 ? "$(output[1])" : "")$(compact ? to_compact(output[2]) : to_SMILES(output[2]))"
             for output in reaction.outputs],
            " + "
        )
    end

    return "$(reactants)$(arrow)$(products)"
end

# Function to neatly print a model
function Base.show(io::IO, model::ReactionNetwork)
    println(io, "Chemical Reaction Network:")
    for (i, reaction) in enumerate(model.reactions)
        print(io, "[$i] ", reaction)
        if i < length(model.reactions)
            print(io, "\n")
        end
    end
end

function to_string(model::ReactionNetwork; compact::Bool = true)
    result = "Chemical Reaction Network:\n"
    for (i, reaction) in enumerate(model.reactions)
        result *= "[$i] $(to_string(reaction, compact=compact))"
        if i < length(model.reactions)
            result *= "\n"
        end
    end
    return result
end

function count_species(network)
    # Count the number of unique species in the network
    unique_species = Set{Molecule}()
    for reaction in network.reactions
        for input in reaction.inputs
            push!(unique_species, input[2])
        end
        for output in reaction.outputs
            push!(unique_species, output[2])
        end
    end

    return length(unique_species)
end

function count_reactions(network)
    return length(network.reactions)
end

function get_bond_inventory(molecules)
    inventory = Dict{String, Int}()
    
    for (number, mol) in molecules
        for bond in mol.bonds
            atom1 = string(mol.atoms[bond.from].name)
            atom2 = string(mol.atoms[bond.to].name)
            
            pair = sort([atom1, atom2])
            bond_str = "$(pair[1])$(bond.bond_type)$(pair[2])"
            
            inventory[bond_str] = get(inventory, bond_str, 0) + number
        end
    end
    return inventory
end

function bonds_changed(reaction)
    input_inventory = get_bond_inventory(reaction.inputs)
    output_inventory = get_bond_inventory(reaction.outputs)
    
    all_bond_types = union(keys(input_inventory), keys(output_inventory))
    
    score = 0
    for bond_type in all_bond_types
        in_count = get(input_inventory, bond_type, 0)
        out_count = get(output_inventory, bond_type, 0)
        
        score += abs(in_count - out_count)
    end
    
    return score ÷ 2
end
