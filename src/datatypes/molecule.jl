@enum BondType single double aromatic triple quadruple

function to_string(bond_type::BondType)::String
    return if bond_type == single
        "-"
    elseif bond_type == double
        "="
    elseif bond_type == aromatic
        ":"
    elseif bond_type == triple
        "≡"
    else
        "≣"
    end
end

struct Atom
    name::String
    is_aromatic::Bool
    Atom(name::String, is_aromatic::Bool = false) = new(name, is_aromatic)
end

abstract type AbstractBond end

struct Bond <: AbstractBond
    from::Int
    to::Int
    bond_type::BondType

    function Bond(from::Int, to::Int, bond_type::BondType)
        from, to = minmax(from, to)
        return new(from, to, bond_type)
    end
end

struct Molecule
    atoms::Vector{Atom}
    bonds::Vector{AbstractBond}
    canonical_smiles::String
    fingerprint::Vector{UInt8}
    morgan_fingerprint::Vector{UInt8}

    function Molecule(atoms::Vector{Atom}, bonds::Vector{Bond}, canonical_smiles::String,
            fingerprint::Vector{UInt8}, morgan_fingerprint::Vector{UInt8} = UInt8[])
        # Create a copy of atoms and sort them
        sorted_indices = sortperm(atoms; by = atom -> atom.name)
        sorted_atoms = atoms[sorted_indices]

        # Create a mapping from old indices to new indices
        index_map = Dict(
            old_idx => new_idx for (new_idx, old_idx) in enumerate(sorted_indices)
        )

        # Adjust bonds to reflect the new atom indices
        adjusted_bonds = map(bonds) do bond
            if bond isa Bond
                new_from = index_map[bond.from]
                new_to = index_map[bond.to]
                Bond(new_from, new_to, bond.bond_type)
            else
                bond  # For other types of bonds
            end
        end

        # Sort the bonds by the new atom indices
        sort!(adjusted_bonds; by = bond -> (bond.from, bond.to))

        return new(
            sorted_atoms, adjusted_bonds, canonical_smiles, fingerprint, morgan_fingerprint)
    end
end

# TODO: Think of a better Molecule struct that doesn't have this many symmetries
import Base.==
function ==(a::Molecule, b::Molecule)
    if length(a.atoms) != length(b.atoms) || length(a.bonds) != length(b.bonds)
        return false
    end

    if a.canonical_smiles != b.canonical_smiles
        return false
    end

    return true
end

import Base.hash
function hash(m::Molecule, h::UInt)
    h = hash(:Molecule, h)
    # println("Hashing molecule: ", m)
    h = hash(m.canonical_smiles, h)
    return h
end

# import Base.hash
# function hash(m::Molecule, h::UInt)
#     h = hash(:Molecule, h)
#     # Hash sorted atom names
#     atom_names = Tuple(atom.name for atom in m.atoms)
#     h = hash(atom_names, h)
#     # Hash sorted bond tuples (from, to, bond_type as Int)
#     bond_tuples = Tuple((b.from, b.to, Int(b.bond_type)) for b in m.bonds if b isa Bond)
#     h = hash(bond_tuples, h)
#     return h
# end

Base.show(io::IO, atom::Atom) = print(io, atom.name)

function Base.show(io::IO, bond::Bond)
    print(io, "Bond(", bond.from, ", ", bond.to, ", ", to_string(bond.bond_type), ")")
end

Base.show(io::IO, molecule::Molecule) = begin
    print(io, "Molecule(atoms=[")
    print(io, join(molecule.atoms, ", "))
    print(io, "], bonds=[")
    print(io, join(molecule.bonds, ", "))
    print(io, "])")
end

function count_atoms(molecule::Molecule)::Dict{String, Int}
    atoms = Dict{String, Int}()
    net_charge = 0
    for atom in molecule.atoms
        # Expected format: [ElementCharge]
        m = match(r"^\[([A-Z][a-z]?)(.*)\]$", atom.name)
        if m !== nothing
            element = m.captures[1]
            charge_str = m.captures[2]

            atoms[element] = get(atoms, element, 0) + 1

            if !isempty(charge_str)
                sign = occursin("-", charge_str) ? -1 : 1
                num_str = filter(isdigit, charge_str)
                val = isempty(num_str) ?
                      length(filter(c -> c == '+' || c == '-', charge_str)) :
                      parse(Int, num_str)
                net_charge += sign * val
            end
        else
            atoms[atom.name] = get(atoms, atom.name, 0) + 1
        end
    end

    if net_charge != 0
        atoms["charge"] = net_charge
    end

    return atoms
end

function get_valences_from_molecules(molecules::Vector{Molecule})::OrderedDict{String, Int}
    valences = OrderedDict{String, Int}()

    for mol in molecules
        # Check if we need to parse this molecule
        needs_parsing = false
        for atom in mol.atoms
            if !haskey(valences, atom.name)
                needs_parsing = true
                break
            end
        end

        if !needs_parsing
            continue
        end

        if isempty(mol.canonical_smiles)
            throw(ArgumentError(string("Invalid input molecule: ", to_SMILES(mol))))
        end

        mf_mol = MoleculeFlow.mol_from_smiles(mol.canonical_smiles)
        if isnothing(mf_mol) || !mf_mol.valid
            throw(ArgumentError(string("Invalid input molecule: ", to_SMILES(mol))))
        end
        mf_mol = MoleculeFlow.add_hs(mf_mol)

        for a in MoleculeFlow.get_atoms(mf_mol)
            symbol = MoleculeFlow.get_symbol(a)
            charge = MoleculeFlow.get_formal_charge(a)

            name = if charge == 0
                "[$symbol]"
            elseif charge == 1
                "[$symbol+]"
            elseif charge == -1
                "[$symbol-]"
            elseif charge > 1
                "[$symbol+$charge]"
            else
                "[$symbol$charge]"
            end
            valences[name] = MoleculeFlow.get_valence(a)
        end
    end
    return valences
end

function to_compact(molecule::Molecule)
    # Get the atom counts
    atoms = count_atoms(molecule)

    # Convert the compact representation to a string
    compact_str = ""
    for (atom_name, count) in sort(collect(atoms); by = x -> x[1])
        if count > 1
            compact_str *= "$atom_name$count"
        else
            compact_str *= atom_name
        end
    end

    return convert_to_subscript(compact_str)
end

function from_SMILES(smiles::String)
    # Extract all atoms
    atoms = Atom[]
    atom_matches = collect(eachmatch(r"\[(.*?)\]", smiles))
    for regmatch::RegexMatch in atom_matches
        atom_name = regmatch.captures[1]  # Get the content inside brackets
        is_aromatic = islowercase(atom_name[1])
        push!(atoms, Atom("[" * uppercasefirst(atom_name) * "]", is_aromatic))
    end

    # Replace atoms with placeholders for easier parsing
    processed_smiles = smiles
    for (i, regmatch) in enumerate(atom_matches)
        submatch::String = Base.String(regmatch.match)
        processed_smiles = replace(processed_smiles, submatch => "[A$i]"; count = 1)
    end

    current_atom_idx::Int = 0
    branch_stack = Int[]  # Stack to keep track of branching points
    ring_connections = Dict{Int, Tuple{Int, BondType}}()  # Store atom idx and bond type
    current_bond_type::BondType = single  # Default to single bond
    bonds = Bond[]  # Store bonds

    # Iterate safely through characters, handling Unicode properly
    i::Int = firstindex(processed_smiles)
    while i <= lastindex(processed_smiles)
        char = processed_smiles[i]

        if char == '['
            # Start of an atom placeholder
            atom_end = findnext(']', processed_smiles, i)
            if isnothing(atom_end)
                error("Unmatched '[' in SMILES string: $smiles")
            end

            atom_placeholder = processed_smiles[i:atom_end]
            atom_idx = parse(
                Int, match(r"A(\d+)", atom_placeholder).captures[1]::SubString{String}
            )

            if current_atom_idx != 0  # If there's a previous atom, create a bond
                # Create bonds in both directions
                bond = Bond(current_atom_idx, atom_idx, current_bond_type)
                push!(bonds, bond)
            end

            current_atom_idx = atom_idx
            i = atom_end
        elseif char == '-'
            # Single bond (already default)
            current_bond_type = single
        elseif char == '='
            # Double bond
            current_bond_type = double
        elseif char == ':'
            current_bond_type = aromatic
        elseif char == '#'
            # Triple bond
            current_bond_type = triple
        elseif char == '$'
            # Quadruple bond
            current_bond_type = quadruple
        elseif char == '≡'
            # Triple bond (Unicode)
            current_bond_type = triple
        elseif char == '≣'
            # Quadruple bond (Unicode)
            current_bond_type = quadruple
        elseif char == '('
            # Start of a branch
            push!(branch_stack, current_atom_idx)
        elseif char == ')'
            # End of a branch, return to parent atom
            if !isempty(branch_stack)
                current_atom_idx = pop!(branch_stack)
            end
        elseif char == '%' || isdigit(char)
            # Ring closure
            if char == '%'
                # Find all subsequent digits
                j = nextind(processed_smiles, i)
                while j <= lastindex(processed_smiles) && isdigit(processed_smiles[j])
                    j = nextind(processed_smiles, j)
                end

                # Extract the multi-digit number
                ring_str = processed_smiles[nextind(
                    processed_smiles, i):prevind(processed_smiles, j)]
                ring_num = parse(Int, ring_str)

                i = prevind(processed_smiles, j)
            else
                ring_num = parse(Int, string(char))
            end

            if haskey(ring_connections, ring_num)
                # Close the ring
                other_atom_idx, ring_bond_type = ring_connections[ring_num]

                # Create bond
                bond = Bond(current_atom_idx, other_atom_idx, current_bond_type)
                push!(bonds, bond)

                # Reset bond type and remove the ring connection
                current_bond_type = single
                delete!(ring_connections, ring_num)
            else
                # Start the ring
                ring_connections[ring_num] = (current_atom_idx, current_bond_type)
                current_bond_type = single  # Reset for next bond
            end
        end
        i = nextind(processed_smiles, i)
    end

    mol = RDKitMinimalLib.get_mol(replace(smiles, "≡" => "#"))
    if isnothing(mol)
        throw(ArgumentError("Invalid molecule: $smiles"))
    end
    canonical_smiles = isnothing(mol) ? smiles : RDKitMinimalLib.get_smiles(mol)
    fingerprint = isnothing(mol) ? UInt8[] : RDKitMinimalLib.get_rdkit_fp_as_bytes(mol)
    morgan_fingerprint = isnothing(mol) ? UInt8[] :
                         RDKitMinimalLib.get_morgan_fp_as_bytes(
        mol, Dict{String, Any}("radius" => 2, "nBits" => 1024))
    return Molecule(atoms, bonds, canonical_smiles, fingerprint, morgan_fingerprint)
end

function to_SMILES(molecule::Molecule)::String
    if isempty(molecule.atoms)
        return ""
    end

    if length(molecule.atoms) == 1
        atom = molecule.atoms[1]
        atom_name_str = atom.is_aromatic ? lowercase(atom.name) : atom.name
        return atom_name_str
    end

    # Create an adjacency and ringbond dict
    adjacency = Dict{Int, Vector{Tuple{Int, String}}}()
    ringbonds = Dict{Int, Vector{Tuple{Int, String}}}()
    ring_digit = 1

    visited_atoms = falses(length(molecule.atoms))
    visited_bonds = Dict{AbstractBond, Bool}()
    for bond in molecule.bonds
        visited_bonds[bond] = false
    end
    bonds_queue = Vector{Bond}()
    for bond in molecule.bonds
        if bond.from == 1
            push!(bonds_queue, bond)
        end
    end

    while !isempty(bonds_queue)
        bond = popfirst!(bonds_queue)
        visited_bonds[bond] = true

        # println("Processing bond: ", bond)

        if visited_atoms[bond.to] && visited_atoms[bond.from]
            if !haskey(ringbonds, bond.to)
                ringbonds[bond.to] = Vector{Int}()
            end
            if !haskey(ringbonds, bond.from)
                ringbonds[bond.from] = Vector{Int}()
            end
            push!(ringbonds[bond.to], (ring_digit, to_string(bond.bond_type)))
            push!(ringbonds[bond.from], (ring_digit, to_string(bond.bond_type)))
            ring_digit += 1
        else
            if !haskey(adjacency, bond.from)
                adjacency[bond.from] = Vector{Int}()
            end
            if !haskey(adjacency, bond.to)
                adjacency[bond.to] = Vector{Int}()
            end

            if visited_atoms[bond.to]
                push!(adjacency[bond.to], (bond.from, to_string(bond.bond_type)))
            else
                push!(adjacency[bond.from], (bond.to, to_string(bond.bond_type)))
            end

            visited_atoms[bond.from] = true
            visited_atoms[bond.to] = true

            for b in molecule.bonds
                if !(b in bonds_queue) &&
                   !visited_bonds[b] &&
                   (
                       b.from == bond.to ||
                       b.to == bond.to ||
                       b.from == bond.from ||
                       b.to == bond.from
                   )
                    push!(bonds_queue, b)
                end
            end
        end
    end

    function to_SMILES(atom_idx)
        atom = molecule.atoms[atom_idx]
        atom_name_str = atom.is_aromatic ? lowercase(atom.name) : atom.name
        result = atom_name_str

        if haskey(ringbonds, atom_idx)
            for ringbond in ringbonds[atom_idx]
                ring_num = ringbond[1]
                bond_str = ringbond[2]

                # Prepend '%' if the ring number is double-digit
                if ring_num > 9
                    result *= bond_str * "%" * string(ring_num)
                else
                    result *= bond_str * string(ring_num)
                end
            end
        end

        if !haskey(adjacency, atom_idx)
            return result
        end
        for (i, neighbour) in enumerate(adjacency[atom_idx])
            bond = neighbour[2]
            if i == length(adjacency[atom_idx])
                result *= bond * to_SMILES(neighbour[1])
            else
                result *= "(" * bond * to_SMILES(neighbour[1]) * ")"
            end
        end

        return result
    end

    return to_SMILES(1)
end

# function to_SMILES(molecule::Molecule)::String
#     if isempty(molecule.atoms)
#         return ""
#     end

#     # Create data structures for traversal
#     adjacency = Dict{Int, Vector{Tuple{Int, String}}}()
#     ringbonds = Dict{Int, Vector{Tuple{Int, String}}}()
#     visited = falses(length(molecule.atoms))
#     ring_digit = 1

#     # Build complete adjacency list from bonds
#     complete_adj = Dict{Int, Vector{Tuple{Int, BondType}}}()
#     for bond in molecule.bonds
#         if !haskey(complete_adj, bond.from)
#             complete_adj[bond.from] = []
#         end
#         if !haskey(complete_adj, bond.to)
#             complete_adj[bond.to] = []
#         end
#         push!(complete_adj[bond.from], (bond.to, bond.bond_type))
#         push!(complete_adj[bond.to], (bond.from, bond.bond_type))
#     end

#     # DFS function to build the traversal-based adjacency map
#     function dfs_build_adjacency(atom_idx, parent_idx=0)
#         visited[atom_idx] = true

#         if !haskey(adjacency, atom_idx)
#             adjacency[atom_idx] = []
#         end

#         # Process all neighbors
#         if haskey(complete_adj, atom_idx)
#             for (neighbor, bond_type) in complete_adj[atom_idx]
#                 # Skip parent we came from
#                 if neighbor == parent_idx
#                     continue
#                 end

#                 bond_str = to_string(bond_type)

#                 if visited[neighbor]
#                     # We found a ring closure
#                     if !haskey(ringbonds, atom_idx)
#                         ringbonds[atom_idx] = []
#                     end
#                     if !haskey(ringbonds, neighbor)
#                         ringbonds[neighbor] = []
#                     end
#                     push!(ringbonds[atom_idx], (ring_digit, bond_str))
#                     push!(ringbonds[neighbor], (ring_digit, bond_str))
#                     ring_digit += 1
#                 else
#                     # Normal bond, add to adjacency and continue DFS
#                     push!(adjacency[atom_idx], (neighbor, bond_str))
#                     dfs_build_adjacency(neighbor, atom_idx)
#                 end
#             end
#         end
#     end

#     # Start DFS from atom 1
#     dfs_build_adjacency(1)

#     # Recursive function to build the SMILES string (unchanged)
#     function to_SMILES(atom_idx)
#         result = "[" * molecule.atoms[atom_idx].name * "]"

#         if haskey(ringbonds, atom_idx)
#             for ringbond in ringbonds[atom_idx]
#                 result *= ringbond[2] * string(ringbond[1])
#             end
#         end

#         if !haskey(adjacency, atom_idx)
#             return result
#         end

#         for (i, neighbour) in enumerate(adjacency[atom_idx])
#             bond = neighbour[2]
#             if i == length(adjacency[atom_idx])
#                 result *= bond * to_SMILES(neighbour[1])
#             else
#                 result *= "(" * bond * to_SMILES(neighbour[1]) * ")"
#             end
#         end

#         return result
#     end

#     return to_SMILES(1)
# end
