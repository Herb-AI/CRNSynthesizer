using MoleculeFlow

# Create a molecule
mol = mol_from_smiles("C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N")
# BRICS fragmentation
brics_frags = brics_decompose(mol)

function parse_bond_type(bond::Bond)
    bond_type = get_bond_type(bond)

    bond_type == "SINGLE" && return "-"
    bond_type == "DOUBLE" && return "="
    bond_type == "TRIPLE" && return "≡"

    throw(ArgumentError(string("Unhandled bond type: ", bond_type)))
end

function is_dummy(atom::Atom)
    get_symbol(atom) == "*"
end

function parse_fragment_to_rule(frag_smiles::String)
    # Add explicit hydrogens for valency constraints
    bric_mol = add_hs(mol_from_smiles(frag_smiles))
    atoms = get_atoms(bric_mol)

    !any(map(is_dummy, atoms)) && return

    # Add explicit edge types instead of aromatic bond types for valency constraints
    kekulize!(bric_mol; clear_aromatic_flags = true)

    visited = falses(length(atoms))
    ring_edges = Set{Tuple{Int, Int}}()
    non_ring_edge_neighbors = [Vector{Int}() for _ in 1:length(atoms)]
    parsed_bond_types = [Vector{String}() for _ in 1:length(atoms)]
    parsed_ringbonds = ["" for _ in 1:length(atoms)]
    ring_counter = 0

    function identify_ring_edges(atom_idx::Int, parent_idx::Int)
        visited[atom_idx] = true
        bonds = get_bonds_from_atom(bric_mol, atom_idx)
        for bond in bonds
            begin_idx = get_begin_atom_idx(bond)
            end_idx = get_end_atom_idx(bond)
            neighbor_idx = begin_idx == atom_idx ? end_idx : begin_idx
            if neighbor_idx == parent_idx
                continue
            end
            if visited[neighbor_idx] # Cycle detected
                if neighbor_idx < atom_idx
                    ring_counter += 1
                    # Mark edge as ring closure
                    push!(ring_edges, (atom_idx, neighbor_idx))
                    push!(ring_edges, (neighbor_idx, atom_idx))
                    # Parse the ringbond for later
                    parsed_ringbond = string(parse_bond_type(bond), ring_counter)
                    parsed_ringbonds[atom_idx] *= parsed_ringbond
                    parsed_ringbonds[neighbor_idx] *= parsed_ringbond
                end
            else
                identify_ring_edges(neighbor_idx, atom_idx)
                if !((atom_idx, neighbor_idx) in ring_edges)
                    push!(non_ring_edge_neighbors[atom_idx], neighbor_idx)
                    push!(parsed_bond_types[atom_idx], parse_bond_type(bond))
                end
            end
        end
    end

    identify_ring_edges(2, 0)

    fill!(visited, false)
    fragment_rule = ""
    function write_fragment_rule_dfs(atom_idx::Int)
        visited[atom_idx] = true

        atom_symbol = get_symbol(atoms[atom_idx])
        fragment_rule *= string('[', atom_symbol, ']', parsed_ringbonds[atom_idx])

        neighbors_indices = non_ring_edge_neighbors[atom_idx]
        wrote_dummy = false
        for (i, other_atom_idx) in enumerate(neighbors_indices)
            if is_dummy(atoms[other_atom_idx])
                if !wrote_dummy
                    wrote_dummy = true
                    fragment_rule *= " * ringbonds * branches * "
                end
                continue
            end

            if visited[other_atom_idx]
                continue
            end

            if i < length(neighbors_indices)
                fragment_rule *= "("
            end

            fragment_rule *= parsed_bond_types[atom_idx][i]
            write_fragment_rule_dfs(other_atom_idx)

            if i < length(neighbors_indices)
                fragment_rule *= ")"
            end
        end
    end

    write_fragment_rule_dfs(2)
    println(fragment_rule)
end

function list_bonds(frag_smiles::String)
    bric_mol = add_hs(mol_from_smiles(frag_smiles))
    kekulize!(bric_mol; clear_aromatic_flags = true)
    atoms = get_atoms(bric_mol)
    for i in eachindex(atoms)
        bonds = get_bonds_from_atom(bric_mol, i)
        if !isempty(bonds)
            atom_symbol = get_symbol(get_atom(bric_mol, i))
            println("\nAtom $i ($atom_symbol) bonds:")

            for bond in bonds
                bond_type = get_bond_type(bond)
                begin_idx = get_begin_atom_idx(bond)
                end_idx = get_end_atom_idx(bond)
                is_aromatic_bond = is_aromatic(bond)
                is_in_ring_bond = is_in_ring(bond)

                other_atom_idx = begin_idx == i ? end_idx : begin_idx
                other_symbol = get_symbol(get_atom(bric_mol, other_atom_idx))

                println("to Atom $other_atom_idx ($other_symbol): $bond_type" *
                        (is_aromatic_bond ? " (aromatic)" : "") *
                        (is_in_ring_bond ? " (in ring)" : ""))
            end
        end
    end
end

for frag_smiles in brics_frags
    println("\nProcessing fragment: ", frag_smiles)
    list_bonds(frag_smiles)
    parse_fragment_to_rule(frag_smiles)
end