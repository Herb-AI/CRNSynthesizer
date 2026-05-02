using MoleculeFlow

# Create a molecule
mol = add_hs(mol_from_smiles("c1ccccc1C(=O)NC"))

# BRICS fragmentation
brics_frags = brics_decompose(mol)

function parse_fragment_to_rule(frag_smiles::String)
    bric_mol = add_hs(mol_from_smiles(frag_smiles))
    atoms = get_atoms(bric_mol)

    fragment_rule = ""
    function fragment_dfs(atom_idx::Int)
        bonds = get_bonds_from_atom(bric_mol, atom_idx)
        if !isempty(bonds)
            atom_symbol = get_symbol(get_atom(bric_mol, atom_idx))
            fragment_rule *= string('[', atom_symbol, ']')
            wrote_dummy = false
            for (bond_idx, bond) in enumerate(bonds)
                bond_type = get_bond_type(bond)
                begin_idx = get_begin_atom_idx(bond)
                end_idx = get_end_atom_idx(bond)
                is_in_ring_bond = is_in_ring(bond)
                other_atom_idx = begin_idx == atom_idx ? end_idx : begin_idx
                other_symbol = get_symbol(get_atom(bric_mol, other_atom_idx))
                if other_symbol == "*"
                    if !wrote_dummy
                        wrote_dummy = true
                        fragment_rule *= " * branches * "
                    end
                    continue
                end
                if other_atom_idx < atom_idx
                    continue
                end
                if bond_idx < length(bonds)
                    fragment_rule *= "("
                end
                if bond_type == "SINGLE"
                    fragment_rule *= "-"
                elseif bond_type == "DOUBLE"
                    fragment_rule *= "="
                elseif bond_type == "TRIPLE"
                    fragment_rule *= "≡"
                end
                fragment_dfs(other_atom_idx)
                if bond_idx < length(bonds)
                    fragment_rule *= ")"
                end
            end
        end
    end

    fragment_dfs(2)
    println(fragment_rule)
end

function list_bonds(frag_smiles::String)
    bric_mol = add_hs(mol_from_smiles(frag_smiles))
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