using MoleculeFlow

function is_bond(character::Char)
    character == '-' || character == '=' || character == '≡'
end

function make_fragment_custom_explicit(frag_smile::String; min_digit::Int = 1)
    frag_smile = replace(
        frag_smile, '#' => '≡', r"\d+" => d -> string(parse(Int, d) + min_digit - 1))

    # Store preceding bond types of digits
    digit_bonds = Dict{Char, Char}()
    for (i, c) in enumerate(frag_smile)
        if isdigit(c) && is_bond(frag_smile[i - 1])
            digit_bonds[c] = frag_smile[i - 1]
        end
    end

    # Add missing square brackets and digit bonds
    result_smile = ""
    in_bracket = false
    i = 1
    while i <= length(frag_smile)
        c = frag_smile[i]
        if c == '['
            in_bracket = true
            result_smile *= c
            i += 1
            continue
        elseif c == ']'
            in_bracket = false
            result_smile *= c
            i += 1
            continue
        end

        if !in_bracket && isletter(c)
            # Check for two-letter atoms like Cl, Br
            if i < length(frag_smile) && isletter(frag_smile[i + 1])
                result_smile *= "[" * frag_smile[i:(i + 1)] * "]"
                i += 2
            else
                result_smile *= "[" * c * "]"
                i += 1
            end
        else
            result_smile *= c
            i += 1
        end

        if i < length(frag_smile) && isdigit(frag_smile[i]) && last(result_smile) == ']'
            result_smile *= get(digit_bonds, frag_smile[i], '-') * frag_smile[i]
            i += 1
        end
    end

    result_smile
end

function make_fragment_rdkit_explicit(frag_smiles::String)
    mol_to_smiles(add_hs(mol_from_smiles(frag_smiles));
        kekule_smiles = true, all_bonds_explicit = true)
end

function has_connection_points(frag_smiles::String)
    '*' in frag_smiles
end

function parse_molecule_to_fragment_rules(mol_smiles::String)
    brics_smiles = brics_decompose(mol_from_smiles(mol_smiles))
    filter!(has_connection_points, brics_smiles)
    map!(make_fragment_custom_explicit ∘ make_fragment_rdkit_explicit, brics_smiles)
    println(brics_smiles)
end

# parse_molecule_to_fragment_rules("C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N")
