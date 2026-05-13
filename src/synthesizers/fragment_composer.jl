function is_bond(character::Char)::Bool
    character == '-' || character == '=' || character == '≡'
end

function make_fragment_custom_explicit(frag_smile::String)::Tuple{Expr, Expr}
    frag_smile = replace(frag_smile, '#' => '≡')

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

    result_smile = replace(result_smile, r"^\[\d+\*\][-=≡]?" => "")
    result_smile = replace(
        result_smile, r"\([-=≡]?\[\d+\*\]\)" => "\" * fragment_X_exit * \"")
    result_smile = replace(result_smile, r"[-=≡]?\[\d+\*\]" => "\" * fragment_X_exit * \"")

    entry_rule = "fragment_X_entry = \"" * result_smile * "\""
    entry_rule = replace(entry_rule, " * \"\"" => "")

    starting_smile = replace(result_smile, r"\](?:[-=≡]?\d+)*" => s -> s * "\" * fragment_X_exit * \""; count=1)
    starting_rule = "starting_fragment = \"" * starting_smile * "\""
    starting_rule = replace(starting_rule, " * \"\"" => "")

    return (Meta.parse(entry_rule), Meta.parse(starting_rule))
end

function make_fragment_rdkit_explicit(frag_smiles::String)::String
    MoleculeFlow.mol_to_smiles(
        MoleculeFlow.add_hs(MoleculeFlow.mol_from_smiles(frag_smiles));
        kekule_smiles = true, all_bonds_explicit = true)
end

function has_connection_points(frag_smiles::String)::Bool
    '*' in frag_smiles
end

function parse_molecule_to_fragment_rules(mol_smiles::String)::Tuple{Vector{Expr}, Vector{Expr}}
    brics_smiles = MoleculeFlow.brics_decompose(
        MoleculeFlow.mol_from_smiles(mol_smiles); min_fragment_size = 2)
    filter!(has_connection_points, brics_smiles)
    map!(make_fragment_rdkit_explicit, brics_smiles)
    tuples = map(x -> make_fragment_custom_explicit(x), brics_smiles)
    return (first.(tuples), last.(tuples))
end