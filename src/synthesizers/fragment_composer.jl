function is_bond(character::Char)::Bool
    character == ':' || character == '-' || character == '=' || character == '≡'
end

function make_smiles_custom_explicit(smiles::String)::String
    smiles_chars = collect(replace(smiles, '#' => '≡'))

    # Store preceding bond types of digits
    digit_bonds_by_index = Dict{Int, Char}()
    digit_open = Dict{Char, Tuple{Int, Char, Char}}()
    active_digits = Set{Int}()
    in_bracket = false
    for (i, c) in enumerate(smiles_chars)
        if c == '['
            in_bracket = true
        elseif c == ']'
            in_bracket = false
        elseif isdigit(c) && !in_bracket
            bond = (i > 1 && is_bond(smiles_chars[i - 1])) ? smiles_chars[i - 1] : '-'
            if haskey(digit_open, c)
                first_i, first_bond, mapped_digit = digit_open[c]
                resolved_bond = '-'
                if bond != '-'
                    resolved_bond = bond
                elseif first_bond != '-'
                    resolved_bond = first_bond
                end
                digit_bonds_by_index[first_i] = resolved_bond
                digit_bonds_by_index[i] = resolved_bond
                smiles_chars[first_i] = mapped_digit
                smiles_chars[i] = mapped_digit

                # Release the mapped digit
                mapped_digit_val = Int(mapped_digit) - Int('0')
                delete!(active_digits, mapped_digit_val)

                delete!(digit_open, c)
            else
                # Find the smallest free digit starting from 1
                d = 1
                while d in active_digits
                    d += 1
                end
                push!(active_digits, d)

                mapped_digit = Char('0' + d)
                digit_open[c] = (i, bond, mapped_digit)
            end
        end
    end

    # Add missing square brackets and digit bonds
    result_smile = ""
    in_bracket = false
    i = 1
    while i <= length(smiles_chars)
        c = smiles_chars[i]
        if c == '['
            in_bracket = true
            result_smile *= string(c)
            i += 1
            continue
        elseif c == ']'
            in_bracket = false
            result_smile *= string(c)
            i += 1
            continue
        end

        if !in_bracket && isletter(c)
            # Check for two-letter atoms like Cl, Br
            if i < length(smiles_chars) && isletter(smiles_chars[i + 1])
                result_smile *= "[" * string(smiles_chars[i], smiles_chars[i + 1]) * "]"
                i += 2
            else
                result_smile *= "[" * string(c) * "]"
                i += 1
            end
        else
            result_smile *= string(c)
            i += 1
        end

        if i <= length(smiles_chars) && isdigit(smiles_chars[i]) &&
           last(result_smile) == ']'
            result_smile *= string(get(digit_bonds_by_index, i, '-')) *
                            string(smiles_chars[i])
            i += 1
        end
    end
    return result_smile
end

function make_fragment_custom_explicit(smiles::String)::Tuple{Int, Expr, Expr}
    result_smile = make_smiles_custom_explicit(smiles)

    # Extract entry digit
    entry_m = match(r"^\[(\d+)\*\]", result_smile)
    entry_digit = entry_m === nothing ? "X" : entry_m.captures[1]

    result_smile = replace(result_smile, r"^\[\d+\*\][-=≡]?" => "")
    result_smile = replace(
        result_smile, r"\([-=≡]?\[(\d+)\*\]\)" =>
            s -> "\" * fragment_" *
                 match(r"\d+", s).match *
                 "_exit * \"")
    result_smile = replace(result_smile,
        r"[-=≡]?\[(\d+)\*\]" =>
            s -> "\" * fragment_" * match(r"\d+", s).match *
                 "_exit * \"")

    entry_rule = "fragment_" * entry_digit * "_entry = \"" * result_smile * "\""
    entry_rule = replace(entry_rule, " * \"\"" => "")

    starting_smile = replace(result_smile,
        r"\](?:[-=≡]?\d+)*" => s -> s * "\" * fragment_" * entry_digit * "_exit * \"";
        count = 1)
    starting_rule = "starting_fragment = \"" * starting_smile * "\""
    starting_rule = replace(starting_rule, " * \"\"" => "")

    return (parse(Int, entry_digit), Meta.parse(entry_rule), Meta.parse(starting_rule))
end

# RDKit represents connections for all dummy atoms as single bonds
# However, the L7 BRICS environment has a fixed double bond
function make_seventh_brics_double_bond_explicit(smiles::String)::String
    s = replace(smiles, r"(?<=[A-Za-z\](])\[7\*\]|-\[7\*\]" => "=[7*]")
    s = replace(s, r"^\[7\*\](?:-|(?=[A-Za-z]))" => "[7*]=")
    
    return s
end

function make_smiles_rdkit_explicit(smiles::String)::String
    fixed_smiles = make_seventh_brics_double_bond_explicit(smiles)
    return MoleculeFlow.mol_to_smiles(
        MoleculeFlow.add_hs(MoleculeFlow.mol_from_smiles(fixed_smiles));
        kekule_smiles = true, all_bonds_explicit = true)
end

function has_connection_points(frag_smiles::String)::Bool
    '*' in frag_smiles
end

function parse_molecule_to_fragment_rules(mol_smiles::String)::Tuple{
        OrderedDict{Int, Vector{Expr}}, Vector{Expr}}
    brics_smiles = MoleculeFlow.brics_decompose(
        MoleculeFlow.mol_from_smiles(mol_smiles); min_fragment_size = 2)
    ismissing(brics_smiles) && return (OrderedDict{Int, Vector{Expr}}(), Expr[])
    filter!(has_connection_points, brics_smiles)
    map!(make_smiles_rdkit_explicit, brics_smiles)
    tuples = map(x -> make_fragment_custom_explicit(x), brics_smiles)
    fragment_rules = OrderedDict{Int, Vector{Expr}}()
    starting_fragments = Expr[]
    for (id, entry_rule, starting_rule) in tuples
        if !haskey(fragment_rules, id)
            fragment_rules[id] = Vector{Expr}()
        end
        if entry_rule ∉ fragment_rules[id]
            push!(fragment_rules[id], entry_rule)
        end
        if starting_rule ∉ starting_fragments
            push!(starting_fragments, starting_rule)
        end
    end
    return (fragment_rules, starting_fragments)
end
