function count_digits(program::AbstractRuleNode, grammar::AbstractGrammar)::Int
    type = grammar.types[get_rule(program)]
    if type == :digit
        return 1
    end
    return sum(count_digits(child, grammar) for child in program.children; init = 0)
end

function rdkit_canonicalize(smiles::String)::Union{String, Nothing}
    mol = RDKitMinimalLib.get_mol(replace(smiles, "≡" => "#"))
    return !isnothing(mol) ? RDKitMinimalLib.get_smiles(mol, Dict{String, Any}("allBondsExplicit" => true, "allHsExplicit" => true)) : nothing
end

function replace_unpaired_digits(smiles::String)
    ringbonds = [m.match for m in eachmatch(r"%?\d+", smiles)]
    counts = Dict{String, Int}()
    for r in ringbonds
        counts[r] = get(counts, r, 0) + 1
    end
    
    new_smiles = smiles
    for (r, c) in counts
        if c == 1
            pat = Regex("[-=≡]$r")
            new_smiles = replace(new_smiles, pat => s -> "(" * replace(s, r => "*") * ")")
        end
    end
    return new_smiles
end

function has_only_first_atom_incomplete_valency(candidate::Molecule)
    atom_valences = Dict("O" => 2, "H" => 1, "C" => 4, "N" => 3)
    bond_orders = Dict(single => 1, double => 2, triple => 3, quadruple => 4)

    atoms = candidate.atoms
    bonds = candidate.bonds

    for (atom_index, atom) in enumerate(atoms)
        atom_str = atom.name
        if !haskey(atom_valences, atom_str)
            return false
        end

        valence = atom_valences[atom_str]
        connected_bonds = filter(b -> b.from == atom_index || b.to == atom_index, bonds)
        total_bond_order = isempty(connected_bonds) ? 0 : sum(bond_orders[b.bond_type] for b in connected_bonds)

        if atom_index == 1
            if total_bond_order >= valence
                return false
            end
        else
            if total_bond_order != valence
                return false
            end
        end
    end

    return true
end

function are_ringbond_digits_strictly_increasing(ringbonds_str::String)::Bool
    if isempty(ringbonds_str)
        return true
    end
    digit_matches = [parse(Int, m.match) for m in eachmatch(r"\d+", ringbonds_str)]
    for i in 2:length(digit_matches)
        if digit_matches[i] <= digit_matches[i-1]
            return false
        end
    end
    return true
end

function interpret_partial_canonical(
        program::AbstractRuleNode, grammar::AbstractGrammar)::Union{String, Nothing}
    min_digit = Ref{Int64}(count_digits(program, grammar) ÷ 2 + 1)
    type = grammar.types[get_rule(program)]

    result = @match type begin
        :molecule => nothing
        :starting_fragment_grammar => begin
            smiles = bu_interpret_fragment_X_entry(program, grammar, min_digit)
            mol = from_SMILES(smiles)
            is_valid(mol) ? mol.canonical_smiles : nothing
        end 
        :structure => begin
            smiles = bu_interpret_structure(program, grammar, min_digit)
            safe_smiles = rdkit_canonicalize(replace_unpaired_digits(smiles))
            if !isnothing(safe_smiles)
                mol = from_SMILES(smiles)
                if has_only_first_atom_incomplete_valency(mol)
                    return safe_smiles * "structure"
                end
            end
            nothing
        end 
        :atom => interpret_atom(program, grammar)
        :digit => interpret_digit(program, grammar)
        :bond => interpret_bond(program, grammar)
        :special_bond => "special_bond"
        :ringbond => interpret_ringbond(program, grammar)
        :ringbonds => begin
            ringbonds = interpret_ringbonds(program, grammar)
            if !are_ringbond_digits_strictly_increasing(ringbonds)
              return nothing
            end
            ringbonds == "" ? "empty_ringbonds" : ringbonds * "ringbonds"
        end
        :branch => begin
            smiles = bu_interpret_branch(program, grammar, min_digit)
            safe_smiles = rdkit_canonicalize("*" * replace_unpaired_digits(smiles))
            !isnothing(safe_smiles) ? safe_smiles * "branch" : nothing
        end 
        :branches => begin
            branches = bu_interpret_branches(program, grammar, min_digit)
            if branches == ""
                "empty_branches"
            else
                safe_branches = rdkit_canonicalize("*" * replace_unpaired_digits(branches))
                !isnothing(safe_branches) ? safe_branches * "branches" : nothing
            end
        end
        _ => begin
            type_str = string(type)
            m = match(r"\d+", type_str)
            num_str = m === nothing ? "" : m.match
            if endswith(type_str, "_entry")
                smiles = interpret_fragment_X_entry(program, grammar, min_digit)
                safe_smiles = replace_unpaired_digits(smiles)
                rdkit_validate("*-" * safe_smiles) ? safe_smiles * "entry" * num_str : nothing
            elseif endswith(type_str, "_exit")
                raw_smiles = interpret_fragment_X_exit(program, grammar, min_digit)
                safe_smiles = replace_unpaired_digits(raw_smiles)
                if raw_smiles[1] == '('
                    rdkit_validate("*" * safe_smiles) ? safe_smiles * "exit" * num_str : nothing
                else
                    safe_smiles * "exit" * num_str
                end
            else
                throw(ArgumentError("Unknown type for partial canonicalization: $type"))
            end
        end
    end

    return result
end

function bu_interpret_molecule(program::AbstractRuleNode, grammar::AbstractGrammar)::Molecule
    rule = grammar.rules[get_rule(program)]

    if rule isa Molecule
        return rule
    end

    min_digit = Ref{Int64}(count_digits(program, grammar) ÷ 2 + 1)

    @match rule begin
        :(atom * ringbonds * branches) => begin
            return from_SMILES(bu_interpret_structure(program, grammar, min_digit))
        end
        :(starting_fragment) => begin
            return from_SMILES(bu_interpret_fragment_X_entry(program.children[1], grammar, min_digit))
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function bu_interpret_fragment_X_entry(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]
    digit_map = Dict{String, String}()

    if rule isa String
        create_digit_mapping(digit_map, rule, min_digit)
        return replace(rule, r"\d" => m -> digit_map[m])
    end

    result = ""
    child_count = 0
    for arg in rule.args[2:end]
        if arg isa String
            create_digit_mapping(digit_map, arg, min_digit)
        end
    end

    for arg in rule.args[2:end]
        if arg isa String
            mapped_arg = replace(arg, r"\d" => m -> digit_map[m])
            result *= mapped_arg
        else
            child_count += 1
            result *= interpret_fragment_X_exit(
                program.children[child_count], grammar, min_digit)
        end
    end
    return result
end

function bu_interpret_fragment_X_exit(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :("(-" * chain * ")") => begin
            return "(-" * interpret_chain(program.children[1], grammar, min_digit) * ")"
        end
        :("-" * digit) => begin
            return "-" * interpret_digit(program.children[1], grammar)
        end
        :("(" * special_bond * $fragment_X_entry * ")") => begin
            return "(-" *
                   interpret_fragment_X_entry(program.children[2], grammar, min_digit) * ")"
        end
        _ => throw(ArgumentError("Unknown fragment exit rule: $rule"))
    end
end

function bu_interpret_structure(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(atom * ringbonds *
        branches) => begin
            atom_str = interpret_atom(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            branches_str = bu_interpret_branches(program.children[3], grammar, min_digit)
            return atom_str * ringbonds_str * branches_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function bu_interpret_branches(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(branch *
        branches) => begin
            branch_str = bu_interpret_branch(program.children[1], grammar, min_digit)
            branches_str = bu_interpret_branches(program.children[2], grammar, min_digit)
            return branch_str * branches_str
        end

        :("") => begin
            return ""
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function bu_interpret_branch(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :("(" * bond * structure *
        ")") => begin
            bond_str = interpret_bond(program.children[1], grammar)
            structure_str = bu_interpret_structure(program.children[2], grammar, min_digit)
            return "(" * bond_str * structure_str * ")"
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end
