function count_digits(program::AbstractRuleNode, grammar::AbstractGrammar)::Int
    type = grammar.types[get_rule(program)]
    if type == :digit
        return 1
    end
    return sum(count_digits(child, grammar) for child in program.children; init = 0)
end

function rdkit_validate(smiles::String)::Bool
    mol = RDKitMinimalLib.get_mol(replace(smiles, "≡" => "#"))
    return !isnothing(mol)
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

function interpret_partial_canonical(
        program::AbstractRuleNode, grammar::AbstractGrammar)::Union{String, Nothing}
    min_digit = Ref{Int64}(count_digits(program, grammar) ÷ 2 + 1)
    type = grammar.types[get_rule(program)]

    result = @match type begin
        :molecule => begin
            mol = interpret_molecule(program, grammar)
            is_valid(mol) ? mol.canonical_smiles : nothing
        end
        :starting_fragment_grammar => begin
            smiles = interpret_fragment_X_entry(program, grammar, min_digit)
            mol = from_SMILES(smiles)
            is_valid(mol) ? mol.canonical_smiles : nothing
        end 
        :chain => begin
            smiles = interpret_chain(program, grammar, min_digit)
            safe_smiles = replace_unpaired_digits(smiles)
            rdkit_validate(safe_smiles) ? safe_smiles * "chain" : nothing
        end
        :structure => begin
            smiles = interpret_structure(program, grammar, min_digit)
            safe_smiles = replace_unpaired_digits(smiles)
            rdkit_validate(safe_smiles) ? safe_smiles * "structure" : nothing
        end 
        :atom => interpret_atom(program, grammar)
        :digit => interpret_digit(program, grammar)
        :bond => interpret_bond(program, grammar)
        :special_bond => "special_bond"
        :ringbond => interpret_ringbond(program, grammar)
        :ringbonds => begin
            ringbonds = interpret_ringbonds(program, grammar)
            ringbonds == "" ? "empty_ringbonds" : ringbonds * "ringbonds"
        end
        :branch => begin
            smiles = interpret_branch(program, grammar, min_digit)
            safe_smiles = replace_unpaired_digits(smiles)
            rdkit_validate("*" * safe_smiles) ? safe_smiles : nothing
        end 
        :branches => begin
            branches = interpret_branches(program, grammar, min_digit)
            if branches == ""
                "empty_branches"
            else
                safe_branches = replace_unpaired_digits(branches)
                rdkit_validate("*" * safe_branches) ? safe_branches * "branches" : nothing
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

function interpret_molecule(program::AbstractRuleNode, grammar::AbstractGrammar)::Molecule
    rule = grammar.rules[get_rule(program)]

    if rule isa Molecule
        return rule
    end

    min_digit = Ref{Int64}(count_digits(program, grammar) ÷ 2 + 1)

    @match rule begin
        :(chain) => begin
            return from_SMILES(interpret_chain(program.children[1], grammar, min_digit))
        end
        :(starting_fragment) => begin
            smiles = interpret_fragment_X_entry(program.children[1], grammar, min_digit)
            return from_SMILES(smiles)
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_chain(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(SMILES_combine_chain(bond,
        structure,
        chain)) => begin
            bond_str = interpret_bond(program.children[1], grammar)
            structure_str = interpret_structure(program.children[2], grammar, min_digit)
            chain_str = interpret_chain(program.children[3], grammar, min_digit)
            return structure_str * bond_str * chain_str
        end

        :(structure * bond *
        chain) => begin
            structure_str = interpret_structure(program.children[1], grammar, min_digit)
            bond_str = interpret_bond(program.children[2], grammar)
            chain_str = interpret_chain(program.children[3], grammar, min_digit)
            return structure_str * bond_str * chain_str
        end

        :(atom *
        ringbonds) => begin
            atom_str = interpret_atom(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            return atom_str * ringbonds_str
        end

        :(structure * "-" * $fragment_X_entry) => begin
            structure_str = interpret_structure(program.children[1], grammar, min_digit)
            fragment_X_entry_str = interpret_fragment_X_entry(
                program.children[2], grammar, min_digit)
            return structure_str * "-" * fragment_X_entry_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function create_digit_mapping(
        digit_map::Dict{String, String}, arg::String, min_digit::Ref{Int64})
    for m in eachmatch(r"\d", arg)
        d = m.match
        if !haskey(digit_map, d)
            digit_map[d] = min_digit[] >= 10 ? "%" * string(min_digit[]) :
                           string(min_digit[])
            min_digit[] += 1 # Increment the global digit counter
        end
    end
end

function interpret_fragment_X_entry(
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

function interpret_fragment_X_exit(
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

function interpret_bond(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_structure(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(atom * ringbonds *
        branches) => begin
            atom_str = interpret_atom(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            branches_str = interpret_branches(program.children[3], grammar, min_digit)
            return atom_str * ringbonds_str * branches_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_atom(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_ringbonds(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(ringbond *
        ringbonds) => begin
            ringbond_str = interpret_ringbond(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            return ringbond_str * ringbonds_str
        end

        :("") => begin
            return ""
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_ringbond(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(bond * digit) => begin
            bond_str = interpret_bond(program.children[1], grammar)
            digit_str = interpret_digit(program.children[2], grammar)
            return bond_str * digit_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_digit(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_branches(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(branch *
        branches) => begin
            branch_str = interpret_branch(program.children[1], grammar, min_digit)
            branches_str = interpret_branches(program.children[2], grammar, min_digit)
            return branch_str * branches_str
        end

        :("") => begin
            return ""
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_branch(
        program::AbstractRuleNode, grammar::AbstractGrammar, min_digit::Ref{Int64})::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :("(" * bond * chain *
        ")") => begin
            bond_str = interpret_bond(program.children[1], grammar)
            chain_str = interpret_chain(program.children[2], grammar, min_digit)
            return "(" * bond_str * chain_str * ")"
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end
