function count_digits(program::AbstractRuleNode, grammar::AbstractGrammar)::Int
    type = grammar.types[get_rule(program)]
    if type == :digit
        return 1
    end
    return sum(count_digits(child, grammar) for child in program.children; init = 0)
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

        :(structure * "=" * $fragment_X_entry) => begin
            structure_str = interpret_structure(program.children[1], grammar, min_digit)
            fragment_X_entry_str = interpret_fragment_X_entry(
                program.children[2], grammar, min_digit)
            return structure_str * "=" * fragment_X_entry_str
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
        :($fixed_bond * chain * ")") => begin
            return fixed_bond * interpret_chain(program.children[1], grammar, min_digit) * ")"
        end
        :($fixed_bond * digit) => begin
            return fixed_bond * interpret_digit(program.children[1], grammar)
        end
        :("(" * special_bond * $fragment_X_entry * ")") => begin
            return "(-" *
                   interpret_fragment_X_entry(program.children[2], grammar, min_digit) *
                   ")"
        end
        :("(" * special_double_bond * $fragment_X_entry * ")") => begin
            return "(=" *
                   interpret_fragment_X_entry(program.children[2], grammar, min_digit) *
                   ")"
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
