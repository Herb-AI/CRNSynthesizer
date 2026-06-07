# TODO: move back to immutable struct when possible
mutable struct GrammarData
    grammar_to_available_atom_connections::Dict{Int, Int}
    grammar_to_bond_connections::Dict{Int, Int}
    ringbond_digit_to_grammar::Dict{Int, Int}
    bond_to_grammar::Dict{String, Int}
end

const bond_orders = Dict(
    "-" => 1,
    "=" => 2,
    "≡" => 3,
    "≣" => 4
)

function update_grammar_data!(
        grammar_data::GrammarData, mapping::AbstractDict{<:Integer, <:Integer}
)
    new_grammar_to_available_atom_connections = Dict{Int, Int}()
    new_grammar_to_bond_connections = Dict{Int, Int}()
    new_ringbond_digit_to_grammar = Dict{Int, Int}()

    for (key, value) in grammar_data.grammar_to_available_atom_connections
        new_key = get(mapping, key, key)
        new_grammar_to_available_atom_connections[new_key] = value
    end

    for (key, value) in grammar_data.grammar_to_bond_connections
        new_key = get(mapping, key, key)
        new_grammar_to_bond_connections[new_key] = value
    end

    for (key, value) in grammar_data.ringbond_digit_to_grammar
        new_value = get(mapping, value, value)
        new_ringbond_digit_to_grammar[key] = new_value
    end

    grammar_data.grammar_to_available_atom_connections = new_grammar_to_available_atom_connections
    grammar_data.grammar_to_bond_connections = new_grammar_to_bond_connections
    grammar_data.ringbond_digit_to_grammar = new_ringbond_digit_to_grammar
end

function digit_to_grammar(grammar_data::GrammarData, digit::Int)::Int
    @assert haskey(grammar_data.ringbond_digit_to_grammar, digit) "key $digit not found in ringbond_digit_to_grammar $(grammar_data.ringbond_digit_to_grammar)"
    return grammar_data.ringbond_digit_to_grammar[digit]
end

function bond_to_grammar(grammar_data::GrammarData, bond::String)::Int
    @assert haskey(grammar_data.bond_to_grammar, bond) "key $bond not found in bond_to_grammar $(grammar_data.bond_to_grammar)"
    return grammar_data.bond_to_grammar[bond]
end

function grammar_to_atom_connections(grammar_data::GrammarData, rule::Int)::Int
    @assert haskey(grammar_data.grammar_to_available_atom_connections, rule) "key $rule not found in grammar_to_available_atom_connections $(grammar_data.grammar_to_available_atom_connections)"
    return grammar_data.grammar_to_available_atom_connections[rule]
end

function grammar_to_bond_connections(grammar_data::GrammarData, rule::Int)::Int
    @assert haskey(grammar_data.grammar_to_bond_connections, rule) "key $rule not found in grammar_to_bond_connections $(grammar_data.grammar_to_bond_connections)"
    return grammar_data.grammar_to_bond_connections[rule]
end

function generate_digit_bond_to_grammar(grammar)
    digit_dict = Dict{Int, Int}()
    bond_dict = Dict{String, Int}()

    # Process the grammar rules
    for (i, rule) in enumerate(grammar.rules)
        if grammar.types[i] == :digit
            for (digit, val) in digits
                if rule == digit
                    digit_dict[val] = i
                    break
                end
            end
        elseif grammar.types[i] == :bond
            for (bond, _) in bond_orders
                if rule == bond
                    bond_dict[bond] = i
                    break
                end
            end
        end
    end

    return digit_dict, bond_dict
end

function generate_atom_bond_dicts(grammar, atom_valences::OrderedDict{String, Int})
    atom_dict = Dict{Int, Int}()
    bond_dict = Dict{Int, Int}()

    # Process the grammar rules
    for (i, rule) in enumerate(grammar.rules)
        rule_str = string(rule)

        # Remove quotes if present in rule string
        clean_rule = replace(rule_str, "\"" => "")

        if grammar.types[i] == :atom
            for (atom, valence) in atom_valences
                if clean_rule == atom
                    atom_dict[i] = valence
                    break
                end
            end
        elseif grammar.types[i] == :bond
            for (bond, order) in bond_orders
                if clean_rule == bond
                    bond_dict[i] = order
                    break
                end
            end
        elseif startswith(string(grammar.types[i]), "fragment_7")
            bond_dict[i] = 2
        elseif startswith(string(grammar.types[i]), "fragment_")
            bond_dict[i] = 1
        end
    end

    return atom_dict, bond_dict
end
