# RDKit's BRICS connection rules
# https://github.com/rdkit/rdkit-orig/blob/57058c886a49cc597b0c40641a28697ee3a57aee/Code/GraphMol/ChemTransforms/MolFragmenter.cpp#L205
#=// L1\n\
1 3 -;!@\n\
1 5 -;!@\n\
1 10 -;!@\n\
// L3 \n\
3 4 -;!@\n\
3 13 -;!@\n\
3 14 -;!@\n\
3 15 -;!@\n\
3 16 -;!@\n\
// L4\n\
4 5 -;!@\n\
4 11 -;!@\n\
// L5\n\
5 12 -;!@\n\
5 14 -;!@\n\
5 16 -;!@\n\
5 13 -;!@\n\
5 15 -;!@\n\
// L6\n\
6 13 -;!@\n\
6 14 -;!@\n\
6 15 -;!@\n\
6 16 -;!@\n\
// L7\n\
7 7 =;!@\n\
// L8\n\
8 9 -;!@\n\
8 10 -;!@\n\
8 13 -;!@\n\
8 14 -;!@\n\
8 15 -;!@\n\
8 16 -;!@\n\
// L9\n\
9 13 -;!@ // not in original paper\n\
9 14 -;!@ // not in original paper\n\
9 15 -;!@\n\
9 16 -;!@\n\
// L10\n\
10 13 -;!@\n\
10 14 -;!@\n\
10 15 -;!@\n\
10 16 -;!@\n\
// L11\n\
11 13 -;!@\n\
11 14 -;!@\n\
11 15 -;!@\n\
11 16 -;!@\n\
// L12\n\
// none left\n\
// L13\n\
13 14 -;!@\n\
13 15 -;!@\n\
13 16 -;!@\n\
// L14\n\
14 14 -;!@ // not in original paper\n\
14 15 -;!@\n\
14 16 -;!@\n\
// L15\n\
15 16 -;!@\n\
// L16\n\
16 16 -;!@ // not in original paper"=#

function SMILES_combine_chain(bond, structure, chain)
    structure * bond * chain
end

function base_grammar(atoms::Vector{String})
    grammar = @csgrammar begin
        molecule = chain

        chain = atom * ringbonds
        chain = SMILES_combine_chain(bond, structure, chain)
        # chain = structure * bond * chain

        structure = atom * ringbonds * branches

        branch = "(" * bond * chain * ")"
        branches = ""
        branches = branch * branches

        ringbond = bond * digit
        ringbonds = ""
        ringbonds = ringbond * ringbonds

        digit = "1" | "2" # | "3" | "4" | "5" | "6" | "7" | "8" | "9"

        bond = "-" | "=" | "≡" # | "≣"

        # This rule makes it so expansion of fragment_X_exit into fragment_X_entry
        # does not have uniform shape with expansion into a ringbond digit
        special_bond = "-"

        fragment_1_exit = "(-" * chain * ")"
        fragment_1_exit = "-" * digit
        fragment_3_exit = "(-" * chain * ")"
        fragment_3_exit = "-" * digit
        fragment_4_exit = "(-" * chain * ")"
        fragment_4_exit = "-" * digit
        fragment_5_exit = "(-" * chain * ")"
        fragment_5_exit = "-" * digit
        fragment_6_exit = "(-" * chain * ")"
        fragment_6_exit = "-" * digit
        fragment_7_exit = "(-" * chain * ")"
        fragment_7_exit = "-" * digit
        fragment_8_exit = "(-" * chain * ")"
        fragment_8_exit = "-" * digit
        fragment_9_exit = "(-" * chain * ")"
        fragment_9_exit = "-" * digit
        fragment_10_exit = "(-" * chain * ")"
        fragment_10_exit = "-" * digit
        fragment_11_exit = "(-" * chain * ")"
        fragment_11_exit = "-" * digit
        fragment_12_exit = "(-" * chain * ")"
        fragment_12_exit = "-" * digit
        fragment_13_exit = "(-" * chain * ")"
        fragment_13_exit = "-" * digit
        fragment_14_exit = "(-" * chain * ")"
        fragment_14_exit = "-" * digit
        fragment_15_exit = "(-" * chain * ")"
        fragment_15_exit = "-" * digit
        fragment_16_exit = "(-" * chain * ")"
        fragment_16_exit = "-" * digit
    end

    for atom in atoms
        atom_str = atom
        grammar = add_rule!(grammar, :(atom = $atom_str))
    end

    # Make the ringbonds list tail ended 
    addconstraint!(grammar, Ordered((@c_rulenode 10{a, 10{b, c}}), [:a, :b]))
    addconstraint!(grammar, Forbidden((@c_rulenode 10{a, 10{a, c}})))

    # Make the branches list tail ended
    addconstraint!(grammar, Ordered((@c_rulenode 7{a, 7{b, c}}), [:a, :b]))
    return grammar
end

function starting_fragment_grammar(starting_fragments::Vector{Expr} = Expr[])
    grammar = @csgrammar begin
        molecule = starting_fragment
    end
    for rule in starting_fragments
        grammar = add_rule!(grammar, rule)
    end
    return grammar
end

function fragment_X_grammar(id::Int, fragment_rules::Vector{Expr})
    @match id begin
        1 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_1_entry
                fragment_3_exit = "(" * special_bond * fragment_1_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_1_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_1_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        3 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_3_entry
                fragment_1_exit = "(" * special_bond * fragment_3_entry * ")"
                fragment_4_exit = "(" * special_bond * fragment_3_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_3_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_3_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_3_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_3_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        4 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_4_entry
                fragment_3_exit = "(" * special_bond * fragment_4_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_4_entry * ")"
                fragment_11_exit = "(" * special_bond * fragment_4_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        5 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_5_entry
                fragment_1_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_4_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_12_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_5_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_5_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        6 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_6_entry
                fragment_13_exit = "(" * special_bond * fragment_6_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_6_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_6_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_6_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        7 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_7_entry
                fragment_7_exit = "(" * special_bond * fragment_7_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        8 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_8_entry
                fragment_9_exit = "(" * special_bond * fragment_8_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_8_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_8_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_8_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_8_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_8_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        9 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_9_entry
                fragment_8_exit = "(" * special_bond * fragment_9_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_9_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_9_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_9_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_9_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        10 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_10_entry
                fragment_1_exit = "(" * special_bond * fragment_10_entry * ")"
                fragment_8_exit = "(" * special_bond * fragment_10_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_10_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_10_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_10_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_10_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        11 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_11_entry
                fragment_4_exit = "(" * special_bond * fragment_11_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_11_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_11_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_11_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_11_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        12 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_12_entry
                fragment_5_exit = "(" * special_bond * fragment_12_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        13 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_13_entry
                fragment_3_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_6_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_8_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_9_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_11_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_13_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_13_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        14 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_14_entry
                fragment_3_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_6_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_8_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_9_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_11_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_14_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_14_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        15 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_15_entry
                fragment_3_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_6_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_8_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_9_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_11_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_15_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_15_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end

        16 => begin
            grammar = @csgrammar begin
                chain = structure * "-" * fragment_16_entry
                fragment_3_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_5_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_6_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_8_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_9_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_10_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_11_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_13_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_14_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_15_exit = "(" * special_bond * fragment_16_entry * ")"
                fragment_16_exit = "(" * special_bond * fragment_16_entry * ")"
            end
            foreach(rule -> add_rule!(grammar, rule), fragment_rules)
            return grammar
        end
    end
end

function SMILES_grammar(
        atom_valences::OrderedDict{String, Int}; settings::SynthesizerSettings = SynthesizerSettings(),
        fragment_rules::OrderedDict{Int, Vector{Expr}} = OrderedDict{Int, Vector{Expr}}(), starting_fragments::Vector{Expr} = Expr[]
)
    grammar = base_grammar(collect(keys(atom_valences)))

    if !isempty(starting_fragments)
        new_grammar = starting_fragment_grammar(starting_fragments)
        merge_grammars!(new_grammar, grammar)
        grammar = new_grammar
    end

    for (id, rules) in fragment_rules
        if !isempty(rules)
            new_grammar = fragment_X_grammar(id, rules)
            merge_grammars!(new_grammar, grammar)
            grammar = new_grammar
        end
    end

    if !(
        haskey(settings.options, :disable_valid_smiles) &&
        settings.options[:disable_valid_smiles]
    )
        atom_dict, bond_dict = generate_atom_bond_dicts(grammar, atom_valences)
        digit_to_grammar, bond_to_grammar = generate_digit_bond_to_grammar(grammar)
        grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar, bond_to_grammar)
        addconstraint!(grammar, ValidSMILES(grammar_data, atom_valences))
    end

    return grammar
end
