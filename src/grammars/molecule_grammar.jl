function SMILES_combine_chain(bond, structure, chain)
    structure * bond * chain
end

function SMILES_grammar(
        atoms::Vector{Atom}; settings::SynthesizerSettings = SynthesizerSettings(),
        fragment_rules::Vector{Expr} = Expr[], starting_fragments::Vector{Expr} = Expr[]
)
    grammar = @csgrammar begin
        molecule = starting_fragment
        molecule = chain

        chain = structure * "-" * fragment_X_entry
        chain = atom * ringbonds
        chain = SMILES_combine_chain(bond, structure, chain)
        # chain = structure * bond * chain

        structure = atom * ringbonds * branches

        # Need to force brackets around non-ringbonds, as brackets around ringbonds are invalid syntax
        fragment_X_exit = "(" * special_bond * fragment_X_entry * ")"
        fragment_X_exit = "(-" * chain * ")"
        fragment_X_exit = "-" * digit

        # fragment_X_entry = "[N]-3-[C](-[H])(-[H])-[C](-[H])(-[H])-[C](-[H])(-[H])-[C]" * fragment_X_exit * "(-[H])-[C]-3(-[H])-[H]"

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
    end

    @assert length(grammar.rules)==21 "Length of static molecule grammar changed. Ensure the indices of static rules in the constraints below are updated accordingly."

    if isempty(starting_fragments)
        addconstraint!(grammar, Forbidden((@c_rulenode 1)))
    else
        for rule in starting_fragments
            grammar = add_rule!(grammar, rule)
        end
    end

    if isempty(fragment_rules)
        addconstraint!(grammar, Forbidden((@c_rulenode 3)))
        addconstraint!(grammar, Forbidden((@c_rulenode 7)))
    else
        for rule in fragment_rules
            grammar = add_rule!(grammar, rule)
        end
    end

    for atom in atoms
        atom_str = "[" * string(atom) * "]"
        grammar = add_rule!(grammar, :(atom = $atom_str))
    end

    if !(
        haskey(settings.options, :disable_valid_smiles) &&
        settings.options[:disable_valid_smiles]
    )
        atom_dict, bond_dict = generate_atom_bond_dicts(grammar)
        digit_to_grammar, bond_to_grammar = generate_digit_bond_to_grammar(grammar)
        grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar, bond_to_grammar)
        addconstraint!(grammar, ValidSMILES(grammar_data))
    end

    # Make the ringbonds list tail ended 
    addconstraint!(grammar, Ordered((@c_rulenode 15{a, 15{b, c}}), [:a, :b]))
    addconstraint!(grammar, Forbidden((@c_rulenode 15{a, 15{a, c}})))

    # Make the branches list tail ended
    addconstraint!(grammar, Ordered((@c_rulenode 12{a, 12{b, c}}), [:a, :b]))
    return grammar
end
