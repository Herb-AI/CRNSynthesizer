function bu_base_grammar(atoms::Vector{Atom})
    grammar = @csgrammar begin
        molecule = atom * ringbonds * branches

        structure = atom * ringbonds * branches

        branch = "(" * bond * structure * ")"
        branches = ""
        branches = branch * branches

        ringbond = bond * digit
        ringbonds = ""
        ringbonds = ringbond * ringbonds

        digit = "1" | "2" # | "3" | "4" | "5" | "6" | "7" | "8" | "9"

        bond = "-" | "=" | "≡" # | "≣"

        #= This rule makes it so expansion of fragment_X_exit into fragment_X_entry
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
        fragment_16_exit = "-" * digit =#
    end
    for atom in atoms
        atom_str = "[" * string(atom) * "]"
        grammar = add_rule!(grammar, :(atom = $atom_str))
    end

    # Make the ringbonds list tail ended 
    addconstraint!(grammar, Ordered(RuleNode(8, [VarNode(:a), RuleNode(8, [VarNode(:b), VarNode(:c)])]), [:a, :b]))
    addconstraint!(grammar, Forbidden(RuleNode(8, [VarNode(:a), RuleNode(8, [VarNode(:a), VarNode(:c)])])))

    # Order the strcutures in connected branches
    # Order the bonds types for same structures in connected branches
    addconstraint!(grammar, Ordered(RuleNode(5, [RuleNode(3, [VarNode(:b1), VarNode(:s1)]), RuleNode(5, [RuleNode(3, [VarNode(:b2), VarNode(:s2)]), VarNode(:c)])]), [:s1, :s2]))
    addconstraint!(grammar, Ordered(RuleNode(5, [RuleNode(3, [VarNode(:b1), VarNode(:s)]), RuleNode(5, [RuleNode(3, [VarNode(:b2), VarNode(:s)]), VarNode(:c)])]), [:b1, :b2]))
    return grammar
end

function bu_SMILES_grammar(
        atoms::Vector{Atom}; settings::SynthesizerSettings = SynthesizerSettings(),
        fragment_rules::Dict{Int, Vector{Expr}} = Dict{Int, Vector{Expr}}(), starting_fragments::Vector{Expr} = Expr[]
)
    grammar = bu_base_grammar(atoms)

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

    return grammar
end
