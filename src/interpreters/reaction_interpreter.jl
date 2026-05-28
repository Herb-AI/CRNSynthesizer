function interpret_reaction(program::AbstractRuleNode, grammar::AbstractGrammar)::Reaction
    rule = grammar.rules[get_rule(program)]

    if rule isa Reaction
        return rule
    end

    @match rule begin
        :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules))) => begin
            inputs_list = interpret_molecule_list(program.children[1], grammar)
            prefilled_inputs = interpret_molecule_list(program.children[2], grammar)
            outputs_list = interpret_molecule_list(program.children[3], grammar)
            prefilled_outputs = interpret_molecule_list(program.children[4], grammar)
            return Reaction(vcat(inputs_list, prefilled_inputs), vcat(outputs_list, prefilled_outputs))
        end
        :(Reaction(molecule_list, molecule_list)) => begin
            inputs_list = interpret_molecule_list(program.children[1], grammar)
            outputs_list = interpret_molecule_list(program.children[2], grammar)
            return Reaction(inputs_list, outputs_list)
        end
        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_molecule_list(
        program::AbstractRuleNode, grammar::AbstractGrammar
)::Vector{Molecule}
    rule = grammar.rules[get_rule(program)]

    if rule isa Vector{Molecule}
        return rule
    end

    @match rule begin
        :(Vector{Molecule}()) => begin
            return Vector{Molecule}()
        end

        :(vcat([molecule],
            molecule_list)) => begin
            molecule::Molecule = interpret_molecule(program.children[1], grammar)
            molecule_list::Vector{Molecule} = interpret_molecule_list(
                program.children[2], grammar
            )
            return vcat(molecule, molecule_list)
        end

        :(vcat([required_molecule],
            molecule_list)) => begin
            required_molecule::Molecule = interpret_molecule(program.children[1], grammar)
            molecule_list::Vector{Molecule} = interpret_molecule_list(
                program.children[2], grammar
            )
            return vcat(required_molecule, molecule_list)
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end
