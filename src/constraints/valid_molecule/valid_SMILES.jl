# TODO: check the place of this definition
# Define chemical properties
digits = Dict(
    "1" => 1, "2" => 2, "3" => 3, "4" => 4, "5" => 5, "6" => 6, "7" => 7, "8" => 8, "9" =>
        9
)

struct ValidSMILES <: AbstractGrammarConstraint
    grammar_data::GrammarData
    atom_valences::OrderedDict{String, Int}
end

function ValidSMILES(grammar::ContextSensitiveGrammar, atom_valences::OrderedDict{
        String, Int})
    atom_dict, bond_dict = generate_atom_bond_dicts(grammar, atom_valences)
    digit_to_grammar, bond_to_grammar = generate_digit_bond_to_grammar(grammar)
    grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar, bond_to_grammar)
    return ValidSMILES(grammar_data, atom_valences)
end

function HerbCore.is_domain_valid(constraint::ValidSMILES, grammar::ContextSensitiveGrammar)
    # TODO: Implement the logic to check if the domain of the ValidSMILES constraint is valid
    return true
end

function HerbCore.update_rule_indices!(
        constraint::ValidSMILES,
        n_rules::Integer,
        mapping::AbstractDict{<:Integer, <:Integer},
        constraints::Vector{<:AbstractConstraint}
)
    # TODO: Implement the logic to update the rule indices of the ValidSMILES constraint
    update_grammar_data!(constraint.grammar_data, mapping)
end

function HerbCore.update_rule_indices!(constraint::ValidSMILES, n_rules::Integer)
    # TODO: Implement the logic to update the rule indices of the ValidSMILES constraint
end

function HerbConstraints.on_new_node(
        solver::Solver, constraint::ValidSMILES, path::Vector{Int}
)
    if solver isa GenericSolver
        node = get_node_at_location(solver, path)
        type = get_node_type(solver.grammar, node)

        if type == :molecule
            constraint = GenericAtomConnections(path, constraint.grammar_data)
            HerbConstraints.post!(solver, constraint)
        end
    else
        node = get_node_at_location(solver, path)
        grammar = solver.grammar

        # Only proceed if node is filled and is a molecule rule
        if !isfilled(node) || grammar.types[HerbCore.get_rule(node)] != :molecule
            return nothing
        end

        # Post the atom constraints
        post_atom_constraints!(solver, path, constraint.grammar_data)

        # Post the ringbond constraints
        post_ringbond_constraints!(solver, path, constraint.grammar_data)
    end
end

function is_valid(candidate::Molecule, constraints::ValidSMILES)
    bond_orders = Dict(
        single => 1.0, aromatic => 1.5, double => 2.0, triple => 3.0, quadruple => 4.0)

    atoms = candidate.atoms
    bonds = candidate.bonds

    for (atom_index, atom) in enumerate(atoms)
        atom_str = atom.name
        if !haskey(constraints.atom_valences, atom_str)
            return false
        end
        valence = constraints.atom_valences[atom_str]
        connected_bonds = filter(b -> b.from == atom_index || b.to == atom_index, bonds)
        if isempty(connected_bonds)
            if valence == 0
                continue
            else
                return false
            end
        end
        total_bond_order = round(Int, sum(bond_orders[b.bond_type]
        for b in connected_bonds))
        if total_bond_order != valence
            return false
        end
    end

    return true
end
