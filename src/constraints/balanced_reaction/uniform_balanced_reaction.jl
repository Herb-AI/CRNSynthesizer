struct LocalUniformBalancedReaction <: AbstractLocalConstraint
    path::Vector{Int}
    input_atom_paths::Vector{Vector{Int}}
    output_atom_paths::Vector{Vector{Int}}
end

import Base.==
function ==(a::LocalUniformBalancedReaction, b::LocalUniformBalancedReaction)
    return a.path == b.path && a.atom_paths == b.atom_paths
end

function Base.hash(a::LocalUniformBalancedReaction, h::UInt)
    return hash(a.output_atom_paths, hash(a.input_atom_paths, hash(a.path, h)))
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::LocalUniformBalancedReaction, path::Vector{Int}
)
    if path in constraint.input_atom_paths || path in constraint.output_atom_paths
        return true
    end

    return false
end

function HerbConstraints.propagate!(
        solver::Solver, constraint::LocalUniformBalancedReaction
)
    input_atoms = Dict{Int, Int}()
    possible_inputs = Dict{Int, Int}()
    unfilled_inputs = 0
    for input_path in constraint.input_atom_paths
        atom = get_node_at_location(solver, input_path)
        if !isfilled(atom)
            unfilled_inputs += 1
            rules = get_rules(atom)
            for rule in rules
                if haskey(possible_inputs, rule)
                    possible_inputs[rule] += 1
                else
                    possible_inputs[rule] = 1
                end
            end
        else
            rule = HerbCore.get_rule(atom)
            if haskey(input_atoms, rule)
                input_atoms[rule] += 1
            else
                input_atoms[rule] = 1
            end
        end
    end

    output_atoms = Dict{Int, Int}()
    possible_outputs = Dict{Int, Int}()
    unfilled_outputs = 0
    for output_path in constraint.output_atom_paths
        if !isfilled(get_node_at_location(solver, output_path))
            unfilled_outputs += 1
            rules = get_rules(get_node_at_location(solver, output_path))
            for rule in rules
                if haskey(possible_outputs, rule)
                    possible_outputs[rule] += 1
                else
                    possible_outputs[rule] = 1
                end
            end
        else
            rule = HerbCore.get_rule(get_node_at_location(solver, output_path))
            if haskey(output_atoms, rule)
                output_atoms[rule] += 1
            else
                output_atoms[rule] = 1
            end
        end
    end

    # println("Input atoms: ", input_atoms)
    # println("Output atoms: ", output_atoms)
    # println("Possible input atoms: ", possible_inputs)
    # println("Possible output atoms: ", possible_outputs)
    # println("Unfilled inputs: ", unfilled_inputs)
    # println("Unfilled outputs: ", unfilled_outputs)
    # println()

    # TODO: calculate the required atoms to be able to restrict more

    # Check if the input and output atoms are balanced
    if unfilled_inputs == 0 && unfilled_outputs == 0
        for (atom, count) in input_atoms
            if !haskey(output_atoms, atom) || output_atoms[atom] != count
                HerbConstraints.set_infeasible!(solver)
                return nothing
            end
        end
    end
end
