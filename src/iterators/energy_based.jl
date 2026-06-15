abstract type AbstractEnergyIterator <: HerbSearch.TopDownIterator end

function HerbSearch.priority_function(
        ::AbstractEnergyIterator,
        ::AbstractGrammar,
        node::HerbCore.AbstractRuleNode,
        parent_value::Union{Real, Tuple{Vararg{Real}}},
        isrequeued::Bool
)

    # Use BFS order as the default priority function
    if isrequeued
        return parent_value
    else
        return parent_value + 1
    end
end

function HerbSearch.derivation_heuristic(::AbstractEnergyIterator, indices::Vector{Int})
    # Use heuristic that prioritizes grammer rules in order of declaration (prioritizing earlier rules)
    return sort(indices)
end

function HerbSearch.hole_heuristic(::AbstractEnergyIterator, node::AbstractRuleNode,
        max_depth::Int)::Union{HerbSearch.ExpandFailureReason, HerbSearch.HoleReference}
    # Use heuristic that prioritizes the leftmost hole (i.e. the hole that is first encountered in a left-to-right traversal of the tree)
    return HerbSearch.heuristic_leftmost(node, max_depth)
end





# ==============================================================================================================================================================
#  __  __              ____                  _   _____ _                 _             
# |  \/  |            |  _ \                | | |_   _| |               | |            
# | \  / | __ ___  __ | |_) | ___  _ __   __| |   | | | |_ ___ _ __ __ _| |_ ___  _ __ 
# | |\/| |/ _` \ \/ / |  _ < / _ \| '_ \ / _` |   | | | __/ _ \ '__/ _` | __/ _ \| '__|
# | |  | | (_| |>  <  | |_) | (_) | | | | (_| |  _| |_| ||  __/ | | (_| | || (_) | |   
# |_|  |_|\__,_/_/\_\ |____/ \___/|_| |_|\__,_| |_____|\__\___|_|  \__,_|\__\___/|_|   
# ==============================================================================================================================================================

@programiterator MaxBondIterator() <: AbstractEnergyIterator



Bond_order = Dict(single => 1, double => 2, triple => 3)
function get_max_bond_order(rulenode::RuleNode, grammar::AbstractGrammar)::Int
    rule_ind = rulenode.ind
    if grammar.types[rule_ind] == :reaction && HerbGrammar.isterminal(grammar, rulenode)
            reaction::Reaction = grammar.rules[rule_ind]
            return maximum((maximum((Bond_order[bond.bond_type] for bond in molecule.bonds), init=0) for (_, molecule) in reaction.inputs), init = 0)
    end
    if grammar.types[rule_ind] == :reaction
        # Only consider the bonds in the reactants to avoid counting bonds that are formed in the products
        return get_max_bond_order(rulenode.children[1], grammar)
    end
    if grammar.types[get_rule(rulenode)] == :molecule && HerbGrammar.isterminal(grammar, rulenode)
        molecule::Molecule = grammar.rules[rulenode.ind]
        return maximum((Bond_order[bond.bond_type] for bond in molecule.bonds), init = 0)
    end
    if grammar.types[rule_ind] == :bond
        # Given that for every grammar rule the smaller bond orders are defined before the larger bond orders
        # we can use the index of the bond rule as smaller orders will have lower indices.
        return rulenode.ind
    end
    return maximum((get_max_bond_order(c, grammar) for c in rulenode.children), init = 0)
end

function get_max_bond_order(hole::AbstractHole, grammar::AbstractGrammar)::Int
    if HerbCore.isfilled(hole)
        rule_ind = findfirst(hole.domain)
        if grammar.types[rule_ind] == :reaction && HerbGrammar.isterminal(grammar, hole)
            reaction::Reaction = grammar.rules[rule_ind]
            return maximum((maximum((Bond_order[bond.bond_type] for bond in molecule.bonds), init=0) for (_, molecule) in reaction.inputs), init = 0)
        end
        if grammar.types[rule_ind] == :reaction
            # Only consider the bonds in the reactants to avoid counting bonds that are formed in the products
            return get_max_bond_order(hole.children[1], grammar)
        end
        if grammar.types[get_rule(hole)] == :molecule && HerbGrammar.isterminal(grammar, hole)
            molecule::Molecule = grammar.rules[rule_ind]
            return maximum((Bond_order[bond.bond_type] for bond in molecule.bonds), init = 0)
        end
        if grammar.types[rule_ind] == :bond
            # Given that for every grammar rule the smaller bond orders are defined before the larger bond orders
            # we can use the index of the bond rule as smaller orders will have lower indices. 
            return rule_ind
        end
        return maximum((get_max_bond_order(c, grammar) for c in hole.children), init = 0)
    else
        return 0
    end
end

function HerbSearch.priority_function(
        iterator::MaxBondIterator,
        grammar::AbstractGrammar,
        node::HerbCore.AbstractRuleNode,
        parent_value::Union{Real, Tuple{Vararg{Real}}},
        isrequeued::Bool
)    
    # Use heuristic that prioritizes grammer rules with the lowest max bond order
    # (i.e. prioritizing reactions that break single bonds over reactions that break double or triple bonds)
    # Get the max bond count for each applicable grammar rule
    max_bond = get_max_bond_order(node, grammar)
    
    if max_bond == 0
        if isrequeued
            return parent_value
        else
            return parent_value + 1000 # Ensure that any node with a bond is prioritized over nodes without bonds
        end
    end
    return max_bond # Prioritizing nodes with fewer holes

end






# ==============================================================================================================================================================
#   _____       _ _          ______                              _____ _                 _             
#  |  __ \     | | |        |  ____|                            |_   _| |               | |            
#  | |  | | ___| | |_ __ _  | |__   _ __   ___ _ __ __ _ _   _    | | | |_ ___ _ __ __ _| |_ ___  _ __ 
#  | |  | |/ _ \ | __/ _` | |  __| | '_ \ / _ \ '__/ _` | | | |   | | | __/ _ \ '__/ _` | __/ _ \| '__|
#  | |__| |  __/ | || (_| | | |____| | | |  __/ | | (_| | |_| |  _| |_| ||  __/ | | (_| | || (_) | |   
#  |_____/ \___|_|\__\__,_| |______|_| |_|\___|_|  \__, |\__, | |_____|\__\___|_|  \__,_|\__\___/|_|   
#                                                   __/ | __/ |                                        
#                                                  |___/ |___/                                         
# ==============================================================================================================================================================

@programiterator DeltaEnergyIterator() <: AbstractEnergyIterator



function get_bond_energy(bond_order::Int)::Int
    if bond_order == 1
        return 361 # The energy of a single bond is approximated to 361 kJ/mol
    elseif bond_order == 2
        return 617 # The energy of a double bond is approximated to 617 kJ/mol
    elseif bond_order == 3
        return 838 # The energy of a triple bond is approximated to 838 kJ/mol
    else
        return 0
    end
end

function get_bond_energy(molecule::Molecule)::Int
    sum = 0
    for bond::Bond in molecule.bonds
        if bond.bond_type == single
            sum += get_bond_energy(1) 
        elseif bond.bond_type == double
            sum += get_bond_energy(2)
        elseif bond.bond_type == triple
            sum += get_bond_energy(3)
        end
    end
    return sum
end

function get_bond_energy(bond::String)::Int
    if bond == "-"
        return get_bond_energy(1)
    elseif bond == "="
        return get_bond_energy(2)
    elseif bond == "≡"
        return get_bond_energy(3)
    else
        return 0
    end
end

function get_bond_energy(reaction::Reaction)::Int
    # Subtract the bond energy of the products from the bond energy of the reactants to get the change in energy for the reaction.
    sum = 0
    for (amount_of_product, product::Molecule) in reaction.inputs
        sum += amount_of_product * get_bond_energy(product);
    end
    for (amount_of_reactant, reactant::Molecule) in reaction.outputs
        sum -= amount_of_reactant * get_bond_energy(reactant);
    end
    return sum
end



function get_delta_E(ruleType, ruleIndx, ruleChildren, grammar::AbstractGrammar, isProduct::Bool = false)::Int
    if ruleType == :reaction && 
        HerbGrammar.isterminal(grammar, ruleIndx)
        reaction = grammar.rules[ruleIndx]
        return get_bond_energy(reaction)
    end
    if ruleType == :reaction
        # Subtract the bond energy of the products from the bond energy of the reactants to get the change in energy for the reaction.
        return get_delta_E(ruleChildren[1], grammar, false) -
               get_delta_E(ruleChildren[2], grammar, true)
    end
    if ruleType == :molecule &&
       HerbGrammar.isterminal(grammar, ruleIndx)
        molecule = grammar.rules[ruleIndx]
        return get_bond_energy(molecule)
    end
    if ruleType == :bond
        bond = grammar.rules[ruleIndx]
        return get_bond_energy(bond)
    end
    return sum((get_delta_E(c, grammar, isProduct) for c in ruleChildren), init = 0)
end


function get_delta_E(
        rulenode::RuleNode, grammar::AbstractGrammar, isProduct::Bool = false)::Int
    ruleType = grammar.types[get_rule(rulenode)]
    ruleIndx = rulenode.ind
    return get_delta_E(ruleType, ruleIndx, rulenode.children, grammar, isProduct)
end

function get_delta_E(
        hole::AbstractHole, grammar::AbstractGrammar, isProduct::Bool = false)::Int
    if HerbCore.isfilled(hole)
        ruleType = grammar.types[get_rule(hole)]
        ruleIndx = findfirst(hole.domain)
        return get_delta_E(ruleType, ruleIndx, hole.children, grammar, isProduct)
    else
        # As large reactions that break many bonds are less likely to occur, we can use a heuristic that prioritizes reactions with less holes.
        # This is done by returning a positive priority for holes in the reactants and a negative priority for holes in the products, 
        # to encourage the search to fill holes and give a rough estimate of the expected energy change for a reaction with a hole.
        return isProduct ? -100 : 100 
    end
end


function HerbSearch.priority_function(
        iterator::DeltaEnergyIterator,
        grammar::AbstractGrammar,
        node::HerbCore.AbstractRuleNode,
        parent_value::Union{Real, Tuple{Vararg{Real}}},
        isrequeued::Bool
)
    # Use heuristic that prioritizes grammer rules with the lowest change in energy (i.e. prioritizing reactions that are more 'exothermic')
    delta_energy = get_delta_E(node, grammar)
    return delta_energy
end