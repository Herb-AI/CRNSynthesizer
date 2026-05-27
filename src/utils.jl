# Function to convert numbers to a subscript format
function convert_to_subscript(name::String)
    subscript_map = Dict(
        '0' => '₀',
        '1' => '₁',
        '2' => '₂',
        '3' => '₃',
        '4' => '₄',
        '5' => '₅',
        '6' => '₆',
        '7' => '₇',
        '8' => '₈',
        '9' => '₉'
    )
    result = ""
    for c in name
        if isdigit(c)
            result *= subscript_map[c]
        else
            result *= string(c)
        end
    end
    return result
end

function remove_all!(a::Vector{Vector{Int}}, item::Vector{Int})
    deleteat!(a, findall(x->x==item, a))
end

#########################################
############ Herb Functions #############
#########################################

# HerbCore.is_domain_valid(c, grammar) = true
# HerbCore.update_rule_indices!(c, new_indices) = nothing

get_rules(rn::RuleNode)::Vector{Int} = [rn.ind]
function get_rules(hole::AbstractHole)::Vector{Int}
    return findall(hole.domain)
end
function get_rules(hole::AbstractUniformHole)::Vector{Int}
    return findall(hole.domain)
end

# function get_node_type(grammar::AbstractGrammar, node)
#     # Get the type of the first node rule
#     rules = get_rules(node)
#     type = grammar.types[rules[1]]

#     # Check if it holds for all rules
#     @assert all(grammar.types[rules[i]] == type for i in eachindex(rules))

#     return type
# end

# TODO: check lazy search
function get_node_type(
        grammar::AbstractGrammar, node::Union{AbstractUniformHole, AbstractHole}
)
    return grammar.types[findfirst(node.domain)]
end
function get_node_type(grammar::AbstractGrammar, node::RuleNode)
    return grammar.types[node.ind]
end

function c_remove_all_but!(set::HerbConstraints.StateSparseSet, vals::Vector{Int})::Bool
    to_remove = filter(v -> v ∉ vals, collect(set))
    if isempty(to_remove)
        return false
    end
    for v in to_remove
        HerbConstraints.remove!(set, v)
    end
    return true
end

function c_remove_all_but!(
        solver::UniformSolver, path::Vector{Int}, rules::Vector{Int}, fix_point::Bool = true
)
    node = get_node_at_location(solver, path)
    if node isa RuleNode
        return nothing
    end

    hole = get_hole_at_location(solver, path)
    if c_remove_all_but!(hole.domain, rules)
        if isempty(hole.domain)
            HerbConstraints.set_infeasible!(solver)
        end
        HerbConstraints.notify_tree_manipulation(solver, path)
        if fix_point
            HerbConstraints.fix_point!(solver)
        end
    end
end

function c_remove_all_but!(
        solver::GenericSolver,
        path::Vector{Int},
        rules_to_keep::Vector{Int},
        fix_point::Bool = true
)
    hole = get_hole_at_location(solver, path)

    bit_to_keep = BitVector(falses(length(hole.domain)))
    bit_to_keep[rules_to_keep] .= true

    updated_domain = hole.domain .& bit_to_keep
    if hole.domain != updated_domain
        hole.domain = updated_domain
        HerbConstraints.simplify_hole!(solver, path)
        HerbConstraints.notify_tree_manipulation(solver, path)
        if fix_point
            HerbConstraints.fix_point!(solver)
        end
    end
end

function c_remove!(
        solver::UniformSolver, path::Vector{Int}, rules::Vector{Int}, fix_point::Bool = true
)
    node = get_node_at_location(solver, path)
    if node isa RuleNode
        return nothing
    end

    #remove the rule_index from the state sparse set of the hole
    hole = get_hole_at_location(solver, path)
    domain_updated = false
    for rule_index in rules
        if HerbConstraints.remove!(hole.domain, rule_index)
            domain_updated = true
        end
    end
    if domain_updated
        if isempty(hole.domain)
            HerbConstraints.set_infeasible!(solver)
        end
        HerbConstraints.notify_tree_manipulation(solver, path)
        if fix_point
            HerbConstraints.fix_point!(solver)
        end
    end
end

# function remove_all_but!(solver::GenericSolver, path::Vector{Int}, rules_to_keep::Vector{Int})
#     hole = get_hole_at_location(solver, path)
#     domain_updated = false
#     for i in eachindex(hole.domain)
#         if hole.domain[i] && !(i in rules_to_keep)
#             domain_updated = true
#             hole.domain[i] = false
#         end
#     end
#     if domain_updated
#         simplify_hole!(solver, path)
#         notify_tree_manipulation(solver, path)
#         fix_point!(solver)
#     end
# end

function c_remove_above!(
        solver::UniformSolver, path::Vector{Int}, rule_index::Int; fix_point::Bool = true
)
    hole = get_hole_at_location(solver, path)
    if HerbConstraints.remove_above!(hole.domain, rule_index)
        if isempty(hole.domain)
            HerbConstraints.set_infeasible!(solver)
        end
        HerbConstraints.notify_tree_manipulation(solver, path)
        if fix_point
            HerbConstraints.fix_point!(solver)
        end
    end
end

# these are the same functions as in herbcore, however called custom to prevent percompile errors
macro c_rulenode(ex::Union{Integer, Expr})
    custom_shorthand2rulenode(ex)
end
function custom_shorthand2rulenode(ex::Int)
    return RuleNode(ex)
end
using MacroTools: @capture, postwalk, prewalk, isexpr
function custom_shorthand2rulenode(ex::Expr)
    # holes with children
    ex = postwalk(ex) do x
        @capture(x, holetype_[Bool[domain__]]{children__}) ||
            @capture(x, holetype_[domain__]{children__}) ||
            return x
        hole_constructor = _get_hole_name(holetype)
        return :($hole_constructor(BitVector([$(domain...)]), {$(children...)}))
    end

    # holes without children
    # need to prewalk here becuase otherwise we try to parse UniformHole[Bool[...]] first as a
    # hole with type Bool, instead of a UniformHole
    ex = prewalk(ex) do x
        @capture(x, holetype_Symbol[Bool[domain__]]) ||
            @capture(x, holetype_Symbol[domain__]) ||
            return x
        hole_constructor = _get_hole_name(holetype)
        return :($hole_constructor(BitVector([$(domain...)])))
    end

    # rulenodes with children
    ex = postwalk(ex) do x
        @capture(x, index_Int{children__}) || return x
        return :(RuleNode($index, {$(children...)}))
    end

    # rulenodes without children
    ex = postwalk(ex) do x
        @capture(x, {children__}) || return x
        children = [isexpr(child, Int, Integer) ? :(RuleNode($child)) : child
                    for child in children]
        return :([$(children...)])
    end

    # varnodes
    ex = postwalk(ex) do x
        @capture(x, name_Symbol) || return x
        if name == :RuleNode || name == :Hole || name == :BitVector || name == :UniformHole
            return :($x)
        else
            return :(VarNode($(QuoteNode(name))))
        end
    end

    return esc(ex)
end

function c_add_rule!(g::AbstractGrammar, e::Expr; recalculate::Bool = false)
    if e.head == :(=) && typeof(e.args[1]) == Symbol
        s = e.args[1]# Name of return type
        rule = e.args[2]# expression?
        rvec = Any[]
        parse_rule!(rvec, rule)
        for r in rvec
            # Only add a rule if it does not exist yet. Check for existance
            # with strict equality so that true and 1 are not considered
            # equal. that means we can't use `in` or `∈` for equality checking.
            if !any(
                s == type && (r === rule || typeof(r)==Expr && r == rule)
            for
            (type, rule) in zip(g.types, g.rules)
            )
                push!(g.rules, r)
                push!(g.iseval, iseval(rule))
                push!(g.types, s)
                g.bytype[s] = push!(get(g.bytype, s, Int[]), length(g.rules))
            end
        end
    else
        throw(ArgumentError("Invalid rule: $e. Rules must be of the form `Symbol = Expr`"))
    end
    alltypes = collect(keys(g.bytype))

    if recalculate
        # is_terminal and childtypes need to be recalculated from scratch, since a new type might 
        # be added that was used as a terminal symbol before.
        g.isterminal = [isterminal(rule, alltypes) for rule in g.rules]
        g.childtypes = [get_childtypes(rule, alltypes) for rule in g.rules]
        g.bychildtypes = [BitVector([g.childtypes[i1] == g.childtypes[i2]
                                     for i2 in 1:length(g.rules)])
                          for i1 in 1:length(g.rules)]
        g.domains = Dict(
            type => BitArray(r in g.bytype[type] for r in 1:length(g.rules))
        for
        type in keys(g.bytype)
        )

        # update grammar constraints to enforce domain correctness
        for c in g.constraints
            HerbCore.update_rule_indices!(c, g)
        end
    end
    return g
end

# # My constraints are not stateless, so I need a deep copy TODO: Check performance hit and if there are better ways
# function Base.copy(state::SolverState) 
#     tree = deepcopy(state.tree)
#     active_constraints = deepcopy(state.active_constraints)
#     SolverState(tree, active_constraints, state.isfeasible)
# end
