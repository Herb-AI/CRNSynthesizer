function catalyst_simulation(
        network::ReactionNetwork, problem::ProblemDefinition; verbose = false, max_iters = 10000
)
    parameters = []

    crn = []
    t = default_t()
    for reaction in network.reactions
        inputs, input_amounts = create_species(reaction.inputs, t)
        outputs, output_amounts = create_species(reaction.outputs, t)

        if reaction.rate === nothing
            sym = Symbol("k$(length(parameters) + 1)")
            parameter = @parameters $(sym)
            push!(parameters, parameter[1])
            push!(
                crn,
                Catalyst.Reaction(
                    parameter[1], inputs, outputs, input_amounts, output_amounts
                )
            )
            continue
        end

        push!(
            crn,
            Catalyst.Reaction(
                reaction.rate, inputs, outputs, input_amounts, output_amounts
            )
        )
    end

    @named crn = ReactionSystem(crn, t)
    crn = complete(crn)

    # Convert ReactionSystem to ODESystem and then to ODEFunction
    # TODO: Investigate the Specialization of the ODEFunction
    sys = ode_model(crn)
    sys = complete(sys)
    f = ODEFunction{true, SciMLBase.NoSpecialize}(sys)

    u0_dummy = fill(1.0, length(species(crn)))
    ps_dummy = fill(0.1, length(parameters))

    species_names, cost_function = catalyst_cost_function(problem)

    idxs = [findfirst(
                s -> string(s) == "var\"" * string(species_name) * "\"(t)", species(crn))
            for species_name in species_names]

    tspan = (problem.time_data[1], problem.time_data[end])
    oprob = ODEProblem{true, SciMLBase.NoSpecialize}(f, u0_dummy, tspan, ps_dummy)
    prob_generator(prob,
        p) = remake(
        prob; u0 = p[1:length(u0_dummy)], p = p[(length(u0_dummy) + 1):end]
    )

    loss_function = build_loss_objective(
        oprob,
        Tsit5(),
        cost_function,
        Optimization.AutoForwardDiff();
        maxiters = max_iters,
        verbose = false,
        save_idxs = idxs,
        prob_generator = prob_generator
    )

    initial_guesses = [u0_dummy; ps_dummy]

    optprob = OptimizationProblem(
        loss_function,
        initial_guesses;
        lb = [fill(0.0, length(u0_dummy)); fill(0.05, length(ps_dummy))],
        ub = [fill(10.0, length(u0_dummy)); ones(length(ps_dummy))]
    )
    optsol = solve(optprob, NLopt.LN_NELDERMEAD(); maxiters = max_iters)

    fitted_u0 = optsol.u[1:length(u0_dummy)]
    fitted_params = optsol.u[(length(u0_dummy) + 1):end]

    oprob_fitted = remake(oprob; u0 = fitted_u0, p = fitted_params)
    sol = solve(oprob_fitted, Tsit5())

    if verbose == true
        display(crn)
        display(reactions(crn))
        display(fitted_params)
        display(fitted_u0)
    end

    return optsol, sol
end

function create_species(species_tuples::Vector{Tuple{Int, Molecule}}, t::Symbolics.Num)
    results = []
    amounts = []
    for (count, molecule) in species_tuples
        # Use a compact representation of the molecule as the species name
        name = to_SMILES(molecule)
        s = @species $(Symbol(name))(t)
        push!(results, s[1])
        push!(amounts, count)
    end
    return results, amounts
end

function catalyst_cost_function(problem::ProblemDefinition)
    species_names = [to_SMILES(mol) for mol in collect(keys(problem.expected_profiles))]
    combined_data = nothing
    for value in values(problem.expected_profiles)
        if isnothing(combined_data)
            combined_data = value
        else
            combined_data = hcat(combined_data, value)
        end
    end
    return (species_names, L2Loss(problem.time_data, Array(combined_data')))
end

function score(network, problem)
    num_species = count_species(network)
    num_reactions = count_reactions(network)

    run = catalyst_simulation(network, problem)
    objective = run[1].objective

    bonds_changed_score = 0
    for reaction in network.reactions
        bonds_changed_score += bonds_changed(reaction)
    end

    return objective * (num_species + num_reactions + bonds_changed_score)
end
