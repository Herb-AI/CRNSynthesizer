using CRNSynthesizer, Catalyst

function ethylene_problem(; selected_known_indices = 1:3, selected_expected_indices = 1:3)
    # Define the ethylene glycol formation reaction network
    rn = @reaction_network begin
        p1, H₂O + C₂H₄O --> C₂H₆O₂
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:H₂O => 4.0, :C₂H₄O => 2.0, :C₂H₆O₂ => 0.0]
    p = [:p1 => 0.2]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob)
    data_sol = solve(prob; saveat = 1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]
    expected_C2H4O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₂H₄O][1:end]
    expected_C2H6O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₂H₆O₂][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[H]-[O]-[H]"),           # H₂O
        from_SMILES("[H]-[C](-[O]-1)(-[H])-[C]-1(-[H])-[H]"),       # Ethylene oxide (three-membered ring)
        from_SMILES("[H]-[O]-[C](-[H])(-[H])-[C](-[H])(-[H])-[O]-[H]")         # Ethylene glycol
    ]

    # All expected profiles
    all_expected = [expected_H2O, expected_C2H4O, expected_C2H6O2]

    # Define the known molecules based on selected indices
    known_molecules = all_molecules[selected_known_indices]
    goal_molecules = Vector{Molecule}()
    for i in eachindex(all_molecules)
        if !(i in selected_known_indices)
            push!(goal_molecules, all_molecules[i])
        end
    end

    atom_valences = get_valences_from_molecules(all_molecules)

    reaction = CRNSynthesizer.Reaction(
        nothing, [(1, all_molecules[1]), (1, all_molecules[2])], [(1, all_molecules[3])]
    )
    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Build expected profiles dictionary
    expected_profiles = Dict{Molecule, Vector{Float64}}()
    for (i, idx) in enumerate(selected_expected_indices)
        if idx in selected_known_indices
            # Find position of this molecule in known_molecules
            known_idx = findfirst(j -> j == idx, selected_known_indices)
            expected_profiles[known_molecules[known_idx]] = all_expected[idx]
        end
    end

    # Define the problem
    problem = ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network;
        expected_profiles = expected_profiles,
        time_data = time_data
    )

    return problem
end
