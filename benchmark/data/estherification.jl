# Estherification Reaction Experiment
using CRNSynthesizer, Catalyst

function estherification_problem(;
        selected_known_indices = 1:6, selected_expected_indices = 1:6
)
    # Define the reaction network
    rn = @reaction_network begin
        p1, CH₂O₂ + C₂H₆O --> C₃H₆O₂ + H₂O
        p2, C₃H₆O₂ + CH₄O --> C₂H₄O₂ + C₂H₆O
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [
        :CH₂O₂ => 2.0,
        :C₂H₆O => 1.0,
        :C₃H₆O₂ => 0.0,
        :H₂O => 0.0,
        :CH₄O => 2.0,
        :C₂H₄O₂ => 0.0
    ]
    p = [:p1 => 0.2, :p2 => 0.1]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob)
    data_sol = solve(prob; saveat = 1.0)

    # Gather the time data and expected values for the known species
    time_data = data_sol.t[1:end]
    expected_CH2O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₂O₂][1:end]
    expected_C2H6O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₂H₆O][1:end]
    expected_C3H6O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₃H₆O₂][1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]
    expected_CH4O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₄O][1:end]
    expected_C2H4O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₂H₄O₂][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[H]-[C](=[O])-[O]-[H]"),
        from_SMILES("[C](-[C](-[H])(-[H])-[O]-[H])(-[H])(-[H])-[H]"),
        from_SMILES("[H]-[C](=[O])-[O]-[C](-[H])(-[H])-[C](-[H])(-[H])-[H]"),
        from_SMILES("[H]-[O]-[H]"),
        from_SMILES("[C](-[H])(-[H])(-[H])-[O]-[H]"),
        from_SMILES("[C](-[H])(=[O])-[O]-[C](-[H])(-[H])-[H]")
    ]

    # All expected profiles
    all_expected = [
        expected_CH2O2,
        expected_C2H6O,
        expected_C3H6O2,
        expected_H2O,
        expected_CH4O,
        expected_C2H4O2
    ]

    # Define the known molecules based on selected indices
    known_molecules = all_molecules[selected_known_indices]
    goal_molecules = Vector{Molecule}()
    for i in eachindex(all_molecules)
        if !(i in selected_known_indices)
            push!(goal_molecules, all_molecules[i])
        end
    end

    atom_valences = get_valences_from_molecules(all_molecules)

    reaction1 = CRNSynthesizer.Reaction(
        nothing,
        [(1, all_molecules[1]), (1, all_molecules[2])],
        [(1, all_molecules[3]), (1, all_molecules[4])]
    )
    reaction2 = CRNSynthesizer.Reaction(
        nothing,
        [(1, all_molecules[3]), (1, all_molecules[5])],
        [(1, all_molecules[6]), (1, all_molecules[2])]
    )

    goal_network = CRNSynthesizer.ReactionNetwork([reaction1, reaction2])

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
