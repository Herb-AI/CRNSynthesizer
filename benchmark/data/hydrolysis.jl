# Methyl Acetate Hydrolysis Reaction Experiment
using CRNSynthesizer, Catalyst, OrdinaryDiffEq


function hydrolysis_problem(;
    selected_known_indices=1:4,
    selected_expected_indices=1:4)
    # Define the Hydrolysis reaction network
    rn = @reaction_network begin
        p1, CH₃COOCH₃ + H₂O --> CH₃COOH + CH₃OH
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:CH₃COOCH₃ => 2.0, :H₂O => 4.0, :CH₃COOH => 0.0, :CH₃OH => 0.0]
    p = [:p1 => 0.2]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob, Tsit5())
    data_sol = solve(prob, Tsit5(), saveat=1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_CH3COOCH3 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃COOCH₃][1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]
    expected_CH3COOH = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃COOH][1:end]
    expected_CH3OH = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃OH][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[O]-[C](-[H])(-[H])-[H]"), # CH₃COOCH₃
        from_SMILES("[H]-[O]-[H]"),                                             # H₂O
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[O]-[H]"),                 # CH₃COOH
        from_SMILES("[C](-[H])(-[H])(-[H])-[O]-[H]")                            # CH₃OH
    ]
    
    # All expected profiles
    all_expected = [
        expected_CH3COOCH3,
        expected_H2O,
        expected_CH3COOH,
        expected_CH3OH
    ]

    # Define the known molecules based on selected indices
    known_molecules = all_molecules[selected_known_indices]

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
        known_molecules=known_molecules,
        expected_profiles=expected_profiles,
        time_data=time_data
    )

    return problem
end


function hydrolysis_network()
    # Define molecules using SMILES
    all_molecules = [
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[O]-[C](-[H])(-[H])-[H]"), # CH₃COOCH₃
        from_SMILES("[H]-[O]-[H]"),                                             # H₂O
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[O]-[H]"),                 # CH₃COOH
        from_SMILES("[C](-[H])(-[H])(-[H])-[O]-[H]")                            # CH₃OH
    ]

    # Define the reaction: CH₃COOCH₃ + H₂O --> CH₃COOH + CH₃OH
    reaction = CRNSynthesizer.Reaction(
        [(1, all_molecules[1]), (1, all_molecules[2])],
        [(1, all_molecules[3]), (1, all_molecules[4])]
    )

    return CRNSynthesizer.ReactionNetwork([reaction])
end
