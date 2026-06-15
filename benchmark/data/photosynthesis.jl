# Photosynthesis Reaction Experiment
using CRNSynthesizer, Catalyst, OrdinaryDiffEq


function photosynthesis_problem(;
    selected_known_indices=1:4,
    selected_expected_indices=1:4)
    # Define the Photosynthesis reaction network
    rn = @reaction_network begin
        p1, 6CO₂ + 6H₂O --> C₆H₁₂O₆ + 6O₂
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:CO₂ => 6.0, :H₂O => 6.0, :C₆H₁₂O₆ => 0.0, :O₂ => 0.0]
    p = [:p1 => 0.2]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob, Tsit5())
    data_sol = solve(prob, Tsit5(), saveat=1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_CO2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CO₂][1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]
    expected_C6H12O6 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₆H₁₂O₆][1:end]
    expected_O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:O₂][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[O]=[C]=[O]"), # CO₂
        from_SMILES("[H]-[O]-[H]"), # H₂O
        from_SMILES("[C](-[H])(-[H])(-[O]-[H])-[C](-1)(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[O]-1"),   # C₆H₁₂O₆
        from_SMILES("[O]=[O]")      # O₂
    ]
    
    # All expected profiles
    all_expected = [
        expected_CO2,
        expected_H2O,
        expected_C6H12O6,
        expected_O2
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


function photosynthesis_network()
    # Define molecules using SMILES
    all_molecules = [
        from_SMILES("[O]=[C]=[O]"), # CO₂
        from_SMILES("[H]-[O]-[H]"), # H₂O
        from_SMILES("[C](-[H])(-[H])(-[O]-[H])-[C](-1)(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[O]-1"),   # C₆H₁₂O₆
        from_SMILES("[O]=[O]")      # O₂
    ]

    # Define the reaction: 6 CO₂ + 6 H₂O --> C₆H₁₂O₆ + 6 O₂
    reaction = CRNSynthesizer.Reaction(
        [(6, all_molecules[1]), (6, all_molecules[2])],
        [(1, all_molecules[3]), (6, all_molecules[4])]
    )

    return CRNSynthesizer.ReactionNetwork([reaction])
end
