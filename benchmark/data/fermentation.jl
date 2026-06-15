# Sugar Fermentation Reaction Experiment
using CRNSynthesizer, Catalyst, OrdinaryDiffEq


function fermentation_problem(;
    selected_known_indices=1:8,
    selected_expected_indices=1:8)
    # Define the Fermentation reaction network
    rn = @reaction_network begin
        p1, C₆H₁₂O₆ + 2NAD --> 2C₃H₄O₃ + 2NADHH         # Glycolysis step: glucose to pyruvate
        p2, C₃H₄O₃ + NADHH --> CH₃CHOHCOOH + NAD        # Lactic acid fermentation
        p3, C₃H₄O₃ --> CH₃CHO + CO₂                     # Acetaldehyde formation
        p4, NADHH + CH₃CHO --> CH₃CH₂OH + NAD           # Ethanol and carbon dioxide formation
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:C₆H₁₂O₆ => 1.0, :NAD => 6.0, :C₃H₄O₃ => 2.0, :NADHH => 2.0, :CH₃CHOHCOOH => 0.0, :CH₃CHO => 0.0, :CO₂ => 0.0, :CH₃CH₂OH => 0.0]
    p = [:p1 => 0.2, :p2 => 0.15, :p3 => 0.1, :p4 => 0.1]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob, Tsit5())
    data_sol = solve(prob, Tsit5(), saveat=1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_C6H12O6 = data_sol[:C₆H₁₂O₆][1:end]
    # expected_C6H12O6 = (0.9 .+ 0.5 * rand(length(time_data))) .* data_sol[:C₆H₁₂O₆][1:end]
    expected_NAD = (0.9 .+ 0.5 * rand(length(time_data))) .* data_sol[:NAD][1:end]
    expected_C3H4O3 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:C₃H₄O₃][1:end]
    expected_NADHH = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:NADHH][1:end]
    expected_CH3CHOHCOOH = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃CHOHCOOH][1:end]
    expected_CH3CHO = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃CHO][1:end]
    expected_CO2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CO₂][1:end]
    expected_CH3CH2OH = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₃CH₂OH][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[C](-[H])(-[H])(-[O]-[H])-[C](-1)(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[O]-1"),       # C6H12O6
        from_SMILES("[NAD]=[NAD]"),                                                             # TODO: use proper formulation, currently this is using a workaround representation to still simulate the hydrogen transfer
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[C](=[O])-[O]-[H]"),                       # C₃H₄O₃
        from_SMILES("[H]-[NAD]-[NAD]-[H]"),                                                     # TODO: use proper formulation, currently this is using a workaround representation to still simulate the hydrogen transfer
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])(-[O]-[H])-[C](=[O])-[O]-[H]"),             # CH₃CHOHCOOH
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])=[O]"),                                     # CH₃CHO
        from_SMILES("[O]=[C]=[O]"),                                                             # CO₂
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])(-[H])-[O]-[H]"),                           # CH₃CH₂OH
    ]
    
    # All expected profiles
    all_expected = [
        expected_C6H12O6,
        expected_NAD,
        expected_C3H4O3,
        expected_NADHH,
        expected_CH3CHOHCOOH,
        expected_CH3CHO,
        expected_CO2,
        expected_CH3CH2OH,
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
        time_data=time_data,
    )

    return problem
end


function fermentation_network()
    # Define molecules using SMILES
    all_molecules = [
        from_SMILES("[C](-[H])(-[H])(-[O]-[H])-[C](-1)(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[C](-[O]-[H])(-[H])-[O]-1"),       # C6H12O6
        from_SMILES("[NAD]=[NAD]"),                                                             # TODO: use proper formulation, currently this is using a workaround representation to still simulate the hydrogen transfer
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](=[O])-[C](=[O])-[O]-[H]"),                       # C₃H₄O₃
        from_SMILES("[H]-[NAD]-[NAD]-[H]"),                                                     # TODO: use proper formulation, currently this is using a workaround representation to still simulate the hydrogen transfer
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])(-[O]-[H])-[C](=[O])-[O]-[H]"),             # CH₃CHOHCOOH
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])=[O]"),                                     # CH₃CHO
        from_SMILES("[O]=[C]=[O]"),                                                             # CO₂
        from_SMILES("[C](-[H])(-[H])(-[H])-[C](-[H])(-[H])-[O]-[H]"),                           # CH₃CH₂OH
    ]

    # Define the reaction: C₆H₁₂O₆ + 2NAD --> 2C₃H₄O₃ + 2NADHH         # Glycolysis step: glucose to pyruvate
    reaction1 = CRNSynthesizer.Reaction(
        [(1, all_molecules[1]), (2, all_molecules[2])],
        [(2, all_molecules[3]), (2, all_molecules[4])]
        )
    # Define the reaction: C₃H₄O₃ + NADHH --> CH₃CHOHCOOH + NAD       # Lactic acid fermentation
    reaction2 = CRNSynthesizer.Reaction(
        [(1, all_molecules[3]), (1, all_molecules[4])],
        [(1, all_molecules[5]), (1, all_molecules[2])]
        )
    # Define the reaction: C₃H₄O₃ --> CH₃CHO + CO₂                     # Acetaldehyde formation
    reaction3 = CRNSynthesizer.Reaction(
        [(1, all_molecules[3])],
        [(1, all_molecules[6]), (1, all_molecules[7])]
        )
    # Define the reaction: NADHH + CH₃CHO --> CH₃CH₂OH + NAD           # Ethanol and carbon dioxide formation
    reaction4 = CRNSynthesizer.Reaction(
        [(1, all_molecules[4]), (1, all_molecules[6])],
        [(1, all_molecules[8]), (1, all_molecules[2])]
        )

    return CRNSynthesizer.ReactionNetwork([reaction1, reaction2, reaction3, reaction4])
end
