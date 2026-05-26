# Estherification Reaction Experiment
using CRNSynthesizer, Catalyst

function water_problem(; selected_known_indices = 1:3, selected_expected_indices = 1:3)
    # Define the water formation reaction network
    rn = @reaction_network begin
        p1, 2H₂ + O₂ --> 2H₂O
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:H₂ => 4.0, :O₂ => 2.0, :H₂O => 0.0]
    p = [:p1 => 0.2]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob)
    data_sol = solve(prob; saveat = 1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_H2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂][1:end]
    expected_O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:O₂][1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[H]-[H]"), from_SMILES("[O]=[O]"), from_SMILES("[H]-[O]-[H]")
    ]

    # All expected profiles
    all_expected = [expected_H2, expected_O2, expected_H2O]

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
        nothing, [(2, all_molecules[1]), (1, all_molecules[2])], [(2, all_molecules[3])]
    )
    

    goal_network =  CRNSynthesizer.ReactionNetwork([reaction])

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
