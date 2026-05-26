using CRNSynthesizer

function syn_molecules()::Vector{Molecule}
    compact_mols = ["CC1(C)OB(c2cn[nH]c2)OC1(C)C", "O=[N+]([O-])c1ccc2c(c1)c(Br)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1", "O=[N+]([O-])c1ccc2c(c1)c(-c1cn[nH]c1)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1", "CC1(C)OB(Br)OC1(C)C"]
    return map(from_SMILES ∘ make_smiles_custom_explicit ∘ make_smiles_rdkit_explicit, compact_mols)
end

function syn_problem(; selected_known_indices = [1,2,3])
    # Time series data is not generated for this problem, as the focus is to
    # just find the balanced reaction

    # All possible molecules
    all_molecules = syn_molecules()
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
        nothing, [(1, all_molecules[1]), (1, all_molecules[2])], [(1, all_molecules[3]), (1, all_molecules[4])]
    )

    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Define the problem
    problem = ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network
    )

    return problem
end
