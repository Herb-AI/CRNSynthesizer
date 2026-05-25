using CRNSynthesizer

function syn_molecules()::Vector{Molecule}
    compact_mols = ["CC(C)(C)OC(=O)CONC(=O)NCc1cccc2ccccc12", "O", "O=C(O)CONC(=O)NCc1cccc2ccccc12", "CC(C)(C)O"]
    return map(from_SMILES ∘ make_smiles_custom_explicit ∘ make_smiles_rdkit_explicit, compact_mols)
end

function syn_problem(; selected_known_indices = [1,3])
    # Time series data is not generated for this problem, as the focus is to
    # just find the balanced reaction

    # All possible molecules
    all_molecules = syn_molecules()
    # Define the known molecules based on selected indices
    known_molecules = all_molecules[selected_known_indices]

    # Define the problem
    problem = ProblemDefinition(;known_molecules = known_molecules)

    return problem
end

function syn_network()
    # Define molecules using SMILES
    all_molecules = syn_molecules()

    # Define the reaction
    reaction = CRNSynthesizer.Reaction(
        nothing, [(1, all_molecules[1]), (1, all_molecules[2])], [(1, all_molecules[3]), (1, all_molecules[4])]
    )

    return CRNSynthesizer.ReactionNetwork([reaction])
end
