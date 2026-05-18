@testitem "Bottom-Up Molecule Synthesizer" begin
    # Create some atoms
    atoms = [Atom("H"), Atom("O"), Atom("C")]

    settings = SynthesizerSettings(
        max_depth = 7, options = Dict{Symbol, Any}(:unique_candidates => true))
    unique_molecules = synthesize_molecules(atoms, settings)
    println("Number of unique molecules synthesized (RDKit canonicalization): ",
        length(unique_molecules))
    for mol in unique_molecules
        println("Unique molecule: ", to_SMILES(mol))
    end

    settings = SynthesizerSettings(
        max_depth = 7,
        iterator = BottomUp,    
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    bottom_up_molecules = synthesize_molecules(atoms, settings)
    println("Number of unique molecules synthesized with bottom-up iterator: ",
        length(bottom_up_molecules))
    for mol in bottom_up_molecules
        println("Bottom-up molecule: ", to_SMILES(mol))
    end
    @test length(bottom_up_molecules) > 0
    @test length(bottom_up_molecules) == length(unique_molecules)
end