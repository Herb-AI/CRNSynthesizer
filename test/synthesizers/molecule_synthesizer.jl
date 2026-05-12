@testitem "Molecule Synthesizer" begin
    # Create some atoms
    atoms = [Atom("H"), Atom("O"), Atom("C")]

    # Synthesize all molecules from root symbol
    settings = SynthesizerSettings(max_depth = 8, options = Dict{Symbol, Any}(:unique_candidates => true))
    molecules = synthesize_molecules(atoms, settings)
    @test length(molecules) > 0
    println("Number of molecules synthesized (original unique check): ", length(molecules))

    # Use RDKit SMILES canonicalization to ensure uniqueness of molecules
    settings = SynthesizerSettings(max_depth = 8, options = Dict{Symbol, Any}(:new_unique_candidates => true))
    unique_molecules = synthesize_molecules(atoms, settings)
    @test length(unique_molecules) > 0
    @test length(unique_molecules) <= length(molecules)
    println("Number of unique molecules synthesized (RDKit canonicalization): ", length(unique_molecules))

    target_molecule = "C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N"
    # println("Target molecule: ", target_molecule)
    fragment_rules = parse_molecule_to_fragment_rules(target_molecule; min_digit = 3)
    # println("Fragment rules: ", fragment_rules)

    settings = SynthesizerSettings(max_depth = 8, options = Dict{Symbol, Any}(:new_unique_candidates => true))
    fragment_molecules = synthesize_molecules(atoms, settings; fragment_rules = fragment_rules)
    @test length(fragment_molecules) > 0
    @test length(fragment_molecules) >= length(unique_molecules)
    println("Number of unique fragment-based molecules synthesized: ", length(fragment_molecules))
end
