@testitem "Molecule Synthesizer" begin
    using DataStructures
    # Create some atom valences
    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2, "[C]" => 4, "[N]" => 3)

    settings = SynthesizerSettings(
        max_depth = 4, options = Dict{Symbol, Any}(:unique_candidates => true))
    unique_molecules = synthesize_molecules(atom_valences; settings = settings)
    @test length(unique_molecules) > 0
    #println("Number of unique molecules synthesized (RDKit canonicalization): ",
    # length(unique_molecules))

    target_molecule = "CC(C)(C)OC(=O)CONC(=O)NCc1cccc2ccccc12"
    # println("Target molecule: ", target_molecule)
    entry_fragments, starting_fragments = parse_molecule_to_fragment_rules(target_molecule)
    #println("Fragment rules: ", entry_fragments)
    #println("Starting fragments: ", starting_fragments)

    
    settings = SynthesizerSettings(
        max_depth = 4, options = Dict{Symbol, Any}(:unique_candidates => true))
    fragment_molecules = synthesize_molecules(
        atom_valences; settings = settings, fragment_rules = entry_fragments,
        starting_fragments = starting_fragments)
    @test length(fragment_molecules) > 0
    @test length(fragment_molecules) > length(unique_molecules)
    #println(fragment_molecules)
    #println("Number of unique fragment-based molecules synthesized: ",
      #length(fragment_molecules))
end

@testitem "SMILES Unique Digit Remapping" begin
    # Test a more complex overlapping ring system to make sure they remain correct
    input_smiles2 = "C1CC2CC1CC2"
    explicit_smiles2 = CRNSynthesizer.make_smiles_custom_explicit(input_smiles2)
    @test occursin("1", explicit_smiles2)
    @test occursin("2", explicit_smiles2)
end

