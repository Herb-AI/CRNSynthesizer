@testitem "Molecule Synthesizer" begin
    using HerbCore
    using HerbGrammar

    # Create some atoms
    atoms = [Atom("H"), Atom("O"), Atom("C")]

    # Synthesize all molecules from root symbol
    settings = SynthesizerSettings(max_depth = 7)
    molecules = synthesize_molecules(atoms, settings)
    @test length(molecules) > 0

    # Test synthesizing from a partial tree
    grammar = SMILES_grammar(atoms; settings = settings)

    # Structure for [H]
    hydrogen_structure = RuleNode(4, [ # structure
        RuleNode(16), # [H]
        RuleNode(9), # empty ringbonds
        RuleNode(6) # empty branches
    ])

    # Structure for [O]
    oxygen_structure = RuleNode(4, [ # structure
        RuleNode(17), # [O]
        RuleNode(9), # empty ringbonds
        RuleNode(6) # empty branches
    ])

    # Inner chain: [O] - [Hole]
    oxygen_tail = RuleNode(3, [ # combine chain
        RuleNode(13), # -
        oxygen_structure,
        Hole(get_domain(grammar, :chain)) # hole for the rest of the chain
    ])

    # Outer chain: [H] - ([O] (- [Hole]))
    partial_tree = RuleNode(1, [ # molecule
        RuleNode(3, [ # combine chain
            RuleNode(13), # -
            hydrogen_structure,
            oxygen_tail
        ])
    ])
    
    molecules_with_substructure = synthesize_molecules(atoms, settings, partial_tree)
    @test length(molecules_with_substructure) > 0

    #All synthesized molecules should contain atoms H and O
    @test all(molecule -> Atom("H") in molecule.atoms && Atom("O") in molecule.atoms, molecules_with_substructure)
    # There should be fewer molecules that contain the substructure than total molecules
    @test length(molecules_with_substructure) < length(molecules)
end
