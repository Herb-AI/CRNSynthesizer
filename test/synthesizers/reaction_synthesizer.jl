#= The synthesis from atoms is not supported with fragment rules
Based on Wijers conclusions that pipelines without molecule step
are not feasible, the synthesis of networks directly from atoms 
is out of scope for a top-down iterator.
 @testitem "Reaction Synthesizer (from Atoms)" begin

    # Create some atoms
    a1 = Atom("H")
    a2 = Atom("O")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 8)
    reactions = synthesize_reactions([a1, a2], settings)
    @test length(reactions) > 0

    # TODO: Add more specific tests for the synthesized reactions
end =#

@testitem "Reaction Synthesizer (from Molecules)" begin

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 8)
    reactions = synthesize_reactions([m1, m2, m3], settings)
    @test length(reactions) > 0

end

@testitem "Reaction Synthesizer (atom -> molecules -> reactions pipeline)" begin
    using DataStructures

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2)
    problem = ProblemDefinition(
        atom_valences,
        [m1, m2, m3],
        Molecule[],
        ReactionNetwork(Reaction[])
    )

    molecule_settings = SynthesizerSettings(max_depth = 4, max_programs = 2)
    reaction_settings = SynthesizerSettings(max_depth = 8)

    reactions, molecules = synthesize_reactions(
        problem,
        molecule_settings,
        reaction_settings;
        initial_molecules_count = 3
    )

    @test length(molecules) >= 3
    @test length(reactions) > 0
end



