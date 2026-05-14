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
