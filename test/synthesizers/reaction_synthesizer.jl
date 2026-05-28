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

@testitem "Reaction Synthesizer (with partial/inbalanced reactions)" begin
    using DataStructures
    using HerbConstraints
    using HerbCore

    # Molecules
    m1 = from_SMILES("[H]-[H]")       # H2
    m2 = from_SMILES("[O]=[O]")       # O2
    m3 = from_SMILES("[H]-[O]-[H]")   # H2O

    # Expected target reaction: 2 H2 + O2 -> 2 H2O
    expected = Reaction(nothing, [(2, m1), (1, m2)], [(2, m3)], false)

    # Synthesizer settings
    settings = SynthesizerSettings(max_depth = 8)

    # 1. H2O missing variation (2 H2 + O2 -> ?)
    # Prefilled inputs: 2 H2 + O2
    # Prefilled outputs: empty
    # Pool has all three. H2 and O2 will be excluded, so only H2O can fill outputs.
    partial_h2o_missing = Reaction(nothing, [(2, m1), (1, m2)], Tuple{Int, Molecule}[], false)
    reactions_h2o_missing = synthesize_reactions([m1, m2, m3], settings; partial_reaction = partial_h2o_missing)
    @test expected in reactions_h2o_missing

    # 2. O2 missing variation (2 H2 + ? -> 2 H2O)
    # Prefilled inputs: 2 H2
    # Prefilled outputs: 2 H2O
    # Pool has all three. H2 and H2O will be excluded, so only O2 can fill inputs.
    partial_o2_missing = Reaction(nothing, [(2, m1)], [(2, m3)], false)
    reactions_o2_missing = synthesize_reactions([m1, m2, m3], settings; partial_reaction = partial_o2_missing)
    @test expected in reactions_o2_missing

    # 3. H2 missing variation (? + O2 -> 2 H2O)
    # Prefilled inputs: O2
    # Prefilled outputs: 2 H2O
    # Pool has all three. O2 and H2O will be excluded, so only H2 can fill inputs.
    partial_h2_missing = Reaction(nothing, [(1, m2)], [(2, m3)], false)
    reactions_h2_missing = synthesize_reactions([m1, m2, m3], settings; partial_reaction = partial_h2_missing)
    @test expected in reactions_h2_missing

    # 4. Verification of dynamic grammar constraints preventing completely empty inputs/outputs
    # If partial reaction outputs are empty (like partial_h2o_missing), the synthesized outputs must not be empty.
    # Let's verify the dynamically generated constraint exists and points to the correct Vector{Molecule}() index.
    pool = [m1, m2, m3]
    grammar = reaction_grammar(pool; settings = settings, partial_reaction = partial_h2o_missing)
    empty_list_idx = findfirst(r -> r == :(Vector{Molecule}()), grammar.rules)
    @test !isnothing(empty_list_idx)

    # Verify that the constraint is indeed Forbidden constraint and is correctly target-indexed
    partial_rule_expr = :(Reaction(vcat(molecule_list, input_molecules), vcat(molecule_list, output_molecules)))
    partial_idx = findfirst(r -> r == partial_rule_expr, grammar.rules)
    @test !isnothing(partial_idx)

    has_empty_constraint = any(grammar.constraints) do constraint
        if constraint isa HerbConstraints.Forbidden
            node = constraint.tree
            return HerbCore.get_rule(node) == partial_idx && length(node.children) == 4 &&
                   node.children[3] isa HerbCore.RuleNode && HerbCore.get_rule(node.children[3]) == empty_list_idx
        end
        return false
    end
    @test has_empty_constraint

    # 5. Synthesis with decoy molecules in the pool
    # Pool has m1(H2), m2(O2), m3(H2O), and m4(CO2)
    # The expected reaction should still be found even with extra molecules in the pool.
    m4 = from_SMILES("O=C=O") # CO2
    reactions_decoy = synthesize_reactions([m1, m2, m3, m4], settings; partial_reaction = partial_h2o_missing)
    @test expected in reactions_decoy

    # 6. Standard reaction synthesis regression check (no partial reaction)
    reactions_std = synthesize_reactions([m1, m2, m3], settings)
    @test expected in reactions_std
end

@testitem "Reaction Synthesizer (Complex partial reactions)" begin
    # Methane combustion: CH4 + 2 O2 -> CO2 + 2 H2O
    ch4 = from_SMILES("[H]-[C](-[H])(-[H])-[H]")
    o2 = from_SMILES("[O]=[O]")
    co2 = from_SMILES("[O]=[C]=[O]")
    h2o = from_SMILES("[H]-[O]-[H]")

    expected = Reaction(nothing, [(1, ch4), (2, o2)], [(1, co2), (2, h2o)], false)
    settings = SynthesizerSettings(max_depth = 11)

    # Missing both O2 on left and H2O on right
    partial_both = Reaction(nothing, [(1, ch4)], [(1, co2)], false)
    reactions_both = synthesize_reactions([ch4, o2, co2, h2o], settings; partial_reaction = partial_both)
    @test expected in reactions_both

    # Missing CO2 on right
    partial_co2 = Reaction(nothing, [(1, ch4), (2, o2)], [(2, h2o)], false)
    reactions_co2 = synthesize_reactions([ch4, o2, co2, h2o], settings; partial_reaction = partial_co2)
    @test expected in reactions_co2
end
