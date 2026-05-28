@testitem "ContainsReactions" begin
    using HerbGrammar
    using DataStructures

    # Create some reactions
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # TODO: Create balanced reactions
    r1 = Reaction([m1, m2], [m3])
    r2 = Reaction([m3], [m1, m2])
    r3 = Reaction([m1], [m2, m3])

    # Settings for the synthesizer without the ContainsReactions constraint
    settings = SynthesizerSettings(
        max_depth = 7, options = Dict{Symbol, Any}(
            :disable_contains_reactions => true, :disable_contains_molecules => true)
    )

    # Define a problem to pass required molecules and atom valences
    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2)
    problem = ProblemDefinition(
        atom_valences,
        Molecule[],
        Molecule[],
        ReactionNetwork(Reaction[]);
        required_molecules = [RequiredMolecule(m1, INPUT), RequiredMolecule(m2, OUTPUT)]
    )

    # Create a network grammar without the ContainsReactions constraint
    without_constraint_grammar = network_grammar([r1, r2, r3], problem; settings = settings)

    # Create a network grammar with the ContainsReactions constraint
    with_constraint_grammar = network_grammar([r1, r2, r3], problem; settings = settings)
    constraint = ContainsReactions(
        with_constraint_grammar,
        [RequiredMolecule(m1, INPUT), RequiredMolecule(m2, OUTPUT)],
        [r1, r2, r3]
    )
    addconstraint!(with_constraint_grammar, constraint)

    # Synthesize both options
    iterator = get_iterator(settings, without_constraint_grammar, :network)
    interpreter = x -> interpret_network(x, without_constraint_grammar)
    without_candidates = Vector{ReactionNetwork}()
    find_programs!(iterator, settings, interpreter, without_candidates)

    iterator = get_iterator(settings, with_constraint_grammar, :network)
    interpreter = x -> interpret_network(x, with_constraint_grammar)
    with_candidates = Vector{ReactionNetwork}()
    find_programs!(iterator, settings, interpreter, with_candidates)

    @test length(without_candidates) > 0
    @test length(with_candidates) > 0
    @test length(with_candidates) < length(without_candidates)
    valid_candidates = filter(x -> is_valid(x, constraint), without_candidates)
    @test length(valid_candidates) == length(with_candidates)
    @test all(x -> is_valid(x, constraint), with_candidates)
    @test all(x -> x in valid_candidates, with_candidates)
    @test all(x -> x in with_candidates, valid_candidates)
end
