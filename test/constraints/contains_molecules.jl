@testitem "ContainsMolecules" begin
    using HerbGrammar
    using DataStructures

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Settings for the synthesizer without the ContainsMolecules constraint
    settings = SynthesizerSettings(
        max_depth = 7, options = Dict{Symbol, Any}(:disable_contains_molecules => true)
    )

    # Define a problem to pass required molecules and atom valences
    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2)
    problem = ProblemDefinition(
        atom_valences,
        Molecule[],
        Molecule[],
        ReactionNetwork(Reaction[])
    )

    # Create a network grammar without the ContainsMolecules constraint
    without_constraint_grammar = network_grammar([m1, m2, m3], problem, settings = settings)

    # Create a network grammar with the ContainsMolecules constraint
    with_constraint_grammar = network_grammar(
        [m2, m3], problem, settings = settings, check_required = false
    )
    constraint = ContainsMolecules(with_constraint_grammar, [(m1, INPUT)])
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
