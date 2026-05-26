@testitem "ValidSMILES" begin
    using HerbGrammar
    using DataStructures

    # Settings for the synthesizer
    settings = SynthesizerSettings(
        max_depth = 4, options = Dict{Symbol, Any}(:disable_valid_smiles => true)
    )

    atom_valences = OrderedDict("[H]" => 1, "[O]" => 2, "[C]" => 4)

    # Create a network grammar without the ValidSMILES constraint
    without_constraint_grammar = SMILES_grammar(atom_valences; settings = settings)
    without_constraint_grammar.constraints

    # Create a network grammar with the ValidSMILES constraint
    with_constraint_grammar = SMILES_grammar(atom_valences; settings = settings)
    constraint = ValidSMILES(with_constraint_grammar, atom_valences)
    addconstraint!(with_constraint_grammar, constraint)

    # Synthesize both options
    iterator = get_iterator(settings, without_constraint_grammar, :molecule)
    interpreter = x -> interpret_molecule(x, without_constraint_grammar)
    without_candidates = Vector{Molecule}()
    find_programs!(iterator, settings, interpreter, without_candidates)

    iterator = get_iterator(settings, with_constraint_grammar, :molecule)
    interpreter = x -> interpret_molecule(x, with_constraint_grammar)
    with_candidates = Vector{Molecule}()
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
