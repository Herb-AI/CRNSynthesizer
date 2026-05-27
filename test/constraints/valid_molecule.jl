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

@testitem "Fragment Trailing Ringbonds Duplicate Digit Prevention" begin
    using HerbGrammar
    using DataStructures
    using HerbSearch

    # Create atom valences
    atom_valences = OrderedDict("[H]" => 1, "[C]" => 4)

    # Define starting fragment expression with matching trailing ringbonds
    starting_fragments = [:(starting_fragment = "[C]-1" * fragment_16_exit * "=[C](-[H])-[C](-[H])=[C](-[H])-[C](-[H])=[C]-1" * fragment_16_exit)]

    # Settings for the synthesizer
    settings = SynthesizerSettings(
        max_depth = 5,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )

    # Build the SMILES grammar with the starting fragment
    grammar = SMILES_grammar(
        atom_valences;
        settings = settings,
        starting_fragments = starting_fragments
    )

    # Synthesize molecules
    iterator = get_iterator(settings, grammar, :molecule)
    interpreter = x -> interpret_molecule(x, grammar)
    candidates = Vector{Molecule}()
    find_programs!(iterator, settings, interpreter, candidates)

    # Ensure that candidates are successfully found
    @test length(candidates) > 0
    
    # Ensure all synthesized candidates are valid Molecules
    for cand in candidates
        @test !isempty(cand.canonical_smiles)
    end
end

