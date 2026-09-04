@testitem "Fragment Trailing Ringbonds Duplicate Digit Prevention" begin
    using HerbGrammar
    using DataStructures
    using HerbSearch

    # Create atom valences
    atom_valences = OrderedDict("[H]" => 1, "[C]" => 4)

    # Define starting fragment expression with matching trailing ringbonds
    starting_fragments = [:(starting_fragment = "[C]-1" * fragment_16_exit *
                                                "=[C](-[H])-[C](-[H])=[C](-[H])-[C](-[H])=[C]-1" *
                                                fragment_16_exit)]

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
