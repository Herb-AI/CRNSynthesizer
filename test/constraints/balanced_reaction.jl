@testitem "BalancedReaction (from Molecules)" begin
    using HerbGrammar

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Settings for the synthesizer without the BalancedReaction constraint
    settings = SynthesizerSettings(
        max_depth = 5, options = Dict{Symbol, Any}(:disable_balanced_reaction => true)
    )

    # Create a reaction grammar without the BalancedReaction constraint
    without_constraint_grammar = reaction_grammar([m1, m2, m3], settings = settings)

    # Create a reaction grammar with the BalancedReaction constraint
    with_constraint_grammar = reaction_grammar([m1, m2, m3], settings = settings)
    constraint = BalancedReaction(complete_grammar = false)
    addconstraint!(with_constraint_grammar, constraint)

    # Synthesize both options
    iterator = get_iterator(settings, without_constraint_grammar, :reaction)
    interpreter = x -> interpret_reaction(x, without_constraint_grammar)
    without_candidates = Vector{Reaction}()
    find_programs!(iterator, settings, interpreter, without_candidates)

    iterator = get_iterator(settings, with_constraint_grammar, :reaction)
    interpreter = x -> interpret_reaction(x, with_constraint_grammar)
    with_candidates = Vector{Reaction}()
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

@testitem "BalancedReaction and Valences for Ions" begin
    using HerbGrammar
    using DataStructures

    # Parse molecules with ions
    m_na = from_SMILES("[Na+]")
    m_cl = from_SMILES("[Cl-]")
    m_h = from_SMILES("[H+]")
    m_nacl = from_SMILES("[Na]-[Cl]")

    # Test count_atoms charge integration
    counts_na = count_atoms(m_na)
    @test counts_na["[Na+]"] == 1

    counts_cl = count_atoms(m_cl)
    @test counts_cl["[Cl-]"] == 1
    @test !haskey(counts_cl, "charge")

    counts_nacl = count_atoms(m_nacl)
    @test counts_nacl["[Na]"] == 1
    @test counts_nacl["[Cl]"] == 1
    @test !haskey(counts_nacl, "charge")

    # Test dynamic valence extraction
    valences = get_valences_from_molecules([m_na, m_cl, m_nacl])
    @test valences["[Na+]"] == 0
    @test valences["[Cl-]"] == 0
    @test valences["[Na]"] == 1  # in NaCl, it has 1 single bond
    @test valences["[Cl]"] == 1  # in NaCl, it has 1 single bond

    # Test aromatic valences (1.5 bond order)
    m_aromatic = from_SMILES("[C]:1:[C]:[C]:[C]:[C]:[C]:1")
    valences_arom = get_valences_from_molecules([m_aromatic])
    @test valences_arom["[C]"] == 3 # each carbon has two aromatic bonds: 1.5 + 1.5 = 3.0 -> 3

    # Test reaction balancing
    constraint = BalancedReaction()

    # Na+ + Cl- -> NaCl (Unbalanced because Na+ vs Na and Cl- vs Cl are different atom types under new semantics)
    r_unbalanced_nacl = Reaction(
        nothing,
        [(1, m_na), (1, m_cl)],
        [(1, m_nacl)]
    )
    @test !is_valid(r_unbalanced_nacl, constraint)

    # Na+ + Cl- -> Na+ + Cl- (Balanced)
    r_balanced = Reaction(
        nothing,
        [(1, m_na), (1, m_cl)],
        [(1, m_na), (1, m_cl)]
    )
    @test is_valid(r_balanced, constraint)

    # Na+ + Cl- -> NaCl + H+ (Unbalanced)
    r_unbalanced = Reaction(
        nothing,
        [(1, m_na), (1, m_cl)],
        [(1, m_nacl), (1, m_h)]
    )
    @test !is_valid(r_unbalanced, constraint)
end
