function synthesize_molecules(
        atoms::Vector{Atom}, settings::SynthesizerSettings,
        starting_element::Union{Symbol, AbstractRuleNode} = :molecule; fragment_rules::Vector{Expr} = Expr[],
        starting_fragments::Vector{Expr} = Expr[]
)::Vector{Molecule}
    grammar = SMILES_grammar(atoms; settings = settings, fragment_rules = fragment_rules,
        starting_fragments = starting_fragments)
    iterator = get_iterator(settings, grammar, starting_element)

    candidates = Vector{Molecule}()
    unique_smiles = Set{String}()
    start_time = time()
    for program in iterator
        molecule = interpret_molecule(program, grammar)
        smiles = to_SMILES(molecule)
        mol = RDKitMinimalLib.get_mol(replace(smiles, "≡" => "#"))
        if isnothing(mol) && !(haskey(settings.options, :disable_valid_smiles) &&
           settings.options[:disable_valid_smiles])
            throw("Synthesized invalid molecule: $smiles")
        end
        if haskey(settings.options, :rdkit_unique_candidates) &&
           settings.options[:rdkit_unique_candidates]
            canon_smiles = RDKitMinimalLib.get_smiles(mol)
            if !(canon_smiles in unique_smiles)
                push!(candidates, molecule)
                push!(unique_smiles, canon_smiles)
            end
        else
            push!(candidates, molecule)
        end

        if check_stop_condition(settings, start_time, candidates, molecule)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end
