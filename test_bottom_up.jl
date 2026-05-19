using CRNSynthesizer

function compare_approaches()
    atoms = [Atom("H"), Atom("O"), Atom("C")]
    max_depth = 8
    
    # Default (Top-Down)
    println("Running Default (Top-Down) with max_depth = $max_depth...")
    settings_td = SynthesizerSettings(
        max_depth = max_depth,
        options = Dict{Symbol, Any}(:unique_candidates => true)
    )
    start_td = time()
    unique_molecules = synthesize_molecules(atoms, settings_td)
    end_td = time()
    println("Top-Down took: ", end_td - start_td, " seconds")
    println("Number of unique molecules (Top-Down): ", length(unique_molecules))
    max_depth += 2
    # Bottom-Up
    println("\nRunning Bottom-Up with max_depth = $max_depth...")
    settings_bu = SynthesizerSettings(
        max_depth = max_depth,
        iterator = BottomUp,    
        options = Dict{Symbol, Any}(:unique_candidates => true, :use_bottom_up => true)
    )
    start_bu = time()
    bottom_up_molecules = synthesize_molecules(atoms, settings_bu)
    end_bu = time()
    println("Bottom-Up took: ", end_bu - start_bu, " seconds")
    println("Number of unique molecules (Bottom-Up): ", length(bottom_up_molecules))

    # Compare sets based on canonical SMILES
    td_smiles = Set(m.canonical_smiles for m in unique_molecules if !isnothing(m.canonical_smiles))
    bu_smiles = Set(m.canonical_smiles for m in bottom_up_molecules if !isnothing(m.canonical_smiles))

    println("\nMolecules in Top-Down but NOT in Bottom-Up:")
    diff_td = setdiff(td_smiles, bu_smiles)
    if isempty(diff_td)
        println("None")
    else
        # Sort by length to find minimal examples
        sorted_diff_td = sort(collect(diff_td), by=length)
        for s in sorted_diff_td
            println(s)
        end
    end
    
    println("\nMolecules in Bottom-Up but NOT in Top-Down:")
    diff_bu = setdiff(bu_smiles, td_smiles)
    if isempty(diff_bu)
        println("None")
    else
        sorted_diff_bu = sort(collect(diff_bu), by=length)
        for s in sorted_diff_bu
            println(s)
        end
    end
end

compare_approaches()
