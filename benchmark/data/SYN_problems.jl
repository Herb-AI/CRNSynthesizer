using CRNSynthesizer

function parse_rxn(rxn_str::SubString{String})::Tuple{Vector{Molecule}, Vector{Molecule}}
    rxn_parts = split(rxn_str, ">>")
    if length(rxn_parts) != 2
        error("Invalid reaction string: $rxn_str")
    end
    reactants_part = strip(rxn_parts[1])
    products_part = strip(rxn_parts[2])

    reactants_smiles = filter(!isempty, map(strip, split(reactants_part, ".")))
    products_smiles = filter(!isempty, map(strip, split(products_part, ".")))

    reactants = [from_SMILES(make_smiles_custom_explicit(make_smiles_rdkit_explicit(String(s))))
                    for s in reactants_smiles]
    products = [from_SMILES(make_smiles_custom_explicit(make_smiles_rdkit_explicit(String(s))))
                for s in products_smiles]

    return reactants, products
end

function parse_syn_problem(reaction_str::String)::ProblemDefinition
    parts = split(reaction_str, ",")
    if length(parts) != 2
        error("Invalid input format: expected two reaction strings separated by a comma")
    end
    incomplete_rxn_str = strip(parts[1])
    target_rxn_str = strip(parts[2])

    incomplete_reactants, incomplete_products = parse_rxn(incomplete_rxn_str)
    target_reactants, target_products = parse_rxn(target_rxn_str)

    # Collect unique molecules from the target reaction
    all_molecules = Molecule[]
    for m in vcat(target_reactants, target_products)
        if !(m in all_molecules)
            push!(all_molecules, m)
        end
    end

    # Construct the target reaction with proper stoichiometry
    reactant_counts = Dict{Molecule, Int}()
    for m in target_reactants
        reactant_counts[m] = get(reactant_counts, m, 0) + 1
    end
    inputs = [(reactant_counts[m], m) for m in unique(target_reactants)]

    product_counts = Dict{Molecule, Int}()
    for m in target_products
        product_counts[m] = get(product_counts, m, 0) + 1
    end
    outputs = [(product_counts[m], m) for m in unique(target_products)]

    reaction = CRNSynthesizer.Reaction(nothing, inputs, outputs, false)
    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Calculate known indices
    incomplete_molecules = Set{Molecule}(vcat(incomplete_reactants, incomplete_products))
    selected_known_indices = Int[]
    for (i, m) in enumerate(all_molecules)
        if m in incomplete_molecules
            push!(selected_known_indices, i)
        end
    end

    known_molecules = all_molecules[selected_known_indices]
    goal_molecules = Molecule[]
    for i in eachindex(all_molecules)
        if !(i in selected_known_indices)
            push!(goal_molecules, all_molecules[i])
        end
    end

    atom_valences = get_valences_from_molecules(all_molecules)

    return ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network
    )
end

function syn_problem()
    DEFAULT_SYN_STR = "CC1(C)OB(c2cn[nH]c2)OC1(C)C.O=[N+]([O-])c1ccc2c(c1)c(Br)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1>>O=[N+]([O-])c1ccc2c(c1)c(-c1cn[nH]c1)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1,CC1(C)OB(c2cn[nH]c2)OC1(C)C.O=[N+]([O-])c1ccc2c(c1)c(Br)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1>>O=[N+]([O-])c1ccc2c(c1)c(-c1cn[nH]c1)nn2C(c1ccccc1)(c1ccccc1)c1ccccc1.CC1(C)OB(Br)OC1(C)C"
    return parse_syn_problem(DEFAULT_SYN_STR)
end
