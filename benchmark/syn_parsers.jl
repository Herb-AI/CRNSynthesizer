using CRNSynthesizer
using DataStructures
using LRUCache

# -------------------------------------------------------------
# Caching wrappers to avoid expensive RDKit/MoleculeFlow calls
# -------------------------------------------------------------
const CACHE_SIZE = 1024

const PARSED_MOLECULE_CACHE = LRU{String, Molecule}(maxsize = CACHE_SIZE)

function get_parsed_molecule_cached(smiles::String)::Molecule
    get!(PARSED_MOLECULE_CACHE, smiles) do
        from_SMILES(make_smiles_custom_explicit(make_smiles_rdkit_explicit(smiles)))
    end
end

function parse_rxn_cached(rxn_str::SubString{String})::Tuple{
        Vector{Molecule}, Vector{Molecule}}
    # Pre-process to convert free hydrogen/oxygen representations to stable molecules
    # We use lookarounds to ensure we only replace isolated atoms, making an exception
    # if the bond is already explicit (e.g. [H][H] or [H]-[H] or [O]=[O]).
    processed_rxn = replace(rxn_str, r"(?<=^|\.|>>)(\[H\])\.(\[H\])(?=$|\.|>>)" => s"\1-\2")
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)(\[O\])\.(\[O\])(?=$|\.|>>)" => s"\1=\2")

    # Replace any leftover standalone [H] and [O]
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)\[H\](?=$|\.|>>)" => "[H]-[H]")
    processed_rxn = replace(processed_rxn, r"(?<=^|\.|>>)\[O\](?=$|\.|>>)" => "[O]=[O]")

    rxn_parts = split(processed_rxn, ">>")
    if length(rxn_parts) != 2
        error("Invalid reaction string: $rxn_str")
    end
    reactants_part = strip(rxn_parts[1])
    products_part = strip(rxn_parts[2])

    reactants_smiles = filter(!isempty, map(strip, split(reactants_part, ".")))
    products_smiles = filter(!isempty, map(strip, split(products_part, ".")))

    reactants = [get_parsed_molecule_cached(String(s)) for s in reactants_smiles]
    products = [get_parsed_molecule_cached(String(s)) for s in products_smiles]

    return reactants, products
end

function parse_syn_problem(reaction_str::AbstractString)::ProblemDefinition
    parts = split(reaction_str, ",")
    if length(parts) != 2
        error("Invalid input format: expected two reaction strings separated by a comma")
    end
    incomplete_rxn_str = strip(parts[1])
    target_rxn_str = strip(parts[2])

    incomplete_reactants, incomplete_products = parse_rxn_cached(incomplete_rxn_str)
    target_reactants, target_products = parse_rxn_cached(target_rxn_str)

    # Collect unique molecules from the target reaction
    all_molecules = Molecule[]
    for m in vcat(target_reactants, target_products)
        if !(m in all_molecules)
            push!(all_molecules, m)
        end
    end

    # Construct the target reaction with proper stoichiometry
    reactant_counts = OrderedDict{Molecule, Int}()
    for m in target_reactants
        reactant_counts[m] = get(reactant_counts, m, 0) + 1
    end
    inputs = [(reactant_counts[m], m) for m in unique(target_reactants)]

    product_counts = OrderedDict{Molecule, Int}()
    for m in target_products
        product_counts[m] = get(product_counts, m, 0) + 1
    end
    outputs = [(product_counts[m], m) for m in unique(target_products)]

    reaction = CRNSynthesizer.Reaction(nothing, inputs, outputs, false)
    goal_network = CRNSynthesizer.ReactionNetwork([reaction])

    # Calculate known indices
    incomplete_molecules = OrderedSet{Molecule}(vcat(incomplete_reactants, incomplete_products))
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

    inc_reactant_counts = OrderedDict{Molecule, Int}()
    for m in incomplete_reactants
        inc_reactant_counts[m] = get(inc_reactant_counts, m, 0) + 1
    end
    partial_inputs = [(inc_reactant_counts[m], m) for m in unique(incomplete_reactants)]

    inc_product_counts = OrderedDict{Molecule, Int}()
    for m in incomplete_products
        inc_product_counts[m] = get(inc_product_counts, m, 0) + 1
    end
    partial_outputs = [(inc_product_counts[m], m) for m in unique(incomplete_products)]

    partial_reaction = if isempty(partial_inputs) && isempty(partial_outputs)
        nothing
    else
        CRNSynthesizer.Reaction(nothing, partial_inputs, partial_outputs, false)
    end

    return ProblemDefinition(
        atom_valences,
        known_molecules,
        goal_molecules,
        goal_network;
        partial_reaction = partial_reaction
    )
end
