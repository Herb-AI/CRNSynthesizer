using CRNSynthesizer
import MoleculeFlow

"""
    is_feasible_problem(problem::ProblemDefinition) -> Tuple{Bool, String, Int, Vector{Tuple{String, Int}}, Vector{Int}}

Classify a synthesis problem as feasible or infeasible.
Returns (is_feasible, reason_string, brics_count, small_molecules_list, unfeasible_goal_sizes).
"""
function is_feasible_problem(problem::ProblemDefinition)::Tuple{
        Bool, String, Int, Vector{Tuple{String, Int}}, Vector{Int}}
    # Decompose all known molecules into BRICS fragments
    known_frags = Set{String}()
    for m in problem.known_molecules
        frags = get_brics_cached(m.canonical_smiles)
        if !ismissing(frags)
            union!(known_frags, frags)
        end
    end

    small_enough_mols = Tuple{String, Int}[]
    brics_count = 0
    unfeasible_goal_sizes = Int[]
    first_failure_reason = ""

    # Check each goal molecule
    for goal_mol in problem.goal_molecules
        total_atoms = length(goal_mol.atoms)

        # A missing molecule is feasible if it has size <= 6 
        # small size due to increased branching factor after introduction of fragments
        if total_atoms <= 6
            push!(small_enough_mols, (goal_mol.canonical_smiles, total_atoms))
            continue
        end

        # Or if it shares a BRICS fragment (i.e. has a substructure match with any known BRICS fragment)
        goal_m = get_mol_cached(goal_mol.canonical_smiles)
        if ismissing(goal_m)
            return false,
            "Goal molecule $(goal_mol.canonical_smiles) could not be parsed by MoleculeFlow",
            0,
            Tuple{String, Int}[],
            Int[]
        end

        goal_atom_count = get_mol_atom_count_cached(goal_mol.canonical_smiles)

        shared_any = false
        for frag_smiles in known_frags
            # Pre-filter by atom count
            frag_atom_count = get_frag_atom_count_cached(frag_smiles)
            if frag_atom_count > goal_atom_count
                continue
            end

            cleaned_m = get_cleaned_frag_cached(frag_smiles)
            if !ismissing(cleaned_m)
                match_res = MoleculeFlow.has_substructure_match(goal_m, cleaned_m)
                if !ismissing(match_res) && match_res
                    shared_any = true
                    break
                end
            end
        end

        if !shared_any
            push!(unfeasible_goal_sizes, total_atoms)
            if isempty(first_failure_reason)
                first_failure_reason = "Goal molecule has > 6 atoms ($total_atoms) and does not share a BRICS fragment"
            end
        else
            brics_count += 1
        end
    end

    if !isempty(unfeasible_goal_sizes)
        return false,
        first_failure_reason,
        0,
        Tuple{String, Int}[],
        unfeasible_goal_sizes
    end

    small_enough_names = [sm for (sm, _) in small_enough_mols]
    small_enough_str = isempty(small_enough_mols) ? "0" :
                       "$(length(small_enough_mols)) ($(join(small_enough_names, ", ")))"
    return true,
    "$small_enough_str goal molecules are small enough, $brics_count share a BRICS fragment",
    brics_count,
    small_enough_mols,
    unfeasible_goal_sizes
end
