function simpson_similarity(m1::Molecule, m2::Molecule)
    fp1 = m1.fingerprint
    fp2 = m2.fingerprint
    intersection = 0
    count1 = 0
    count2 = 0

    for (b1, b2) in zip(fp1, fp2)
        intersection += count_ones(b1 & b2)
        count1 += count_ones(b1)
        count2 += count_ones(b2)
    end

    min_count = min(count1, count2)
    return iszero(min_count) ? 0.0 : intersection / min_count
end

function tanimoto_similarity(m1::Molecule, m2::Molecule)
    fp1 = m1.morgan_fingerprint
    fp2 = m2.morgan_fingerprint
    intersection = 0
    count1 = 0
    count2 = 0

    for (b1, b2) in zip(fp1, fp2)
        intersection += count_ones(b1 & b2)
        count1 += count_ones(b1)
        count2 += count_ones(b2)
    end

    union_count = count1 + count2 - intersection
    return iszero(union_count) ? 0.0 : intersection / union_count
end

function combined_similarity(m1::Molecule, m2::Molecule)
    tan = tanimoto_similarity(m1, m2)
    simp = simpson_similarity(m1, m2)
    return 0.5 * tan + 0.5 * simp
end

function get_similarity_metric(metric_name::Symbol)
    if metric_name == :tanimoto
        return tanimoto_similarity
    elseif metric_name == :both
        return combined_similarity
    elseif metric_name == :simpson
        return simpson_similarity
    else
        return (m1, m2) -> 0.0
    end
end

function sort_molecules_by_similarity(
        molecules::Vector{Molecule},
        known_molecules::Vector{Molecule};
        metric_name::Symbol = :simpson,
        cache::Union{Nothing, Dict{Molecule, Float64}} = nothing
)
    (isempty(known_molecules) || metric_name == :none) && return molecules

    metric = get_similarity_metric(metric_name)

    scores = map(
        m -> get_molecule_similarity_score(m, known_molecules, metric, cache), molecules)

    return molecules[sortperm(scores, rev = true)]
end

function get_molecule_similarity_score(
        m::Molecule,
        known_molecules::Vector{Molecule},
        metric::Function,
        molecule_cache::Union{Nothing, Dict{Molecule, Float64}}
)
    if !isnothing(molecule_cache) && haskey(molecule_cache, m)
        return molecule_cache[m]
    end
    mol_max = maximum(
        (metric(m, km) for km in known_molecules),
        init = 0.0
    )
    if !isnothing(molecule_cache)
        molecule_cache[m] = mol_max
    end
    return mol_max
end

function score_reaction_by_similarity(
        reaction::Reaction,
        known_molecules::Vector{Molecule};
        metric::Function = tanimoto_similarity,
        molecule_cache::Union{Nothing, Dict{Molecule, Float64}} = nothing
)
    total_score = 0.0
    for (_, m) in Iterators.flatten((reaction.inputs, reaction.outputs))
        total_score += get_molecule_similarity_score(
            m, known_molecules, metric, molecule_cache)
    end
    n_total = max(length(reaction.inputs) + length(reaction.outputs), 1)
    return total_score / n_total
end

function sort_reactions_by_similarity(
        reactions::Vector{Reaction},
        known_molecules::Vector{Molecule};
        metric_name::Symbol = :simpson,
        cache::Union{Nothing, Dict{Reaction, Float64}} = nothing,
        molecule_cache::Union{Nothing, Dict{Molecule, Float64}} = nothing
)
    (isempty(known_molecules) || metric_name == :none) && return reactions

    metric = get_similarity_metric(metric_name)
    scores = zeros(Float64, length(reactions))
    for (i, r) in enumerate(reactions)
        if !isnothing(cache) && haskey(cache, r)
            scores[i] = cache[r]
        else
            score = score_reaction_by_similarity(r, known_molecules; metric = metric,
                molecule_cache = molecule_cache)
            scores[i] = score
            if !isnothing(cache)
                cache[r] = score
            end
        end
    end

    return reactions[sortperm(scores, rev = true)]
end
