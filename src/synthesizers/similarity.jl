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

function tango_similarity(m1::Molecule, m2::Molecule)
    tan = tanimoto_similarity(m1, m2)
    simp = simpson_similarity(m1, m2)
    return 0.5 * tan + 0.5 * simp
end

function get_similarity_metric(metric_name::Symbol)
    if metric_name == :tanimoto
        return tanimoto_similarity
    elseif metric_name == :tango
        return tango_similarity
    elseif metric_name == :simpson
        return simpson_similarity
    else
        return (m1, m2) -> 0.0
    end
end

function sort_molecules_by_similarity(molecules::Vector{Molecule}, known_molecules::Vector{Molecule}; metric_name::Symbol = :simpson)
    isempty(known_molecules) && return molecules

    metric = get_similarity_metric(metric_name)
    scores = zeros(Float64, length(molecules))
    for (i, m) in enumerate(molecules)
        max_score = 0.0
        for km in known_molecules
            score = metric(m, km)
            max_score = max(max_score, score)
        end
        scores[i] = max_score
    end

    return molecules[sortperm(scores, rev=true)]
end

function score_reaction_by_similarity(
    reaction::Reaction,
    known_molecules::Vector{Molecule};
    metric::Function = tanimoto_similarity,
    combine::Symbol = :multiplicative  # :multiplicative, :additive, :pooled, or :sum
)
    if combine == :pooled
        total_score = 0.0
        for (_, m) in Iterators.flatten((reaction.inputs, reaction.outputs))
            mol_max = maximum(
                (metric(m, km) for km in known_molecules),
                init=0.0
            )
            total_score += mol_max
        end
        n_total = max(length(reaction.inputs) + length(reaction.outputs), 1)
        return total_score / n_total
    end

    input_score = 0.0
    for (_, m) in reaction.inputs
        mol_max = maximum(
            (metric(m, km) for km in known_molecules),
            init=0.0
        )
        input_score += mol_max
    end

    output_score = 0.0
    for (_, m) in reaction.outputs
        mol_max = maximum(
            (metric(m, km) for km in known_molecules),
            init=0.0
        )
        output_score += mol_max
    end

    if combine == :multiplicative
        # Normalize by count to avoid biasing toward reactions with more molecules
        n_in = max(length(reaction.inputs), 1)
        n_out = max(length(reaction.outputs), 1)
        return (input_score / n_in) * (output_score / n_out)
    elseif combine == :additive
        n_in = max(length(reaction.inputs), 1)
        n_out = max(length(reaction.outputs), 1)
        return (input_score / n_in) + (output_score / n_out)
    else # :sum
        return input_score + output_score # The original unnormalized sum
    end
end

function sort_reactions_by_similarity(reactions::Vector{Reaction}, known_molecules::Vector{Molecule}; metric_name::Symbol = :simpson, combine::Symbol = :multiplicative)
    isempty(known_molecules) && return reactions

    metric = get_similarity_metric(metric_name)
    scores = zeros(Float64, length(reactions))
    for (i, r) in enumerate(reactions)
        scores[i] = score_reaction_by_similarity(r, known_molecules; metric=metric, combine=combine)
    end

    return reactions[sortperm(scores, rev=true)]
end
