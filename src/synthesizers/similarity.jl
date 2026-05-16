
function simpson_similarity(fp1::Vector{UInt8}, fp2::Vector{UInt8})
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

function sort_molecules_by_similarity(molecules::Vector{Molecule}, known_molecules::Vector{Molecule})
    isempty(known_molecules) && return molecules

    scores = zeros(Float64, length(molecules))
    for (i, m) in enumerate(molecules)
        max_score = 0.0
        for km in known_molecules
            score = simpson_similarity(m.fingerprint, km.fingerprint)
            max_score = max(max_score, score)
        end
        scores[i] = max_score
    end

    return molecules[sortperm(scores, rev=true)]
end

function sort_reactions_by_similarity(reactions::Vector{Reaction}, known_molecules::Vector{Molecule})
    isempty(known_molecules) && return reactions

    scores = zeros(Float64, length(reactions))
    for (i, r) in enumerate(reactions)
        for (_, m) in Iterators.flatten((r.inputs, r.outputs))
            mol_max_score = 0.0
            for km in known_molecules
                score = simpson_similarity(m.fingerprint, km.fingerprint)
                mol_max_score = max(mol_max_score, score)
            end
            scores[i] += mol_max_score
        end
    end

    return reactions[sortperm(scores, rev=true)]
end
