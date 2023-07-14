module EmbeddingsTest
using FromFile, Test
@from "../utils/Embeddings.jl" using Embeddings

function embeddingextractiontest()

    # The iterator iter corresponds to the iterator used to index a matrix of dimension 2x3x2x2x3x2.
    iter = Iterators.product(1:2, 1:3, 1:2, 1:2, 1:3, 1:2)
    tosizes = (2, 3, 2)
    fromsizes = tosizes
    embeddedtest = zeros(prod(fromsizes), prod(tosizes))

    # We now transform this iterator into that of a (2*3*2)x(2*3*2) matrix.
    sum = 1
    for (i, j, k, l, m, n) in iter
        rowcol = embedinmatrix((i, j, k), (l, m, n), tosizes, fromsizes)
        embeddedtest[rowcol[1], rowcol[2]] = sum
        sum += 1
    end

    #=
    Checking whether the transformation was done correctly corresponds to ensuring
    that the mapping from the 2x3x2x2x3x2 iterator to the (2*3*2)x(2*3*2) iterator
    is simply T(iter[i]) == T(iter)[i] (if we transformed iter to an array!). If
    this relation holds then the transformation is well behaved. In this case the 
    sum keeps track of the index of the untransformed iterator.
    =#
    embeddedexp = zeros(prod(fromsizes), prod(tosizes))
    sum = 1
    for i in eachindex(embeddedexp)
        embeddedexp[i] = sum
        sum += 1
    end

    if embeddedexp != embeddedtest
        return false
    end

    # Now we turn to extracting the original 2x3x2x2x3x2 iterator from
    # the (2*3*2)x(2x3x2) iterator.
    recoveredidx = []
    for i in CartesianIndices(embeddedtest)
        ij = Tuple(i)
        idx = extractfromatrix(ij[1], ij[2], tosizes, fromsizes)
        push!(recoveredidx, idx)
        println(idx)
    end

    expidx = collect(Iterators.product(1:2, 1:3, 1:2, 1:2, 1:3, 1:2))
    for i in eachindex(expidx)

        cond = (recoveredidx[i] == expidx[i])
        if (!cond)
            return false
        end
    end

    return true
end

@test embeddingextractiontest()
end