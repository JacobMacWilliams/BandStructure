module EmbeddingsTest
using FromFile, Test
@from "../utils/Embeddings.jl" using Embeddings

function embeddingextractiontest()

    # The iterator iter corresponds to the iterator used to index a matrix of dimension 2x3x2x2x3x2.
    tostates = Iterators.product(1:2, 1:3, 1:2)
    fromstates = tostates
    iter = Iterators.product(tostates, fromstates)
    sizes = (2, 3, 2)
    embeddedtest = zeros(prod(sizes), prod(sizes))

    # We now transform this iterator into that of a (2*3*2)x(2*3*2) matrix.
    sum = 1
    for (tostate, fromstate) in iter
        rowcol = embedinmatrix(tostate, fromstate, sizes)
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
    embeddedexp = zeros(prod(sizes), prod(sizes))
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
    recovedidx = []
    for i in CartesianIndices(embeddedtest)
        ij = Tuple(i)
        push!(recovedidx, extractfromatrix(ij[1], ij[2], sizes))
    end

    expidx = collect(Iterators.product(tostates, fromstates))
    for i in eachindex(expidx)
        cond = (recovedidx[i] == expidx[i])
        if (!cond)
            return false
        end
    end

    return true
end

@test embeddingextractiontest()
end