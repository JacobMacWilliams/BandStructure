module Embeddings
export embedinlist,
       embedinmatrix,
       extractfromlist,
       extractfromatrix
       
# Assuming the first index of state changes the fastest along the rows
# the second the second fastest along the rows...etc.
function embedinlist(state, sizes)
    idx = 1
    size = 1
    for (i, qnum) in enumerate(state)
        idx += (qnum - 1) * size 
        size *= sizes[i]
    end
    return idx
end

function embedinmatrix(tostate, fromstate, tosizes, fromsizes)
    row = embedinlist(tostate, tosizes)
    col = embedinlist(fromstate, fromsizes)
    return (row, col)
end

function extractfromlist(state::Int, sizes)
    totalstates = prod(sizes)
    if state > totalstates
        error("The state number provided is incompatible with the total" *
              "degrees of freedom of the system.")
    end
    
    state = state - 1
    smallstep = 1
    counted = 0
    
    extstate = zeros(Int, length(sizes))
    for (i, n) in enumerate(sizes)
        
        largestep = n * smallstep
        stepsleft = mod(state, largestep) - counted
        if stepsleft == 0
            smallstep = largestep
            continue
        end

        steps = Int(stepsleft / smallstep)
        extstate[i] = steps
        counted += stepsleft

        if counted == state
            break
        end

        smallstep = largestep
    end

    extstate = extstate .+ 1
    return Tuple(extstate)
end

function extractfromatrix(tostate::Int, fromstate::Int, tosizes, fromsizes)
    extractedtostate = extractfromlist(tostate, tosizes)
    extractedfromstate = extractfromlist(fromstate, fromsizes)
    idx = (extractedtostate..., extractedfromstate...)
    return idx
end

end