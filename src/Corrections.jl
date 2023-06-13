module Corrections
using FromFile
@from "../utils/Embeddings.jl" using Embeddings: embedinmatrix
export gethopdelta,
       getspinflipdelta,
       getpotdelta

    function gethopdelta(tostate, fromstate, sizes, potential::Matrix, correlator::Matrix)
        row, col = embedinmatrix(tostate, fromstate, sizes)
        v = potential[row, col]
        g = conj(correlator[row, col])
        hopdelta = - v * g
        return hopdelta
    end

    function getspinflipdelta(tostate, fromstate, sizes, pairing::Float64, correlator::Matrix)
        if (tostate[1] == fromstate[1]) || (tostate[2:end] != fromstate[2:end])
            error("These states don't participate in spin flips.")
        end
        row, col = embedinmatrix(tostate, fromstate, sizes)
        g = conj(correlator[tostate..., fromstate...])
        spinflipdelta = - g * pairing
        return spinflipdelta
    end

    function getpotdelta(state, sizes, pairing::Float64, potential::Matrix, correlator::Matrix)
        if state[end] != 1
            error("This function can only calculate potential terms for sites located at the origin.")
        end

        spin = state[1] # Assuming labels 1, 2
        site = state[2:end]

        oppositespin = mod(spin + 1, 2) == 0 ? 2 : 1
        oppositespinstate = [oppositespin site[1] site[2]]
        pairingterm = pairing * correlator[oppositespinstate..., oppositespinstate...]

        potentialterm = 0
        sitesizes = sizes[2:end]
        for i in CartesianIndices(potential)
            if Tuple(site) == i
                continue
            end

            siterow, sitecol = embedinmatrix(i, site, sitesizes)
            vij = potential[siterow, sitecol]

            for s in 1:sizes[1]
                element = [(s, i[2], 1) (s, i[2], 1)]
                srow, scol = embedinmatrix(element..., sizes)
                potentialterm += vij * correlator[srow, scol]
            end
        end

        return pairingterm + potentialterm
    end

end