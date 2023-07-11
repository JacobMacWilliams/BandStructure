module BravaisLattice
using FromFile
@from "../utils/TOMLIO.jl" using TOMLIO
export Bravais, getName, getDimension, getVecs, getLatticePoints

const CLASSKEY = "bravais"

struct Bravais
    name::String
    dimension::UInt8
    vecs::Matrix{Float64}
end

function Bravais(name::String, file::String)
    config = parseconfig(file, CLASSKEY, name)
    dim = get(config, "dim", nothing)
    vecs = parsematrix(config, "v")
    Bravais(name, dim, vecs)
end

#=
function Bravais(name::String, file::String)
    dim, vecs = getLattice(file, name)
    Bravais(name, dim, vecs)
end
=#

function getName(lattice::Bravais)
    return lattice.name
end

function getDimension(lattice::Bravais)
    return lattice.dimension
end

function getVecs(lattice::Bravais)
    return lattice.vecs
end

function getLatticePoints(lattice::Bravais, order)
    bravais = getVecs(lattice)
    dim = size(bravais, 1)
    relvecs = hcat(-bravais, bravais)
    latticepoints = zeros(dim)
    
    for i = 1:order
        candidatepoints = []
        for point in eachcol(latticepoints)
            surroundingpoints = hcat(Iterators.map(v -> point + v, eachcol(relvecs))...)
            if size(candidatepoints, 1) == 0
                candidatepoints = surroundingpoints
                continue
            end
            candidatepoints = hcat(candidatepoints, surroundingpoints)
        end
        latticepoints = hcat(unique(eachcol(hcat(latticepoints, candidatepoints)))...)
    end
    return latticepoints
end

end
