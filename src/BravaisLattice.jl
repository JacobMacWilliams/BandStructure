module BravaisLattice
using FromFile
@from "../utils/CrystalIO.jl" using CrystalIO: parseLatticeConf as getLattice
export Bravais, getName, getDimension, getVecs

struct Bravais
    name::String
    dimension::UInt8
    vecs::Matrix{Float64}
end

function Bravais(name::String, file::String)
    dim, vecs = getLattice(file, name)
    Bravais(name, dim, vecs)
end

function getName(lattice::Bravais)
    return lattice.name
end

function getDimension(lattice::Bravais)
    return lattice.dimension
end

function getVecs(lattice::Bravais)
    return lattice.vecs
end

end
