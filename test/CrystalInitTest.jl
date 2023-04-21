module CrystalInitTest
using FromFile
using Test
@from "..\\utils\\CrystalIO.jl" using CrystalIO: writeToBravaisConf, writeToCrystalConf
@from "..\\src\\CrystalLattice.jl" using CrystalLattice

bravaisconf = joinpath("conf", "bravais.ini")
if !ispath(bravaisconf)
    bravaisconf = joinpath("conf", "bravais.default.ini")
end

crystalconf = joinpath("conf", "crystal.ini")
if !ispath(crystalconf)
    crystalconf = joinpath("conf", "bravais.default.ini")
end

function crystalBravaisInit(name::String, lattice::String, N::Integer, basisvecs::Matrix{Float64}, vecs::Matrix{Float64}, bfile::String, cfile::String)
    open(bfile, truncate=true)
    writeToBravaisConf(lattice, vecs, bfile)
    open(cfile, truncate=true)
    writeToCrystalConf(name, N, lattice, basisvecs, cfile)
    crystal = Crystal(name, cfile, bfile)
    c1 = (name == getName(crystal))
    c2 = (lattice == getName(getLattice(crystal)))
    c3 = (N == getSize(crystal))
    parsedBasisVecs, parsedVecs = getVecs(crystal)
    c4 = (basisvecs == parsedBasisVecs)
    c5 = (vecs == parsedVecs)
    return c1 && c2 && c3 && c4 && c5
end

@test crystalBravaisInit("testCrystal", "testLattice", 200, rand(Float64, 2, 3), rand(Float64, 2, 2), bravaisconf, crystalconf)

end