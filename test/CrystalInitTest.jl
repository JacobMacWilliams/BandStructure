module CrystalInitTest
using FromFile
using Test
@from "..\\utils\\CrystalIO.jl" using CrystalIO: writeToBravaisConf, writeToCrystalConf
@from "..\\src\\CrystalLattice.jl" using CrystalLattice
@from "..\\src\\ElectronLattice.jl" using ElectronLattice

bravaisconf = joinpath("conf", "bravais.ini")
if !ispath(bravaisconf)
    bravaisconf = joinpath("conf", "bravais.default.ini")
end

crystalconf = joinpath("conf", "crystal.ini")
if !ispath(crystalconf)
    crystalconf = joinpath("conf", "crystal.default.ini")
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

function electronInit(ename::String, cname::String, states::Int, population::Int, cfile::String, bfile::String)
    ecrystal = ElectronCrystal(ename, cname, states, population, cfile, bfile)
    ename_2 = getName(ecrystal)
    cname_2 = getName(getCrystal(ecrystal))
    states_2 = getStatesPerSite(ecrystal)
    population_2 = getElectronNumber(ecrystal)

    c1 = (ename == ename_2)
    c2 = (cname == cname_2)
    c3 = (states == states_2)
    c4 = (population == population_2)

    return c1 && c2 && c3 && c4
end

# These tests must be carried out in the order as they appear here...
# TODO: Decouple test cases.
@test crystalBravaisInit("testCrystal", "testLattice", 200, rand(Float64, 2, 3), rand(Float64, 2, 2), bravaisconf, crystalconf)
@test electronInit("testEcrystal", "testCrystal", 2, 1000,  crystalconf, bravaisconf)

end