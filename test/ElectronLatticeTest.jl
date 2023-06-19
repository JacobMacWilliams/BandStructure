module ElectronLatticeTest
using FromFile, Test
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/ElectronLattice.jl" using ElectronLattice

function electronInit(efile::String, cfile::String, bfile::String)
    model = "hubbard model"
    crystal = "graphene"
    N = 2

    ecrystal = ElectronCrystal(model, efile, cfile, bfile)
    ename = lowercase(getName(ecrystal))
    cname = lowercase(getName(getCrystal(ecrystal)))
    states= getStatesPerSite(ecrystal)

    c1 = (ename == model)
    c2 = (cname == crystal)
    c3 = (states == N)

    return c1 && c2 && c3
end

function constructHoppingMatrixTest(efile::String, cfile::String, bfile::String)
    model = "hubbard model"
    ecrystal = ElectronCrystal(model, efile, cfile, bfile)
    nnhopping = [1 1 15 5; 1 1 5 15; 1 1 1 1; 1 1 1 1]
    localhopping = [7 1 15 5; 1 7 5 15; 15 5 7 1; 5 15 1 7]
    nnarrays = getPrimitiveNearestNeighbors(getCrystal(ecrystal))
    nnhoppingtest = getHoppingMatrixStructure(ecrystal, 1, nnarrays)
    localhoppingtest = getHoppingMatrixStructure(ecrystal, 0, nnarrays)

    c1 = (nnhopping == nnhoppingtest)
    c2 = (localhopping == localhoppingtest)
    return c1 && c2
end

bravaisconf = joinpath("conf", "bravais.default.toml")
crystalconf = joinpath("conf", "crystal.default.toml")
electronconf = joinpath("conf", "elattice.default.toml")
@test electronInit(electronconf, crystalconf, bravaisconf)
@test constructHoppingMatrixTest(electronconf, crystalconf, bravaisconf)
end