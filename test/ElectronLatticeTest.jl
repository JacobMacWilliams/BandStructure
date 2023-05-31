module ElectronLatticeTest
using FromFile, Test
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/ElectronLattice.jl" using ElectronLattice

bravaisconf = joinpath("conf", "bravais.default.ini")
crystalconf = joinpath("conf", "crystal.default.ini")

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

function constructHoppingMatrixTest(cfile::String, bfile::String)
    ecrystal = ElectronCrystal("egraphene", "graphene", 2, 1000, cfile, bfile)
    nnhopping = [1 1 15 5; 1 1 5 15; 1 1 1 1; 1 1 1 1]
    localhopping = [7 1 15 5; 1 7 5 15; 15 5 7 1; 5 15 1 7]
    nnarrays = getPrimitiveNearestNeighbors(getCrystal(ecrystal))
    nnhoppingtest = getHoppingMatrixStructure(ecrystal, 1, nnarrays)
    localhoppingtest = getHoppingMatrixStructure(ecrystal, 0, nnarrays)

    c1 = (nnhopping == nnhoppingtest)
    c2 = (localhopping == localhoppingtest)
    return c1 && c2
end

@test electronInit("testEcrystal", "graphene", 2, 1000,  crystalconf, bravaisconf)
@test constructHoppingMatrixTest(crystalconf, bravaisconf)
end