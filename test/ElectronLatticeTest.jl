module ElectronLatticeTest
using FromFile, Test
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/ElectronLattice.jl" using ElectronLattice

bravaisconf = joinpath("conf", "bravais.default.toml")
crystalconf = joinpath("conf", "crystal.default.toml")

function electronInit(ename::String, cname::String, states::Int, cfile::String, bfile::String)
    ecrystal = ElectronCrystal(ename, cname, states, cfile, bfile)
    ename_2 = getName(ecrystal)
    cname_2 = getName(getCrystal(ecrystal))
    states_2 = getStatesPerSite(ecrystal)

    c1 = (ename == ename_2)
    c2 = (cname == cname_2)
    c3 = (states == states_2)

    return c1 && c2 && c3
end

function constructHoppingMatrixTest(cfile::String, bfile::String)
    ecrystal = ElectronCrystal("egraphene", "graphene", 2, cfile, bfile)
    nnhopping = [1 1 15 5; 1 1 5 15; 1 1 1 1; 1 1 1 1]
    localhopping = [7 1 15 5; 1 7 5 15; 15 5 7 1; 5 15 1 7]
    nnarrays = getPrimitiveNearestNeighbors(getCrystal(ecrystal))
    nnhoppingtest = getHoppingMatrixStructure(ecrystal, 1, nnarrays)
    localhoppingtest = getHoppingMatrixStructure(ecrystal, 0, nnarrays)

    c1 = (nnhopping == nnhoppingtest)
    c2 = (localhopping == localhoppingtest)
    return c1 && c2
end

@test electronInit("testEcrystal", "graphene", 2, crystalconf, bravaisconf)
@test constructHoppingMatrixTest(crystalconf, bravaisconf)
end