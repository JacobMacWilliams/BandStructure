module ElectronLatticeTest
using FromFile, Test
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/ElectronLattice.jl" using ElectronLattice

bravaisconf = joinpath("conf", "elatticeconf", "bravais.ini")
crystalconf = joinpath("conf", "elatticeconf", "crystal.ini")

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

@test electronInit("testEcrystal", "graphene", 2, 1000,  crystalconf, bravaisconf)
end