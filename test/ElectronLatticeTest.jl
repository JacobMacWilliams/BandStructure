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

bravaisconf = joinpath("conf", "bravais.default.toml")
crystalconf = joinpath("conf", "crystal.default.toml")
electronconf = joinpath("conf", "elattice.default.toml")
@test electronInit(electronconf, crystalconf, bravaisconf)
end