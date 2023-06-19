module MeanFieldTest
# This is not a true test module, it is simply here to test whether or not the code runs.
using FromFile, Test
@from "../src/MeanField.jl" using MeanField
@from "../src/ElectronLattice.jl" using ElectronLattice
@from "../src/CrystalLattice.jl" using CrystalLattice

#=
function initest(cfile, bfile)
    ecrystal = ElectronCrystal("egraphene", "graphene", 2, 500, cfile, bfile)
    nnarrays = getPrimitiveNearestNeighbors(getCrystal(ecrystal))
    ferrigreenfunc = initMeanFieldCorrellator(ecrystal, nnarrays, "ferri")
    aferrigreenfunc = initMeanFieldCorrellator(ecrystal, nnarrays, "aferri")
    randgreenfunc = initMeanFieldCorrellator(ecrystal, nnarrays, "rand")
    println(typeof(ferrigreenfunc))
    println(size(ferrigreenfunc))
    return true
end
=#

function initstructest()
    model = "hubbard model"
    crystalconf = joinpath("conf", "crystal.default.toml")
    bravaisconf = joinpath("conf", "bravais.default.toml")
    econf = joinpath("conf", "elattice.default.toml")
    ecrystal = ElectronCrystal(model, econf, crystalconf, bravaisconf)
    b, _ = getVecs(ecrystal)
    nbasis = size(b, 1)
    pmap = hubbardprimeslicemapgraphene(nbasis)
    #pvmap = Dict([(2, (0.15, 0.75)), (3, (0.15, 0.25))])
    mf = FieldCorrelator(ecrystal, pmap, true)
    initmfstruct!(mf)

    expstruct = [2 3 1 1; 3 2 1 1; 1 1 2 3; 1 1 3 2]
    correlator = getcorrelator(mf)
    return (expstruct == correlator)
end

function initvaluestest()
    model = "hubbard model"
    crystalconf = joinpath("conf", "crystal.default.toml")
    bravaisconf = joinpath("conf", "bravais.default.toml")
    econf = joinpath("conf", "elattice.default.toml")
    ecrystal = ElectronCrystal(model, econf, crystalconf, bravaisconf)
    b, _ = getVecs(ecrystal)
    nbasis = size(b, 1)
    pmap = hubbardprimeslicemapgraphene(nbasis)
    
    mf = FieldCorrelator(ecrystal, pmap, true)
    initmfstruct!(mf)
    pvmap = Dict([(1, 0.0), (2, (0.15, 0.75)), (3, (0.15, 0.25))])
    initmfvalues!(mf, pvmap)
    correlator = getcorrelator(mf)
    println(correlator)
    println(typeof(correlator))
end

@test initstructest()
#initvaluestest()
end