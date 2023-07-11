module CrystalInitTest
using FromFile, Test
using Plots
@from "..\\src\\BravaisLattice.jl" using BravaisLattice
@from "..\\src\\CrystalLattice.jl" using CrystalLattice

function configureGraphene()
    vecs = zeros(2, 2)
    vecs[1,1] = 1/2 * sqrt(3)
    vecs[2,1] = 1/2 
    vecs[1,2] = 1/2 * sqrt(3)
    vecs[2,2] = -1/2
  
    basevecs = zeros(2,2)
    basevecs[1,1] = 0
    basevecs[2,1] = 0
    basevecs[1,2] = 1 / sqrt(3)
    basevecs[2,2] = 0
    return basevecs, vecs
end

function generateLatticeSanityCheck()
    orderstrings = ["first", "second", "third"]
    bravaisconf = joinpath("conf", "bravais.default.toml")
    lattice = Bravais("triangle", bravaisconf)
    for i in 1:3
        points = getLatticePoints(lattice, i)
        plot = scatter(points[1, :], points[2, :])
        pngname = orderstrings[i] * "orderpoints.png"
        savepath = joinpath("plots", pngname)
        savefig(plot, savepath)
    end
end

function crystalBravaisInit(name::String)
    bravaisconf = joinpath("conf", "bravais.default.toml")
    crystalconf = joinpath("conf", "crystal.default.toml")
    crystal = Crystal(name, crystalconf, bravaisconf)
    _name = getName(crystal)
    _basisvecs, _bravaisvecs = getVecs(crystal)
    basisvecs, bravaisvecs = configureGraphene()
    c1 = (basisvecs == _basisvecs)
    c2 = (bravaisvecs == _bravaisvecs)
    c3 = (lowercase(_name) == "graphene")
    
    return c1 && c2 && c3
end

#generateLatticeSanityCheck()
@test crystalBravaisInit("graphene")
end