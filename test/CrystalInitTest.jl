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

    for i in 1:3

        ordername = orderstrings[i]
        latticename = "triangle-" * ordername * "-order"

        lattice = Bravais(latticename, bravaisconf)
        latticepoints = getlatticepoints(lattice)
        plot = scatter(latticepoints[1, :], latticepoints[2, :])

        pngname = ordername * "orderpoints.png"
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

function nearestNeighborsGrapheneSanityCheck()
    bravaisconf = joinpath("conf", "bravais.default.toml")
    crystalconf = joinpath("conf", "crystal.default.toml")
    crystal = Crystal("graphene", crystalconf, bravaisconf)

    # Only the unitcell at the origin an all surrounding unitcells are generated.
    # The nearest neighbors of each of the atoms of the unit cell at the origin should
    # be contained in this set.
    latticepoints, crystalpoints = getPoints(crystal)
    basisvecs, _ = getVecs(crystal)
    natoms = size(basisvecs, 2)

    # Each atom of the unit cell has exactly 3 nearest neighbors. The set of points in which
    # we are searching for the nearest neighbors of a given site contains the site itself. since
    # each site will register as its own nearest neighbor we have to perform a 4-nn search.
    _, dist = getNearestNeighbors(crystal, crystalpoints, 4)
    
    # c1 confirms that each of the basis vectors for which we are performing the 4-nn search 
    # is contained in the set points being searched through. c2 confirms that, for each atom
    # in the unit cell, we find three nearest neighbors.
    c1, c2 = true, true
    for i in 1:natoms
        c2 = isapprox(dist[1, i], dist[2, i]) && isapprox(dist[2, i], dist[3, i])
    end

    return c1 && c2
end

#generateLatticeSanityCheck()
#@test nearestNeighborsGrapheneSanityCheck()
@test crystalBravaisInit("graphene")
end