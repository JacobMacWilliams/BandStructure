module BandStructureTest
using FromFile
using Test
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/BandStructure.jl" using BandStructure

function getHoppingMatrixTest()
    bravaisconf = joinpath("conf", "bravais.default.toml")
    crystalconf = joinpath("conf", "crystal.default.toml")
    crystal = Crystal("graphene", crystalconf, bravaisconf)
    nnarrays = getPrimitiveNearestNeighbors(crystal)
    hopping = getHoppingMatrix(crystal, 1, nnarrays)
    expected = zeros(2, 2)
    expected[1, 2] = -1 # only hopping from B to A in primitive cell a1
    return (expected == hopping)
end

@test getHoppingMatrixTest()

end