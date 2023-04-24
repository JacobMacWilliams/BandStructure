module BandStructureTest
using FromFile
using Test
@from "../utils/CrystalIO.jl" using CrystalIO
@from "../src/CrystalLattice.jl" using CrystalLattice
@from "../src/BandStructure.jl" using BandStructure

bravaisconf = joinpath("conf", "bravais.ini")
if !ispath(bravaisconf)
    bravaisconf = joinpath("conf", "bravais.default.ini")
end

crystalconf = joinpath("conf", "crystal.ini")
if !ispath(crystalconf)
    crystalconf = joinpath("conf", "bravais.default.ini")
end

open(bravaisconf, truncate=true)
open(crystalconf, truncate=true)
configureGraphene(crystalconf, bravaisconf)

function getHoppingMatrixTest()
    crystal = Crystal("graphene", crystalconf, bravaisconf)
    nnarrays = getPrimitiveNearestNeighbors(crystal)
    println(nnarrays)
    hopping = getHoppingMatrix(crystal, 1, nnarrays)
    expected = zeros(2, 2)
    expected[1, 2] = -1 # only hopping from B to A in primitive cell a1
    println(hopping)
    println(expected)
    return (expected == hopping)
end

@test getHoppingMatrixTest()

end