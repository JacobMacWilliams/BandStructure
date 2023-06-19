module ChemicalPotentialTest
using FromFile, Test
using LinearAlgebra
using Plots
@from "../src/ElectronLattice.jl" using ElectronLattice
@from "../src/ChemicalPotential.jl" using ChemicalPotential
@from "../utils/Rotations.jl" using Rotations
@from "../utils/Reciprocal.jl" using Reciprocal

function highband(k, r)      
    total = 3.0
    for j in 0:5
        a = rotate2d(r, j*pi / 3)
        term = cos(dot(k, a))
        total += term
    end

    if isapprox(total, 0.0, atol=1e-10, rtol=0)
        total = 0
    end
    energy = sqrt(total)
    return energy
end

function highbandTest()
    crystalconf = joinpath("conf", "crystal.default.toml")
    bravaisconf = joinpath("conf", "bravais.default.toml")
    electronconf = joinpath("conf", "elattice.default.toml")
    ecrystal = ElectronCrystal("hubbard model", electronconf, crystalconf, bravaisconf)
    _, bravais = getVecs(ecrystal)

    points = getgraphenepath()
    println(bravais[:, 1])
    energies = collect(Iterators.map(p -> highband(p, bravais[:,1]), points))
    return energies
end
        
function plothighbandTest()
    energies = highbandTest()
    plot([i for i in 1:length(energies)], energies)
    savepath = joinpath("test", "plots", "analytichighband.png")
    savefig(savepath)
end

function chemicalpotentialtest()
    #TODO: WE ONLY NEED THE BRAVAIS VECTORS THEREFORE WE CAN JUST INSTANTIATE THE BRAVAIS LATTICE.
    crystalconf = joinpath("conf", "crystal.default.toml")
    bravaisconf = joinpath("conf", "bravais.default.toml")
    electronconf = joinpath("conf", "elattice.default.toml")
    ecrystal = ElectronCrystal("hubbard model", electronconf, crystalconf, bravaisconf)

    # Calculate energy of each k point in the BZ for both bands.
    _, bravais = getVecs(ecrystal)
    points = discretehexagon(0.01, "constantstep")
    energyhigh = collect(Iterators.map(x -> highband(x, bravais[:, 1]), points))
    energylow = collect(Iterators.map(x -> -highband(x, bravais[:, 1]), points))
    energies = []
    push!(energies, energyhigh...)
    push!(energies, energylow...)

    # find chemical potential
    mu = findchempot(2, 0.5, 100.0, energies)
    return isapprox(mu, 0.0)
end
#plothighbandTest()
#plotpath()
@test chemicalpotentialtest()
end