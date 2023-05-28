module ElectronLattice
using FromFile
@from "CrystalLattice.jl" using CrystalLattice
@from "CrystalLattice.jl" import CrystalLattice.getName, CrystalLattice.getVecs,
                                 CrystalLattice.getSize
export ElectronCrystal, getName, getCrystal, getStatesPerSite, getElectronNumber

# The true states per site number is 2 * statespersite corresponding
# to the spin degree of freedom.
struct ElectronCrystal
    name::String
    crystal::Crystal
    statespersite::Int32
    population::Int64
end

function ElectronCrystal(name::String, crystal::String, statespersite::Int, population::Int, crystalfile::String, bravaisfile::String)
    crystal = Crystal(crystal, crystalfile, bravaisfile)
    ElectronCrystal(name, crystal, statespersite, population)
end

function getName(ecrystal::ElectronCrystal)
    return ecrystal.name
end

function getCrystal(ecrystal::ElectronCrystal)
    return ecrystal.crystal
end

function getStatesPerSite(ecrystal::ElectronCrystal)
    return ecrystal.statespersite
end

function getElectronNumber(ecrystal::ElectronCrystal)
    return ecrystal.population
end

function getSize(ecrystal::ElectronCrystal)
    crystal = getCrystal(ecrystal)
    return getSize(crystal)
end

function getVecs(ecrystal::ElectronCrystal)
    crystal = getCrystal(ecrystal)
    basisvecs, bravaisvecs = getVecs(crystal)
    return basisvecs, bravaisvecs
end

end