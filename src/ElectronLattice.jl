module ElectronLattice
using FromFile
@from "CrystalLattice.jl" using CrystalLattice
@from "CrystalLattice.jl" import CrystalLattice.getName, CrystalLattice.getVecs,
                                 CrystalLattice.getSize
export ElectronCrystal,
       getName, 
       getCrystal, 
       getVecs, 
       getStatesPerSite,
       getSize,
       getElectronNumber, 
       getHoppingMatrixStructure

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

# This function constructs the hopping matrix between the unit cell located at the origin and a neighboring
# unit cell assuming a tight binding model (i.e. only nearest neighbor interactions).
function getHoppingMatrixStructure(crystal::ElectronCrystal, cell::Int, nnarrays::Matrix{Bool})
  
    basisvecs, _ = getVecs(crystal)
    basis = size(basisvecs, 2)
    n = getStatesPerSite(crystal)
    hopping = ones(basis*n, basis*n)
    ijkl = Iterators.product(1:n, 1:n, 1:basis, 1:basis)
    
    # i labels the spin state "to which the electron is transitioning"
    # j labels the spin state "from which the electron is transitioning"
    # k labels the crystal lattice site to which "the electron is hopping"
    # l labels the crystal lattice site from which "the electron is hopping"
    for (i, j, k, l) in ijkl
        tosite = n * (k - 1)
        fromsite = n * (l - 1)
        tostate = tosite + i
        fromstate = fromsite + j
        isnearest = nnarrays[cell*basis + k, l]

        # Hopping term
        if (isnearest)
            hopping[tostate, fromstate] *= 5 # Mean field correction
            if (i == j)
                hopping[tostate, fromstate] *= 3 # hopping term
            end
        end

        # Interaction correction
        if ((i == j) && (k == l) && (cell == 0))
            hopping[tostate, fromstate] *= 7
        end
    end
    return hopping
  end

end