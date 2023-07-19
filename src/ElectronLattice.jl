module ElectronLattice
using FromFile
@from "../utils/TOMLIO.jl" using TOMLIO
@from "CrystalLattice.jl" using CrystalLattice
@from "CrystalLattice.jl" import CrystalLattice.getName, 
                                 CrystalLattice.getVecs
export ElectronCrystal,
       getName, 
       getCrystal, 
       getVecs, 
       getStatesPerSite,
       getLocalInteraction,
       getNonLocalInteraction,
       getHoppingMatrixStructure

const CLASSKEY = "electronlattice"
struct ElectronCrystal
    name::String
    crystal::Crystal
    U::Float64
    V::Vector{Float64}
    statespersite::Int32
end

function ElectronCrystal(name::String, configfile::String, crystalfile::String, bravaisfile::String)
    config = parseconfig(configfile, CLASSKEY, name)
    crystalname = get(config, "crystal", nothing)
    U = get(config, "U", nothing)
    V = get(config, "V", nothing)
    statespersite = get(config, "statespersite", nothing)
    crystal = Crystal(crystalname, crystalfile, bravaisfile)
    ElectronCrystal(name, crystal, U, V, statespersite)
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

function getVecs(ecrystal::ElectronCrystal)
    crystal = getCrystal(ecrystal)
    basisvecs, bravaisvecs = getVecs(crystal)
    return basisvecs, bravaisvecs
end

function getLocalInteraction(ecrystal::ElectronCrystal)
    return ecrystal.U
end

function getNonLocalInteraction(ecrystal::ElectronCrystal)
    return ecrystal.V
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