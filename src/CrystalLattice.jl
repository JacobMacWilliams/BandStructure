module CrystalLattice
using FromFile
# include(joinpath("..", "utils", "CrystalIO.jl"))
@from "BravaisLattice.jl" using BravaisLattice: Bravais
@from "BravaisLattice.jl" import BravaisLattice.getName
@from "BravaisLattice.jl" import BravaisLattice.getVecs
@from "../utils/CrystalIO.jl" using CrystalIO: parseCrystalConf as getCrystal
export Crystal, getName, getSize, getLattice, getVecs

struct Crystal
  name::String
  size::Integer
  basis::Matrix{Float64}
  lattice::Bravais
end

function Crystal(name::String, crystalfile::String, bravaisfile::String)
  N, basisvecs, latticename = getCrystal(crystalfile, name)
  lattice = Bravais(latticename, bravaisfile)
  Crystal(name, N, basisvecs, lattice)
end

function getName(crystal::Crystal)
  return crystal.name
end

function getSize(crystal::Crystal)
  return crystal.size
end

function getLattice(crystal::Crystal)
  return crystal.lattice
end

function getVecs(crystal::Crystal)
  lattice = getLattice(crystal)
  bravaisvecs = getVecs(lattice)
  basisvecs = crystal.basis
  return basisvecs, bravaisvecs
end

end

