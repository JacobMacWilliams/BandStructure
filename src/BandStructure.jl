module BandStructure
using FromFile
@from "CrystalLattice.jl" using CrystalLattice
export getHoppingMatrix
# This function assumes a tight binding model with only nearest neighbor hopping.
# Furthermore we assume that the nearest neighbors are of an equal displacement 
# from the point as every other point.
function getHoppingMatrix(crystal::Crystal, cell::Int, nnarrays::Matrix{Bool})
  basisvecs, _ = getVecs(crystal)
  basis = size(basisvecs, 2)
  hopping = zeros(basis, basis)
  ij = Iterators.product(1:basis, 1:basis)
  for (i, j) in ij
    hopping[i, j] = nnarrays[cell*basis + i, j] ? -1 : 0
  end
  return hopping
end
end