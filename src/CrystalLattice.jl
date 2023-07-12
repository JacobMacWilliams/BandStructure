module CrystalLattice
using FromFile
using LinearAlgebra: norm
using NearestNeighbors
@from "BravaisLattice.jl" using BravaisLattice: Bravais, getLatticePoints
@from "BravaisLattice.jl" import BravaisLattice.getName
@from "BravaisLattice.jl" import BravaisLattice.getVecs
@from "../utils/TOMLIO.jl" using TOMLIO

export Crystal,
       getName,
       getSize,
       getLattice,
       getVecs,
       getPrimitiveNearestNeighbors,
       getPoints,
       getNearestNeighborsTest

const CLASSKEY = "crystal"

struct Crystal
  name::String
  basis::Matrix{Float64}
  lattice::Bravais
end

function Crystal(name::String, crystalfile::String, bravaisfile::String)
  config = parseconfig(crystalfile, CLASSKEY, name)
  bravaisname = get(config, "lattice", nothing)
  basisvecs = parsematrix(config, "v")
  lattice = Bravais(bravaisname, bravaisfile)
  Crystal(name, basisvecs, lattice)
end

function getName(crystal::Crystal)
  return crystal.name
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

function getCellPoints(basisvecs, bravaisvec)
  cellpoints = hcat(Iterators.map(b -> b + bravaisvec, eachcol(basisvecs))...)
  return cellpoints
end

function getPoints(crystal::Crystal, order)
  lattice = getLattice(crystal)
  latticepoints = getLatticePoints(lattice, order)
  basisvecs, _ = getVecs(crystal)

  crystalpoints = []
  for v in eachcol(latticepoints)
    cellpoints = getCellPoints(basisvecs, v)
    if size(crystalpoints, 1) == 0
      crystalpoints = cellpoints
      continue
    end
    crystalpoints = hcat(crystalpoints, cellpoints)
  end

  return latticepoints, crystalpoints
end
 
function getNearestNeighborsTest(crystal, points, k)
  basisvecs, _ = getVecs(crystal)
  tree = KDTree(points, Euclidean())
  idxs, dist = knn(tree, basisvecs, k, true)

  # Remove the entries for which zero distance between two points
  # is found. This corresponds to a distance measurement between
  # a point and itself.
  idxs = [idx[2:end] for idx in idxs]
  dist = [d[2:end] for d in dist]
  return idxs, dist
end

# This function has a number of parameters simply because it only serves
# as a helper function to getPrimitiveNearestNeighbors, and the paramters shouldn't
# be recalculated if they don't change across multiple calls.
function getNearestNeighbors(latticeiter, currentvec, dim, basis)
  #=
  The nearest neighbor list here is a boolean list containing true at position
  n[j*basis + i] if basisvecs[:, i] + bravaisvecs[:, j] is a nearest neighbor 
  of the currentvector.
  =#
  nnarray = Vector{Bool}(undef, (2*dim + 1)*basis)
  minidx = 0
  min::Float64 = NaN64
  for (idx, (basisvec, primitive)) in enumerate(latticeiter)
      position = primitive + basisvec
      d = norm(position - currentvec)

      if isapprox(d, 0) # case: position is the currentvec
          nnarray[idx] = false
      elseif (isnan(min) || (!isapprox(min, d) && d < min)) # case: position is true new nn
          min = d
          minidx = idx
          nnarray[idx] = true
      elseif isapprox(min, d) # case: position is equally distant from currentvec as last nn
          nnarray[idx] = true
      else # case: position is further from currentvec as last nn
          nnarray[idx] = false
      end
  end

  # minidx stores the index of the first true nearest neighbor since findNearestNeighbor 
  # updates minidx every time a point is found whos distance to the point stored at bidx
  # is strictly smaller then the last smallest distance found. For this reason all true
  # assignments to the nnarray list coming after minidx are correct and reflect true 
  # nearest neighbors while all true assignments coming before this index are incorrect
  # and must be changed to false.
  nnview = @view nnarray[1:minidx - 1]
  fill!(nnview, false)
  return nnarray
end

function getPrimitiveNearestNeighbors(crystal::Crystal)::Matrix{Bool}

  basisvecs, bravaisvecs = getVecs(crystal)
  dim = size(basisvecs, 1)
  basis = size(basisvecs, 2)
  bravaisvecs = cat(zeros(dim), bravaisvecs, -bravaisvecs, dims = 2)
  latticeiter = Iterators.product(eachcol(basisvecs), eachcol(bravaisvecs))
  nnarrays = []
  for basisvec in eachcol(basisvecs)
    nnarray = getNearestNeighbors(latticeiter, basisvec, dim, basis)
    push!(nnarrays, nnarray)
  end
  nnarrays = cat(nnarrays..., dims = 2)

  return nnarrays
end

end

