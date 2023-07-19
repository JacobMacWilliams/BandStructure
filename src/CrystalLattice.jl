module CrystalLattice
using FromFile
using LinearAlgebra: norm
using NearestNeighbors
@from "BravaisLattice.jl" using BravaisLattice: Bravais, getlatticepoints
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

function getPoints(crystal::Crystal)
  lattice = getLattice(crystal)
  latticepoints = getlatticepoints(lattice)
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
  idxs, dists = knn(tree, basisvecs, k, true)

  # Remove the entries for which zero distance between two points
  # is found. This corresponds to a distance measurement between
  # a point and itself.
  idxs = [idx[2:end] for idx in idxs]
  nnlabels = [labelneighbors(d[2:end]) for d in dists]
  
  return idxs, nnlabels
end

function labelneighbors(dist)

  len = length(dist)
  label = 1
  neighborlabels = [label]

  dlast = dist[1]
  for i in 2:len
    dnext = dist[i]
    !isapprox(dnext, dlast) ? label+=1 : nothing
    push!(neighborlabels, label)
    dlast = dnext
  end

  return neighborlabels
end

end

