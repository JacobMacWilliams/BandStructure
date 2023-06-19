module BandStructure
using FromFile
using Plots
using LinearAlgebra
@from "CrystalLattice.jl" using CrystalLattice
@from "../utils/Reciprocal.jl" using Reciprocal
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

function getBlochMatrix(k::Vector{Float64}, crystal::Crystal, nnarrays::Matrix{Bool})
  
  basisvecs, bravaisvecs = getVecs(crystal)
  dim = size(basisvecs, 1)
  basis = size(bravaisvecs, 2)
  bravaisvecs = cat(zeros(dim), bravaisvecs, -bravaisvecs, dims=2)
  bloch = zeros(basis, basis)
  
  for (i, bravaisvec) in enumerate(eachcol(bravaisvecs))
    kdotR = dot(k, bravaisvec)
    expkdotR = exp(im*kdotR)
    hopping = getHoppingMatrix(crystal, i - 1, nnarrays)
    bloch += expkdotR * hopping
  end
  
  return eigvals(bloch)
end

function calculateGrapheneBandStructure(crystal::Crystal)
  if getName(crystal) != "graphene"
    return
  end

  nnarrays = getPrimitiveNearestNeighbors(crystal)
  points = getgraphenepath()
  
  x = []
  z1 = []
  z2 = []
  for (i, point) in enumerate(points)
    push!(x, i)
    eigvalues = getBlochMatrix(point, crystal, nnarrays)
    push!(z1, eigvalues[1])
    push!(z2, eigvalues[2])
  end

  return x, z1, z2

end

function plotGrapheneBandStructure(crystal::Crystal)
  if getName(crystal) != "graphene"
    return
  end
  x, z1, z2 = calculateGrapheneBandStructure(crystal)
  plot(x, z1)
  plt = plot!(x, z2)
  pltpath = joinpath("plots", "graphenebands.png")
  savefig(plt, pltpath)
end

function main()
  bravaisconf = joinpath(@__DIR__, "..", "conf", "bravais.toml")
  if !ispath(bravaisconf)
      bravaisconf = joinpath(@__DIR__, "..", "conf", "bravais.default.toml")
  end

  crystalconf = joinpath(@__DIR__, "..", "conf", "crystal.toml")
  if !ispath(crystalconf)
      crystalconf = joinpath(@__DIR__, "..", "conf", "crystal.default.toml")
  end

  crystal = Crystal("graphene", crystalconf, bravaisconf)
  plotGrapheneBandStructure(crystal)
end

main()
end