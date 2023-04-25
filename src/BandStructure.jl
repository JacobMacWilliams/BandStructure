module BandStructure
using FromFile
using Plots
using LinearAlgebra
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

function getGraphenePath()
  k1 = zeros(2)
  k1[1] = 1
  k1[2] = sqrt(3)
  k1 = 2 * pi / sqrt(3) * k1
  k2 = zeros(2)
  k2[1] = 1
  k2[2] = -sqrt(3)
  k2 = 2 * pi / sqrt(3) * k2

 # CONSTRUCTING HIGH SYMMETRY PATH

 # from bottom face to origin
  start = k2 / 2
  stop = zeros(2)
  steps = 300
  delta = (stop - start) / steps
  path1 = [start + i*delta for i in 1:steps]

  # from origin to corner
  start = stop
  rot90 = [0 -1; 1 0]
  stop = 1 / sqrt(3) * rot90 * k2
  steps = 300
  delta = (stop - start) / steps
  path2 = [start + i * delta for i in 1:steps]

  # from corner to corner
  start = stop
  stop = start + 1 / sqrt(3) * norm(k2) * [0; -1]
  steps = 300
  delta = (stop - start) / steps
  path3 = [start + i * delta for i in 1:steps]

  # from corner back to origin
  start = stop
  stop = zeros(2)
  steps = 300
  delta = (stop - start) / steps
  path4 = [start + i * delta for i in 1:steps]

  # from origin to top face
  start = stop
  stop = k1 / 2
  steps = 300
  delta = (stop - start) / steps
  path5 = [start + i * delta for i in 1:steps]

  points = cat(path1, path2, path3, path4, path5, dims = 1)
  return points
end

function plotGraphenePath()
  points = getGraphenePath()
  x = []
  y = []
  for point in points
    push!(x, point[1])
    push!(y, point[2])
  end
  plot = scatter(x, y)
  savefig(plot, "path.png")
end

function calculateGrapheneBandStructure(crystal::Crystal)
  if getName(crystal) != "graphene"
    return
  end

  nnarrays = getPrimitiveNearestNeighbors(crystal)
  points = getGraphenePath()
  
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
  bravaisconf = joinpath(@__DIR__, "..", "conf", "bravais.ini")
  if !ispath(bravaisconf)
      bravaisconf = joinpath(@__DIR__, "..", "conf", "bravais.default.ini")
  end

  crystalconf = joinpath(@__DIR__, "..", "conf", "crystal.ini")
  if !ispath(crystalconf)
      crystalconf = joinpath(@__DIR__, "..", "conf", "crystal.default.ini")
  end

  crystal = Crystal("graphene", crystalconf, bravaisconf)
  plotGrapheneBandStructure(crystal)
end

main()
end