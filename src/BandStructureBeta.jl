module BandStructureBeta
using FromFile
using LinearAlgebra
using Plots
@from "ElectronLattice.jl" using ElectronLattice
@from "CrystalLattice.jl" using CrystalLattice
@from "MeanField.jl" using MeanField
@from "Corrections.jl" using Corrections
@from "ChemicalPotential.jl" using ChemicalPotential
@from "../utils/Embeddings.jl" using Embeddings
@from "../utils/Reciprocal.jl" using Reciprocal

function main(modelname::String, configfiles::Vector{String})
    
    ecrystal = ElectronCrystal(modelname, configfiles...)
    
    #=
    Generate the crystal points up to the n-th order whereby n is given by how far out from a point
    the coulomb interaction is specfied. If it is specified for just nearest neighbors the order is 
    one, if it is specified for nearest and next nearest neighbors the order is 2. It is assumes this
    guarantees that the set of all neighbors for which there is a non-zero coloumb interaction is
    present in the set of all points generated to the nth-order. This is not generally true, and is
    a potential point of failure.
    =#
    coulomb = getNonLocalInteraction(ecrystal)
    
    crystal = getCrystal(ecrystal)
    latticepoints, crystalpoints = getPoints(crystal)
    
    #TODO: Setting k to 4 is only appropriate for the hubbard model
    idxs_i, nnlabels = getNearestNeighborsTest(crystal, crystalpoints, 4)

    #=
    With the column indices of the crystalpoints which constitute the nearest
    neighbors of each of the atoms in the unit cell positioned at the origin,
    we can retreive the corresponding bravais vector and basis vector indices
    with our knowledge of the number of atoms in each unit cell (in other words
    knowing how many points are generated for each point of the bravais lattice)
    =#
    atomspercell = length(idxs_i)
    cellnumber = size(latticepoints, 2)
    dof = [atomspercell, cellnumber]
    
    idxs_ij = []
    for idxs in idxs_i
        ij = collect(Iterators.map(idx -> extractfromlist(idx, dof), idxs))
        push!(idxs_ij, ij)
    end

  
    
    mu = 0.0
    beta = 100.0
    fill = 0.5
    correlator = meanfielditeration(ecrystal, latticepoints, idxs_ij, nnlabels, mu, beta, fill)

    hoppingmatrices = []
    for i in 1:length(eachcol(latticepoints))
        hopping = gethoppingmatrix(ecrystal, correlator, i, idxs_ij, nnlabels)
        push!(hoppingmatrices, hopping)
    end

    points = getgraphenepath()
    es = []
    for k in points
        bloch = getblochmatrix(atomspercell, hoppingmatrices, k, latticepoints)
        eigvals = eigvals(bloch)
        push!(es, eigvals)
    end
    
    plt = scatter(1, 0.0) ## DUMMY POINT
    for (i, e) in enumerate(es)
        scatter!(i, e[1])
        scatter!(i, e[2])
    end

    savepath = joinpath("plots", "hubbardmodel.png")
    savefig(plt, savepath)
    # TEST CODE
    #=
    c = true
    for (basisoriginidx, idxs) in enumerate(eachcol(idxs_ij))
        for (n, idx) in enumerate(idxs)
            bravaisidx = Int(idx[2])
            baseidx = Int(idx[1])
            println(basisoriginidx)
            println(n)
            latticeidx = idxs_i[n, basisoriginidx]
            bravaisvec = latticepoints[:, bravaisidx]
            basevec = basis[:, baseidx]
            point = bravaisvec + basevec

            c = c && (crystalpoints[:, latticeidx] == point)
            println(c)
        end
    end
    println(c)
    =#
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

function meanfielditeration(ecrystal, latticepoints, nnidxs, nnlabels, mu, beta, fillfactor)
      
    mf = FieldCorrelator(ecrystal, hubbardslices)
    pvmap = Dict([(1, 0.0), (2, (0.15, 0.75)), (3, (0.15, 0.25))])
    initmfvalues!(mf, pvmap)
    correlator = getfieldcorrelator(mf)
    diffs = [10.0]

    while (diff[end] - diff[end - 1]) > 0.00001 && length(diffs) < 10000
        diff, correlator, energies = meanfielditerationstep(ecrystal, correlator, latticepoints, nnidxs, nnlabels, mu, beta)
        mu = chempot(2, fillfactor, beta, mu, energies)
        push!(diffs, diff)
    end

    return correlator
end
function meanfielditerationstep(ecrystal, correlator, latticepoints, nnidxs, nnlabels, mu, beta)
    
    hoppingmatrices = []
    for i in 1:length(eachcol(latticepoints))
        hopping = gethoppingmatrix(ecrystal, correlator, i, nnidxs, nnlabels)
        push!(hoppingmatrices, hopping)
    end

    # BEGIN MEAN FIELD ITERATION
    nextcorrelator = zeros(sizes(correlator))
    kpoints = discretehexagontest(0.05, "constantstep")
    energies = []
    for k in eachcol(kpoints)
       bloch = getblochmatrix(atomspercell, hoppingmatrices, k, latticepoints)
       blochfactors = eigen(bloch)
       eks = blochfactors.values
       push!(energies, eks)

       fromdiagonalbase = blochfactors.vectors
       todiagonalbase = adjoint(fromdiagonalbase) # Unnecessary transformation
       correlatorupdatestep!(ecrystal, k, correlator, nextcorrelator, eks, mu, beta, todiagonalbase)
    end
    energies = hcat(energies...)
    diff = norm(correlator - nextcorrelator)
    return diff, nextcorrelator, energies
end

function gethoppingmatrix(ecrystal::ElectronCrystal, correlator, siteidx::Int, nnidxs::Matrix{Int}, nnlabels::Matrix{Int})
    basisvecs, _ = getVecs(ecrystal)
    atomspercell = size(basisvecs, 2)

    U = getLocalInteraction(ecrystal)
    V = getNonLocalInteraction(ecrystal)
    
    hopping = zeros(2, atomspercell, 2, atomspercell)

    for (s2, b2, s1, b1) in CartesianIndex(hopping)
        # Handling spin flip and potential correction        
        if (siteidx == 1)

            # SPINFLIP
            if (b2 == b1) && (s2 != s1)
                hopping[s2, b2, s1, b1] += -U * conj(correlator[s2, b2, siteidx, s1, b1])
            end

            # POTENTIAL CORRECTIONS
            # PAIRING POTENTIAL CORRECTIONS
            if (b2 == b1) && (s2 == s1)
                hopping[s2, b2, s1, b1] += U * sum([s != s1 ? correlator[s, b1, siteidx, s, b1] : 0 for s in 1:2])
            end

            # COLOUMB POTENTIAL CORRECTIONS
            if (b2 == b1) && (s2 == s1)
                nnidx = nnidxs[b1]
                for (i, idx) in enumerate(nnidx)
                    neighbororder = nnlabels[b1][i]
                    if neighbororder > length(V)
                        break
                    end

                    coulomb = V[neighbororder]
                    corsum = sum([correlator[s, idx[1], siteidx, s, idx[1]] for s in 1:2])
                    hopping[s2, b2, s1, b1] += coulomb * corsum
                end
            end
        end

        # HANDLING HOPPING AND HOPPING CORRECTIONS
        atomidx = (b2, siteidx)
        idx =  findnext(nnidxs[b1][1:end], atomidx)
        if idx === nothing
            hopping[s2, b2, s1, b1] = 0
            continue
        else
            neighbororder = nnlabels[idx]
            if neighbororder <= length(V)
                coloumb = V[neighbororder]
                if correlator[s2, b2, siteidx, s1, b1] == 0.0
                    error("Sanity check failed.")
                end
                hopping[s2, b2, s1, b1] += - coloumb * conj(correlator[s2, b2, siteidx, s1, b1])
            end

            # HANDLING NEAREST NEIGHBOR HOPPING
            if (neighbororder == 1) && (s2 == s1)
                hopping[s2, b2, s1, b1] += - 1
            end
        end
    end

    hopping = reshape(hopping, (2*atomspercell, 2*atomspercell))
    return hopping
end

function getblochmatrix(atomspercell::Int, hopmats::Matrix{Float64}, k::Vector{Float64}, latticepoints::Matrix{Float64})
  
    bloch = zeros(2*atomspercell, 2*atomspercell)
    for (i, point) in enumerate(eachcol(latticepoints))
      kdotR = dot(k, point)
      expkdotR = exp(im*kdotR)
      hopping = hopmats[i]
      bloch += expkdotR * hopping
    end
    
    return bloch
end

function correlatorupdatestep!(ecrystal::ElectronCrystal, k::Vector{Float64}, correlator::Matrix{Float64}, newcorrelator::Matrix{Float64}, eigenvalues::Vector{Float64}, mu::Float64, beta::Float64, diagonaltrafo::Matrix{Float64})
    crystal = getCrystal(ecrystal)
    latticepoints, crystalpoints = getPoints(crystal)
    basisvecs, _ = getVecs(ecrystal)
    
    atomspercell = size(basisvecs, 1)
    latticesize = size(latticepoints, 2)
    latticedof = (atomspercell, latticesize)
    sitedof = (2, atomspersite)

    for (s2, b2, r2, s1, b1) in CartesianIndices(correlator)

        if correlator[s2, b2, r2, s1, b1] == 0.0
            continue;
        end

        tositeidx = (b2, r2)
        fromsiteidx = (b1, 1)
        
        toidx, fromidx = embedinmatrix(tositeidx, fromsiteidx, latticedof, latticedof)
        
        topoint = crystalpoints[:, toidx]
        frompoint = crystalpoints[:, fromidx]

        for (i, e) in enumerate(eigenvalues)
            exptopoint = exp(i * dot(k, topoint))
            expfrompoint = exp(- i * dot(k, frompoint))
            n = fermidistribution(e, mu, beta)
            tostateidx = (s2, b2)
            fromstateidx = (s1, b1)

            tostateidx, fromstateidx = embedinmatrix(tostateidx, fromstateidx, sitedof, sitedof)
            newcorrelator[s2, b2, r2, s1, b1] += exptopoint * expfrompoint * n * diagonaltrafo[n, fromstateidx] * adjoint(diagonaltrafo)[tostateidx, n]
        end
    end  
end

#=
function nearestneigborsplay()
    electronconf = joinpath("..", "conf", "elattice.default.toml")
    bravaisconf = joinpath("..", "conf", "bravais.default.toml")
    crystalconf = joinpath("..", "conf", "crystal.default.toml")
    ecrystal = ElectronCrystal(model, electronconf, crystalconf, bravaisconf)
    crystal = getCrystal(ecrystal)

    coulomb = getNonLocalInteraction(ecrystal)
    order = size(coulomb, 1) #TODO: Doesn't work so long as coulomb is integer
    bravais, crystal = getPoints(crystal, order)
    idx, dist = getNearestNeighborsTest(crystal)
    idx = parseknearest(idx, dist)
    
    #=
    EXTRACT APPROPRIATE NEAREST NEIGHBORS FROM LIST

    FROM EACH IDX FROM LIST WITH SIZES SIZE(BRAVAIS), SIZE(BASIS) THIS WILL GIVE ME THE INDEX OF BRAVAIS AND 
    BASIS USED TO GET THE point

    I WILL THEN HAVE THE CORRESPONDING IJ INDICES 

    FOR CONSTRUCTING THE HOPPING MATRIX I CAN JUST SEE IF IJ IS ONE OF THE NEIGHBORS, IF IT IS ADD APPROPRIATE CORRECTIONS
    =#

    

end
=#
bravaisconf = joinpath("conf", "bravais.default.toml")
crystalconf = joinpath("conf", "crystal.default.toml")
electronconf = joinpath("conf", "elattice.default.toml")
configpaths = [electronconf, crystalconf, bravaisconf]
main("hubbard model", configpaths)
end