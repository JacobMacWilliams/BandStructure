module BandStructureBeta
using FromFile
using LinearAlgebra
using Plots
@from "ElectronLattice.jl" using ElectronLattice
@from "CrystalLattice.jl" using CrystalLattice
@from "MeanField.jl" using MeanField
@from "Corrections.jl" using Corrections
@from "ChemicalPotential.jl" using ChemicalPotential
@from "../utils/Distributions.jl" using Distributions
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
    beta = 0.01
    fill = 0.499999999

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
        eigvalues = eigvals(bloch)
        push!(es, eigvalues)
    end
    
    plotbandstructure(es)
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

function plotbandstructure(energies)
    nbands = length(energies[1])
    steps = [i for i in 1:length(energies)]

    plt = plot()
    for i in 1:nbands
        band = [e[i] for e in energies]
        plot!(plt, steps, band)
    end

    savepath = joinpath("plots", "hubbardmodel.png")
    savefig(plt, savepath)
end

function meanfielditeration(ecrystal, latticepoints, nnidxs, nnlabels, mu, beta, fillfactor)
      
    mf = FieldCorrelator(ecrystal)
    pvmap = Dict([(1, 0.0), (2, (0.15, 0.75)), (3, (0.15, 0.25))])
    initmfvalues!(mf, pvmap)
    correlator = getcorrelator(mf)
    diffs = [0.0, 10.0]

    iteration = 1
    while abs(diffs[end] - diffs[end - 1]) > 0.00001 && length(diffs) < 10000
        diff, correlator, energies = meanfielditerationstep(ecrystal, correlator, latticepoints, nnidxs, nnlabels, mu, beta)
        mu = findchempot(2, fillfactor, beta, energies)
        push!(diffs, diff)
        println("Iteration " * string(iteration) * ":" * " difference = " * string(diff))
        iteration += 1
    end

    return correlator
end
function meanfielditerationstep(ecrystal, correlator, latticepoints, nnidxs, nnlabels, mu, beta)
    
    atomspercell = size(getVecs(ecrystal)[1], 2)

    hoppingmatrices = []
    for i in 1:length(eachcol(latticepoints))
        hopping = gethoppingmatrix(ecrystal, correlator, i, nnidxs, nnlabels)
        push!(hoppingmatrices, hopping)
    end

    # BEGIN MEAN FIELD ITERATION
    nextcorrelator = zeros(Complex{Float64}, size(correlator))
    kpoints = discretehexagon(0.05, "constantstep")
    energies = []
    for k in kpoints
       bloch = getblochmatrix(atomspercell, hoppingmatrices, k, latticepoints)
       blochfactors = eigen(bloch)
       eks = blochfactors.values
       push!(energies, eks)

       fromdiagonalbase = blochfactors.vectors
       todiagonalbase = adjoint(fromdiagonalbase) # Unnecessary transformation
       correlatorupdatestep!(ecrystal, k, correlator, nextcorrelator, eks, mu, beta, todiagonalbase)
    end
    energies = hcat(energies...)

    diff = sum(broadcast(abs, correlator - nextcorrelator))
    return diff, nextcorrelator, energies
end

function correlatorupdatestep!(ecrystal::ElectronCrystal, k::Vector{Float64}, correlator, newcorrelator, eigenvalues::Vector{Float64}, mu::Float64, beta::Float64, diagonaltrafo)
    crystal = getCrystal(ecrystal)
    latticepoints, crystalpoints = getPoints(crystal)
    basisvecs, _ = getVecs(ecrystal)
    
    atomspercell = size(basisvecs, 1)
    latticesize = size(latticepoints, 2)
    latticedof = (atomspercell, latticesize)
    sitedof = (2, atomspercell)

    for i in CartesianIndices(correlator)
        (s2, b2, r2, s1, b1) = Tuple(i)
        if correlator[s2, b2, r2, s1, b1] == 0.0
            continue
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
            newcorrelator[s2, b2, r2, s1, b1] += exptopoint * expfrompoint * n * diagonaltrafo[i, fromstateidx] * adjoint(diagonaltrafo)[tostateidx, i]
        end
    end  
end

function gethoppingmatrix(ecrystal::ElectronCrystal, correlator, siteidx::Int, nnidxs, nnlabels)
    basisvecs, _ = getVecs(ecrystal)
    atomspercell = size(basisvecs, 2)

    U = getLocalInteraction(ecrystal)
    V = getNonLocalInteraction(ecrystal)
    
    hopping = zeros(ComplexF64, 2, atomspercell, 2, atomspercell)

    for i in CartesianIndices(hopping)
        (s2, b2, s1, b1) = Tuple(i)
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
        isatomidx(x::Tuple) = (x == atomidx)
        idx =  findnext(isatomidx, nnidxs[b1][1:end], 1)

        if idx === nothing
            hopping[s2, b2, s1, b1] = 0
            continue
        else
            neighbororder = nnlabels[b1][idx]
            if neighbororder <= length(V)
                coloumb = V[neighbororder]
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

function getblochmatrix(atomspercell::Int, hopmats, k::Vector{Float64}, latticepoints::Matrix{Float64})
  
    bloch = zeros(2*atomspercell, 2*atomspercell)
    for (i, point) in enumerate(eachcol(latticepoints))
      kdotR = dot(k, point)
      expkdotR = exp(im*kdotR)
      hopping = hopmats[i]
      bloch += expkdotR * hopping
    end
    
    return bloch
end



bravaisconf = joinpath("conf", "bravais.default.toml")
crystalconf = joinpath("conf", "crystal.default.toml")
electronconf = joinpath("conf", "elattice.default.toml")
configpaths = [electronconf, crystalconf, bravaisconf]
main("hubbard model", configpaths)
end