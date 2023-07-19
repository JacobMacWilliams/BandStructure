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
    kpoints = discretehexagon(0.05, "constantstep")
    
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
    fill = 0.499999999

    correlator = meanfielditeration(ecrystal, latticepoints, kpoints, idxs_ij, nnlabels, mu, beta, fill)

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
        imeigvalues = broadcast(imag, eigvalues)
        if any(imeigvalues .> 1e-12)
            error("Non-negligible imaginary energy value for hermitian operator was found.")
        end
        eigvalues = broadcast(real, eigvalues)
        push!(es, eigvalues)
    end
    
    plotbandstructure(es)
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

function meanfielditeration(ecrystal, latticepoints, kpoints, nnidxs, nnlabels, mu, beta, fillfactor)
      
    mf = FieldCorrelator(ecrystal)
    pvmap = Dict([(1, -1.0), (2, 0.0), (3, 1.0)])
    initmfvalues!(mf, pvmap)
    correlator = getcorrelator(mf)
    diffs = [0.0, 10.0]

    iteration = 1
    while abs(diffs[end] - diffs[end - 1]) > 0.00001 && length(diffs) < 10000
        energies, trafos = getEigenEnergies(ecrystal, correlator, latticepoints, kpoints, nnidxs, nnlabels)
        earray = hcat(energies...)
        mu = findchempot(2, fillfactor, beta, earray)

        newcorrelator = meanfielditerationstep(ecrystal, correlator, kpoints, energies, trafos, mu, beta)
        diff = sum(broadcast(abs, correlator - newcorrelator))
        push!(diffs, diff)
        println("Iteration " * string(iteration) * ":" * " difference = " * string(diff))
        
        correlator = newcorrelator
        iteration += 1
    end

    return correlator
end

function getEigenEnergies(ecrystal, correlator, latticepoints, kpoints, nnidxs, nnlabels)
    
    atomspercell = size(getVecs(ecrystal)[1], 2)
    hoppingmatrices = []
    for i in 1:length(eachcol(latticepoints))
        hopping = gethoppingmatrix(ecrystal, correlator, i, nnidxs, nnlabels)
        push!(hoppingmatrices, hopping)
    end

    energies = []
    fromdiagonaltrafos = []
    
    for k in kpoints
        bloch = getblochmatrix(atomspercell, hoppingmatrices, k, latticepoints)
        blochfactors = eigen(bloch)
        eks = blochfactors.values
        
        imeks = broadcast(imag, eks)
        if any(imeks .> 1e-12)
            error("Non-negligible imaginary energy value for hermitian operator was found.")
        end

        eks = broadcast(real, eks)
        push!(energies, eks)

        fromdiagonalbase = blochfactors.vectors
        push!(fromdiagonaltrafos, fromdiagonalbase)
    end

    return energies, fromdiagonaltrafos
end

function meanfielditerationstep(ecrystal, correlator, kpoints::Vector{Any}, energies::Vector{Any}, trafos, mu, beta)
    crystal = getCrystal(ecrystal)
    latticepoints, crystalpoints = getPoints(crystal)
    basisvecs, _ = getVecs(ecrystal)
    
    atomspercell = size(basisvecs, 1)
    latticesize = size(latticepoints, 2)
    latticedof = (atomspercell, latticesize)
    sitedof = (2, atomspercell)
    nk = length(kpoints)

    newcorrelator = zeros(Complex{Float64}, size(correlator))

    for i in CartesianIndices(correlator)
        (s2, b2, r2, s1, b1) = Tuple(i)
        if correlator[s2, b2, r2, s1, b1] == -1.0
            newcorrelator[s2, b2, r2, s1, b1] = -1.0
            continue
        end

        tositeidx = (b2, r2)
        fromsiteidx = (b1, 1)
        
        toidx, fromidx = embedinmatrix(tositeidx, fromsiteidx, latticedof, latticedof)
        
        topoint = crystalpoints[:, toidx]
        frompoint = crystalpoints[:, fromidx]

        for (j, k) in enumerate(kpoints)
            es = energies[j]
            for (l, e) in enumerate(es)
                diagonaltrafo = trafos[j]
                exptopoint = exp(im * dot(k, topoint))
                expfrompoint = exp(- im * dot(k, frompoint))
                n = fermidistribution(e, mu, beta)
                tostateidx = (s2, b2)
                fromstateidx = (s1, b1)

                tostateidx, fromstateidx = embedinmatrix(tostateidx, fromstateidx, sitedof, sitedof)
                newcorrelator[s2, b2, r2, s1, b1] += exptopoint * expfrompoint * n * diagonaltrafo[fromstateidx, l] * conj(diagonaltrafo[tostateidx, l]) / nk
            end
        end
    end  
return newcorrelator
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
                for s in 1:2
                    if (s != s1)
                        hopping[s2, b2, s1, b1] += U *  correlator[s, b1, siteidx, s, b1]
                    end 
                end
            end

            # COLOUMB POTENTIAL CORRECTIONS
            #=
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
            =#
        end

        # HANDLING HOPPING AND HOPPING CORRECTIONS
        atomidx = (b2, siteidx)
        isatomidx(x::Tuple) = (x == atomidx)
        idx =  findnext(isatomidx, nnidxs[b1][1:end], 1)

        if idx === nothing
            continue
        end

        neighbororder = nnlabels[b1][idx]
        #=
        if neighbororder <= length(V)
            coloumb = V[neighbororder]
            hopping[s2, b2, s1, b1] += - coloumb * conj(correlator[s2, b2, siteidx, s1, b1])
        end
        =#

        # HANDLING NEAREST NEIGHBOR HOPPING
        if (neighbororder == 1) && (s2 == s1)
            hopping[s2, b2, s1, b1] += - 1
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