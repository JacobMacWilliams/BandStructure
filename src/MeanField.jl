module MeanField
using FromFile
@from "ElectronLattice.jl" using ElectronLattice
@from "../utils/Embeddings.jl" using Embeddings
export FieldCorrelator,
       getcorrelator,
       getlocalinteraction,
       getnonlocalinteraction,
       getprimeslicemap,
       initmfstruct!,
       initmfvalues!,
       hubbardprimeslicemapgraphene

struct FieldCorrelator
    localint::Bool
    nonlocalint::Bool
    primeslicemap::Dict # maps primes to values
    correlator::Array
end

function FieldCorrelator(ecrystal::ElectronCrystal, primeslicemap::Dict, localint::Bool)
    FieldCorrelator(ecrystal, primeslicemap, localint, [])
end

function FieldCorrelator(ecrystal::ElectronCrystal, primeslicemap::Dict, localint::Bool, potential::Array)
    
    basisvecs, _ = getVecs(ecrystal)
    nbasis = size(basisvecs, 2)
    nspins = getStatesPerSite(ecrystal)
    
    npot = size(potential)
    nonlocal = (npot[1] != 0)

    totsize = 1
    for d in npot
        totsize *= d
    end

    if !nonlocal
        ncells = 1
    else
        ncells = Int64(totsize / nbasis ^ 2)
    end

    nstates = nspins * nbasis * ncells
    correlator = ones(nstates, nstates)
    FieldCorrelator(localint, nonlocal, primeslicemap, correlator)
end

function getcorrelator(mf::FieldCorrelator)
    return mf.correlator
end

function getlocalinteraction(mf::FieldCorrelator)
    return mf.localint
end

function getnonlocalinteraction(mf::FieldCorrelator)
    return mf.nonlocalint
end

function getprimeslicemap(mf::FieldCorrelator)
    return mf.primeslicemap
end

function initmfstruct!(mf::FieldCorrelator)
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    correlator = getcorrelator(mf)
    primeslicemap = getprimeslicemap(mf)
    for p in primes
        if !haskey(primeslicemap, p)
            break
        end
        for ij in get(primeslicemap, p, nothing)
            correlator[ij[1], ij[2]] *= p
        end
    end
end

function initmfvalues!(mf::FieldCorrelator, primevaluemap::Dict)
    correlator = getcorrelator(mf)
    for (i, p) in enumerate(correlator)
        if !haskey(primevaluemap, p)
            error("The correlator contains an undefined prime code.")
        end
        
        v = primevaluemap[p]
        if typeof(v) == Tuple{Float64, Float64}
            v = randn() * 2 * v[1] ^ 2  + v[2] # gaussian with mean v[1] and std v[2]
        end
        correlator[i] = v
    end
end

function hubbardprimeslicemapgraphene(basisize::Int)
    allstates = Iterators.product(1:2, 1:basisize, 1:2, 1:basisize)
    pspin(i) = (i[1] == i[3]) && (i[2] == i[4]) ? true : false
    apspin(i) = (i[1] != i[3]) && (i[2] == i[4]) ? true : false
    sizes = (2, basisize)
    #map2d(i) = (2 * (i[2] - 1) + i[1], 2 * (i[4] - 1) + i[3])

    pspiniter = Iterators.filter(pspin, allstates)
    pspiniter = Iterators.map(x -> embedinmatrix((x[1], x[2]), (x[3], x[4]), sizes), pspiniter)

    apspiniter = Iterators.filter(apspin, allstates)
    apspiniter = Iterators.map(x -> embedinmatrix((x[1], x[2]), (x[3], x[4]), sizes), apspiniter)

    primeslicemap = Dict([(2, pspiniter), (3, apspiniter)])
    return primeslicemap
end
end