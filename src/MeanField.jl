module MeanField
using FromFile
@from "ElectronLattice.jl" using ElectronLattice
@from "../utils/Embeddings.jl" using Embeddings
@from "CrystalLattice.jl" using CrystalLattice: getLattice
@from "BravaisLattice.jl" using BravaisLattice: getlatticepoints
export FieldCorrelator,
       getcorrelator,
       getlocalinteraction,
       getnonlocalinteraction,
       getprimeslicemap,
       initmfstruct!,
       initmfvalues!,
       hubbardslices

struct FieldCorrelator
    primeslicemap::Dict # maps primes to slices
    correlator::Array
end

function FieldCorrelator(ecrystal::ElectronCrystal)
    
    basisvecs, _ = getVecs(ecrystal)
    latticepoints = getlatticepoints(getLattice(getCrystal(ecrystal)))
    ncells = size(latticepoints, 2)
    nbasis = size(basisvecs, 2)
    nspins = getStatesPerSite(ecrystal)

    # Since the diagonalization of the mean field requires that the two point correlator describing
    # the correlations between unit cell i and unit cell j depend only on their relative displacement.
    # The entire information of the two point correlator is contained in the set of two point correlators
    # between the unit cell at the origin and all other unit cells in the crystal.
    correlator = ones(ComplexF64, nspins, nbasis, ncells, nspins, nbasis)
    idx = CartesianIndices(correlator)
    primemap = hubbardslices(idx)
    initmfstruct!(correlator, primemap)
    FieldCorrelator(primemap, correlator)
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

function initmfstruct!(correlator::Array, primeslicemap::Dict)
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    for p in primes
        if !haskey(primeslicemap, p)
            break
        end
        
        itr = get(primeslicemap, p, nothing)
        for i in itr
            correlator[Tuple(i)...] *= p
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

#=
function hubbardslices(idxitr)
    pspin(i) = (i[1] == i[4]) && (i[2] == i[5]) && (i[3] == 1)  ? true : false
    apspin(i) = (i[1] != i[4]) && (i[2] == i[5]) && (i[3] == 1) ? true : false

    pspiniter = Iterators.filter(pspin, idxitr)
    apspiniter = Iterators.filter(apspin, idxitr)
    primeslicemap = Dict([(2, pspiniter), (3, apspiniter)])
    return primeslicemap
end
=#

function hubbardslices(idxitr)
    unoccupied(i) = ((2, 1, 1, 1, 1) == Tuple(i)) || ((1, 2, 1, 2, 2) == Tuple(i))
    occupied(i) = ((1, 1, 1, 1, 1) == Tuple(i)) || ((2, 2, 1, 2, 2) == Tuple(i))

    unoccupiedidxs = filter(unoccupied, idxitr)
    occupiedidxs = filter(occupied, idxitr)

    primeslicemap = Dict([(2, unoccupiedidxs), (3, occupiedidxs)])
    return primeslicemap
end

end