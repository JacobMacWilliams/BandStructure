module ChemicalPotential
using FromFile
using Roots
@from "../utils/Reciprocal.jl" using Reciprocal
@from "../utils/Distributions.jl" using Distributions
export findchempot

function chempot(bands::Int, fill::Float64, beta::Float64, mu::Float64, energies::Array{Float64})
    size = length(energies) / bands
    N = 0
    for e in energies
        N += fermidistribution(e, mu, beta)
    end
    N = N / size
    return N - fill * bands 
end

function findchempot(bands::Int, fill::Float64, beta::Float64, energies::Array{Float64})
    mumax = max(energies...)
    mumin = min(energies...)
    mufinder(mu) = chempot(bands, fill, beta, mu, energies)
    return find_zero(mufinder, (mumin, mumax), Bisection())
end

end