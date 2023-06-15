module Distributions
export fermidistribution,
       bosedistribution

function fermidistribution(e::Float64, mu::Float64, beta::Float64, )
    erel = e - mu
    ex = exp(-beta * erel)
    n = 1 / (ex + 1)
    return n
end

function bosedistribution(e::Float64, mu::Float64, beta::Float64)
    erel = e - mu
    ex = exp(-beta * erel)
    n = 1 / (ex - 1)
    return n
end

end