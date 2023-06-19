module TOMLIO
using TOML
export parseconfig,
       parsematrix

function parseconfig(filename, cmptype, cmpname)
    config = TOML.parsefile(filename)
    if !haskey(config, cmptype)
        error("An appropriate entry for " * cmptyp * "Is not available.")
    end

    configarray = get(config, cmptype, nothing)
    for cmp in configarray
        name = get(cmp, "name", nothing)
        if name == cmpname
            return cmp
        end
    end
    return nothing
end

function parsematrix(config::Dict, prefix::String)
    i = 1
    vectors = []
    while haskey(config, string(prefix, i))
        vec = get(config, string(prefix, i), nothing)
        push!(vectors, vec)
        i += 1
    end
    vectors = cat(vectors..., dims=2)
    return vectors
end

end