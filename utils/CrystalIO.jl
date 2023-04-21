module CrystalIO
using FromFile
using ConfParser
using Parsers: parse
export parseLatticeConf, parseCrystal, 
writeToBravaisConf, writeToCrystalConf

function parseLatticeConf(file::String, sec::String)
    conf = ConfParse(file, "ini")
    parse_conf!(conf)
  
    # get and store config parameters
    dim = parse(Int, retrieve(conf, sec, "dim"))
    vecs = zeros(dim, dim)
    for i = 1:dim
      v = retrieve(conf, sec, string("v",i))
      for j in eachindex(v)
        vecs[j, i] = parse(Float64, v[j])
      end
    end
    return dim, vecs
end

function parseCrystalConf(file::String, sec::String)
  conf = ConfParse(file, "ini")
  parse_conf!(conf)

  # get and store config parameters
  N = parse(Int, retrieve(conf, sec, "N"))
  lattice = retrieve(conf, sec, "lattice")
  dim = parse(Int, retrieve(conf, sec, "dim"))
  basis = parse(Int, retrieve(conf, sec, "basis"))

  basisvecs = zeros(dim, basis)
  for i = 1:basis
    v = retrieve(conf, sec, string("v",i))
    for j in eachindex(v)
      basisvecs[j, i] = parse(Float64, v[j])
    end
  end
  return N, basisvecs, lattice
end


function writeToBravaisConf(name, vecs::Matrix{Float64}, path::String)
    conf = ConfParse(path, "ini")
    parse_conf!(conf)

    dim = size(vecs, 1)
    commit!(conf, name, "dim", string(dim))

    for i in 1:dim
      keyName = string("v", i)
      values = Array{String, 1}(undef, dim)
      for j in 1:dim
          values[j] = string(vecs[j, i])
      end
      commit!(conf, name, keyName, values)
    end
    save!(conf, path)
end

function writeToCrystalConf(name, N::Integer, lattice::String, basisvecs::Matrix{Float64}, path::String)
    dim = size(basisvecs, 1)
    basis = size(basisvecs, 2)

    conf = ConfParse(path, "ini")
    parse_conf!(conf)

    commit!(conf, name, "N", string(N))
    commit!(conf, name, "lattice", lattice)
    commit!(conf, name, "dim", string(dim))
    commit!(conf, name, "basis", string(basis))


    for i in 1:basis
        key = string("v", i)
        values = Array{String, 1}(undef, dim)
        for j in 1:dim
            values[j] = string(basisvecs[j, i])
        end
        commit!(conf, name, key, values)
    end
    save!(conf, path)
end

end