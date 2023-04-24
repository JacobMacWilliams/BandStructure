module CrystalWriterTest
using FromFile
@from "../utils/CrystalIO.jl" using CrystalIO: 
writeToBravaisConf as bravaisWriter,
writeToCrystalConf as crystalWriter, 
parseLatticeConf as bravaisParser, 
parseCrystalConf as crystalParser

using ConfParser, Test

function bravaisReaderWriterTest(confpath::String)
    open(confpath, truncate = true)
    vecs = zeros(2, 2)
    vecs[1,1] = 1/2 * sqrt(3)
    vecs[2,1] = 1/2 
    vecs[1,2] = 1/2 * sqrt(3)
    vecs[2,2] = -1/2
    bravaisWriter("Triangle", vecs, confpath)
    dim, vecs_read = bravaisParser(confpath, "Triangle")
    
    id = (dim == 2) && (vecs_read == vecs)
    return id
end

function crystalReaderWriterTest(confpath::String)
    open(confpath, truncate = true)
    basevecs = zeros(2,2)
    basevecs[1,1] = 0
    basevecs[2,1] = 0
    basevecs[1,2] = 1
    basevecs[2,2] = 0
    crystalWriter("graphene", 400, "triangle", basevecs, confpath)
    N, basevecs_read, lattice = crystalParser(confpath, "graphene")

    id = ((N == 400) && (lattice == "triangle") && (basevecs_read == basevecs))
    return id

end

bravaisconf = joinpath("conf", "bravais.ini")
crystalconf = joinpath("conf", "crystal.ini")
@test bravaisReaderWriterTest(bravaisconf)
@test crystalReaderWriterTest(crystalconf)

end