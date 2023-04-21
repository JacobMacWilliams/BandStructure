module BandStructure
println("test")
#=
using ConfParser
using Base.Filesystem: joinpath

include("BravaisLattice.jl")
using .BravaisLattice: BravaisLattice
global CONFIGPATH = "configs.ini"
conf = ConfParse(CONFIGPATH)
parse_conf!(conf)
confroot = retrieve(conf, "Configurations", "root")
confpath = joinpath(confroot, "bravais.ini")
b = BravaisLattice("Triangle", confpath)
=#
end