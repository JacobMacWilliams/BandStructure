tests = ["CrystalIOTest.jl", "CrystalInitTest.jl", "BandStructureTest.jl"]
for test in tests
    include(test)
end