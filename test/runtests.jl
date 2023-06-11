tests = ["CrystalIOTest.jl", 
         "CrystalInitTest.jl", 
         "BandStructureTest.jl", 
         "ElectronLatticeTest.jl",
         "MeanFieldTest.jl"]
         
for test in tests
    include(test)
end