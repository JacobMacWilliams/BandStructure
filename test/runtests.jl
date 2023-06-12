tests = ["CrystalIOTest.jl", 
         "CrystalInitTest.jl", 
         "BandStructureTest.jl", 
         "ElectronLatticeTest.jl",
         "MeanFieldTest.jl",
         "EmbeddingsTest.jl"]
         
for test in tests
    include(test)
end