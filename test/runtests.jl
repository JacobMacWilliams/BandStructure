tests = ["CrystalInitTest.jl", 
         "BandStructureTest.jl", 
         "ElectronLatticeTest.jl",
         "MeanFieldTest.jl",
         "EmbeddingsTest.jl",
         "ChemicalPotentialTest.jl"]
         
for test in tests
    include(test)
end