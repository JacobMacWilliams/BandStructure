# BandStructure
Calculate and plot the band structure of a simple graphene lattice.

## Usage
Simply create a plots folder under the root of this project and run BandStructure.jl in src.

## Extendability
This package can be extended to graph the band structure of all crystal lattices with or without lattice imbalances between the lattice sites so long as the number of points in the crystals basis is not too large (The exact upper limit is unknown by the author of this paper but standard crystals with a number of points in there basis, on the order of ten, is surely non-problematic). As of now this package can plot the band structure of any such crystal without lattice imbalances and an isotropic lattice constant so long the real and reciprocal lattices are known (though methods must be added to create the reciprocal lattice and trace out a highly symmetric path, after which calculateGraphenBandStructure can be used simply be replacing the call to getGraphenePath with the just created method).

## Contact
For more information please email: jacob.macwilliams@uni-konstanz.de
