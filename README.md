# wave-equation

Simulation of acoustic wave propagation through a water-earth velocity model. The 2D synthetic model consists of a water layer and a sediment layer. The water layer has a velocity of 1.5km/s. The sea floor is at 500m. Below, there is a layer of sediment of velocity 2.0km/s. The simulation extends to a depth of 3.0km and a width of 6.0km.

This program uses the 2D scalar acoustic wave equation, with a second-order time and fourth-order space discretisation.

The seismic source is located at the surface, halfway along the model. The source term is a Ricker wavelet with centre frequency 20Hz.

Purely reflecting boundary conditions are assumed.

The code can easily be adapted for different velocity models, see `wave.f90`. Compile with `gfortran -Wall -Werror -O3 -fopenmp -o wave wave.f90`.

<p align="center">
  <img src="https://github.com/sam-marsh/wave-equation/blob/master/output.gif?raw=true">
</p>

In addition, see *Mathematica* notebook `GFD-Coefficients.nb` for derivation of the generalised finite difference formula and coefficients used in this program.
