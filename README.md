Simulation of acoustic wave propagation through a water-earth velocity model. The 2D synthetic model consists of a water layer and a sediment layer. The water layer has a velocity of 1.5km/s. The sea floor is at 500m. Below, there is a layer of sediment of velocity 2.0km/s. The simulation extends to a depth of 3.0km and a width of 6.0km.

The 2D scalar acoustic wave equation is used, with a second-order time and fourth-order space discretisation.

The seismic source is located at the surface, halfway along the model. The source term is a Ricker wavelet with centre frequency 20Hz.

Reflecting boundary conditions are used.

![simulation](https://github.com/sam-marsh/wave-equation/blob/master/output.gif?raw=true)
