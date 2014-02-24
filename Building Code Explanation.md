# VEFIE Building Code Explanation

1. *run_solution.m*  
Calls *variables.m*, *create_shape.m*, *create_vefie_elements.m* and *vefie_solution.m*. It plots the Total Electric Field, the Incident Field and generates a movie using multiplication of the wave by $$$e^{jwt}$$$.

2. *variables.m*  
Most (probably all) of the variables used are defined in here.
	* The frequency of the wave
	* The various material properties, for free space, concrete, glass and wood
	* Number of discretisations per lambda
	* Size of the shape in x and y lengths
		+ Wall and door sizes are computed as a fraction of these
	* The antenna location

3. *create_shape.m*  
This sets up an outline for the shape and computes some variables like N and M before calling *create_building.m*. *create_building.m* puts in the walls around the outside of the building, creates a few rooms and adds gaps for doors. Various vectors are then created with the position of each point in the building and its wave number before the building is plotted.

4. *create_vefie_elements.m*  
Simply this creates the V, G and D vectors for the problem.

5. *vefie_solution.m*  
This solves the equation $$$V = ZE$$$ where $$$Z = I + GD$$$ using the CGNE-FFT and the reduced forward operator. At various points this calls *conv_fft.m* which performs convolution using the Fast Fourier Transform.

There are various comments throughout the code to which should clarify a bit more in-depth what is going on. If anything isn't clear let me know and I'll have a look.