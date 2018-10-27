# oceanmapper
generate 3D maps and animations of ocean bathymetry and data using mayavi

The module map3d takes ETOPO1 land surface topography and ocean bathymetry (doi:10.7289/V5C8276M) at 30 arc minute resolution (or an alternative user input topography product) and three dimensional data and maps it in three possible different 3D proejection modes: 'sphere', 'cylinder' or 'rectangle'.

The key functionality of the code is to take the input x,y,z data values and convert them to the correct coordinate system for the chosen mode, before implementing a mayavi plotting function. 

