from map3d import map3d_surface
import numpy as np




#user defined inputs

mode='sphere'
lon_min = -80 
lon_max = 0
lat_min = 10
lat_max = 90
tlimits=[lon_min, lon_max, lat_min, lat_max]
zs  = 30000
vmin=150
vmax=210

#optional input arguments

mlab=map3d_surface(mode,zscale=zs,land_constant=True)
mlab.show()


