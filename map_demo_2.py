from map3d import map3d_surface
import numpy as np
from IPython import embed




#user defined inputs
mode='sphere'

#optional input arguments
lon_min = -80 
lon_max = 0
lat_min = 10
lat_max = 90
tlimits=[lon_min, lon_max, lat_min, lat_max]
zs  = 40000
tvmin = -100 
tvmax = 7000
lcolor = (0.5, 0.6, 0.5)

mlab=map3d_surface(mode,zscale=zs,topo_vmin = tvmin, topo_vmax= tvmax, topo_cmap_reverse=True,land_constant=True,land_color=lcolor)
embed()
mlab.show()

