from map3d import topo_surface3d
import numpy as np



# Load in sample data surface to plot
fn = 'i09s_20041223_oxygen.npz'
d=np.load(fn)
xdata=d['lon']
ydata=d['lat']
zdata = d['depth']
scalardata = d['oxygen']


#user defined inputs
mode='rectangle'
vmin=150
vmax=210
lat_min = -55
lat_max = -25
lon_min = 100
lon_max = 150

#run function
mlab=map3d_surface(mode,xdata,ydata,zdata,scalardata,vmin=vmin,vmax=vmax,topo_limits=[lon_min, lon_max, lat_min, lat_max],topo_cmap_reverse=True)
mlab.show()


