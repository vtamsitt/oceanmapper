from map3d import map3d_surface
import numpy as np



# Load in surface to plot
cruise='33RR20090320'
fn = cruise+'_oxygen.npz'
d=np.load(fn)
xdata=d['lon']
ydata=d['lat']
zdata = d['depth']
scalardata = d['oxygen']


#user defined inputs

mode='rectangle'
vmin=150
vmax=210

#optional input arguments

map3d_surface(xdata,ydata,zdata,scalardata,mode,vmin,vmax)



