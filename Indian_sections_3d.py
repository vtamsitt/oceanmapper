#Henri F. Drake and Veronica Tamsitt 2017
#User inputs

#set input mode (cylinder,sphere or rectangle)
mode = 'rectangle'

#camera settings
azimuth = -59+180
elevation = 40
distance = 90

#slice of antarctic pie to plot (phi=latitude, theta=longitude)
phi_1 = -50
phi_2 = -25
theta_1 = 100
theta_2 = 150

#set vertical scaling
zscale=500

# set topo colormap
bathy_cmap = 'bone'

##############
import sys
import numpy as np
from mayavi import mlab
#from mpl_toolkits.basemap import Basemap
from scipy.io import loadmat
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib as mpl
import matplotlib.cm as cm
import IPython

#load custom colormap
d = loadmat('custom_PuYl.mat')
C = d['customPuYl']
cmap_PuYl = col.ListedColormap(C/255.0)

cmap_grey = col.LinearSegmentedColormap.from_list('own2',
                                                  ['lightgrey','lightgrey'])
cm.register_cmap(cmap=cmap_grey)


## Load in subsampled bathymetry from ETOPO1.
data = np.load('etopo1_30min.npz')
xraw = data['x']
yraw = data['y']
zraw = np.swapaxes(data['z'][:,:],0,1)
zraw[zraw>0]=0.
## Smoothing filter (scaling is how many grid cells to smooth over in each of two directions, so total smoothing factor is scaling**2)
#scaling = 2.0
#oceandepth = np.zeros((int(np.ceil(np.size(xraw)/scaling)),int(np.ceil(np.size(yraw)/scaling))))
#yt_ocean = np.zeros(int(np.ceil(np.size(yraw)/scaling)))
#xt_ocean = np.zeros(int(np.ceil(np.size(xraw)/scaling)))

# Smooth the bathymetry
#for i in range(np.size(xt_ocean)-1):
#    xt_ocean[i] = xraw[int(i*scaling)]

#for j in range(np.size(yt_ocean)):
#    yt_ocean[j] = yraw[int(j*scaling)]

#for j in range(np.size(yt_ocean)):
#    for i in range(np.size(xt_ocean)-1):
#        oceandepth[i,j] = np.average(zraw[int(i*scaling):int((i+1)*scaling),int(j*scaling):int((j+1)*scaling)])

phi = (yraw[:]*np.pi*2)/360.+np.pi/2.
theta = (xraw[:]*np.pi*2)/360.
c = zraw

phi_deg = data['y']
theta_deg = data['x']
theta=np.append(theta,theta[0])
c = np.concatenate((c,np.expand_dims(c[0,:],axis=0)),axis=0)


phi_ind1 = np.argmin(np.abs(phi_deg-phi_1))
phi_ind2 = np.argmin(np.abs(phi_deg-phi_2))
theta_ind1 = np.argmin(np.abs(theta_deg-theta_1))
theta_ind2 = np.argmin(np.abs(theta_deg-theta_2))


phi=phi[phi_ind1:phi_ind2]
theta=theta[theta_ind1:theta_ind2]
c = c[theta_ind1:theta_ind2:,phi_ind1:phi_ind2]
phi, theta = np.meshgrid(phi,theta)

# Scale back height of mountains for clarity
#c[c>-30]=c[c>-30]/3


# Create variable dimensions
if mode == 'sphere':
    x = np.sin(phi) * np.cos(theta[::-1]) * (1 + c/30000.)
    y = np.sin(phi) * np.sin(theta[::-1]) * (1 + c/30000.)
    z = np.cos(phi) * (1 + c/30000.)
    maxdep_lon = np.cos(phi[:,0]) * (1-(7000/30000.))
    maxdep_lat = np.cos(phi[0,:]) * (1-(7000/30000.))
elif mode == 'cylinder':
    x = np.sin(phi) * np.cos(theta[::-1])
    y = np.sin(phi) * np.sin(theta[::-1])
    z = c/15000.
    maxdep_lat = -7000/15000.
    maxdep_lon = -7000/15000.

elif mode == 'rectangle':
    y, x = np.meshgrid(phi_deg[phi_ind1:phi_ind2],theta_deg[theta_ind1:theta_ind2])
    z = c/500.
    maxdep = -7500/100.
    #y = y
    maxdep_lon = -7500/30000.
    maxdep_lat = -7500/30000.


mlab.figure(size = (1024,768),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
mlab.clf()
# Plot Bathymetry mesh
m = mlab.mesh(x, y, z, colormap=bathy_cmap,vmin=-7000/zscale,vmax=500/zscale)
#lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
#ilut=lut[::-1]
#m.module_manager.scalar_lut_manager.lut.table = ilut

#optional: plot constant color on land
sl = mlab.mesh(x, y, z,mask = z<0,color =(.7, .7, .71))
#Optional: plot bathymetry mesh at constant lon/lat

# Plot Lon / Dep bathymetry Mesh
#x_lon_dep = np.tile(x[:,-1],(2,1))
#y_lon_dep = np.tile(y[:,-1],(2,1))
#z_lon_dep = np.zeros_like(x_lon_dep)
#z_lon_dep[0,:] = z[:,-1]
#z_lon_dep[1,:] = maxdep_lon
#m = mlab.mesh(x_lon_dep, y_lon_dep, z_lon_dep, color=(0.5,0.5,0.5))
#lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
#ilut=lut[::-1]
#m.module_manager.scalar_lut_manager.lut.table = ilut

# Plot Lat / Dep bathymetry Meshes
#y_lat_dep = np.tile(y[-1,:],(2,1))
#x_lat_dep = np.tile(x[-1,:],(2,1))
#z_lat_dep = np.zeros_like(x_lat_dep)
#z_lat_dep[0,:] = z[-1,:]
#z_lat_dep[1,:] = maxdep_lat
#m = mlab.mesh(x_lat_dep, y_lat_dep, z_lat_dep, colormap='Greys')
#lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
#ilut=lut[::-1]
#m.module_manager.scalar_lut_manager.lut.table = ilut

#y_lat_dep = np.tile(y[0,:],(2,1))
#x_lat_dep = np.tile(x[0,:],(2,1))
#z_lat_dep[0,:] = z[0,:]
#m = mlab.mesh(x_lat_dep, y_lat_dep, z_lat_dep, colormap ='Greys')
#lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
#ilut=lut[::-1]
#m.module_manager.scalar_lut_manager.lut.table = ilut


cmappu='YlGnBu'
# Load in surface to plot
cruise='s05_120E'
fn = cruise+'_oxygen.npz'
d=np.load(fn)
xt_ocean=d['lon']
yt_ocean=d['lat']
depth_h = d['depth']
scalar = d['oxygen'] #2D field to color on mesh surface

phi_iso, theta_iso = np.meshgrid(((yt_ocean*np.pi*2)/360.)+np.pi/2.,(xt_ocean*np.pi*2)/360.)

# Create variable dimensions
if mode == 'sphere':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
elif mode == 'cylinder':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
    z_iso = -depth_h/15000.

elif mode == 'rectangle':
    y_iso,z_iso = np.meshgrid(yt_ocean,depth_h)
    x_iso,z_iso = np.meshgrid(xt_ocean,depth_h) 
    #y_iso = y_iso*3
    z_iso =-z_iso/500.
#    z_iso = -np.tile(depth_h,(depth_h.shape[0],y_iso.shape[0]))/100.
print scalar
print cmappu
m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalar,colormap=cmappu,vmin =150,vmax=210,opacity=1,mask=np.isnan(scalar))
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

# Load in surface to plot
cruise='s05_132E'
fn = cruise+'_oxygen.npz'
d=np.load(fn)
xt_ocean=d['lon']
yt_ocean=d['lat']
depth_h = d['depth']
scalar = d['oxygen'] #2D field to color on mesh surface

phi_iso, theta_iso = np.meshgrid(((yt_ocean*np.pi*2)/360.)+np.pi/2.,(xt_ocean*np.pi*2)/360.)

# Create variable dimensions
if mode == 'sphere':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
elif mode == 'cylinder':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
    z_iso = -depth_h/15000.

elif mode == 'rectangle':
    y_iso,z_iso = np.meshgrid(yt_ocean,depth_h)
    x_iso,z_iso = np.meshgrid(xt_ocean,depth_h) 
    #y_iso = y_iso*3
    z_iso =-z_iso/500.
#    z_iso = -np.tile(depth_h,(depth_h.shape[0],y_iso.shape[0]))/100.
m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalar,colormap=cmappu,vmin =150,vmax=210,opacity=1)
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]


# Load in surface to plot
cruise='i09s_20041223'
fn = cruise+'_oxygen.npz'
d=np.load(fn)
xt_ocean=d['lon']
yt_ocean=d['lat']
depth_h = d['depth']
scalar = d['oxygen'] #2D field to color on mesh surface

phi_iso, theta_iso = np.meshgrid(((yt_ocean*np.pi*2)/360.)+np.pi/2.,(xt_ocean*np.pi*2)/360.)

# Create variable dimensions
if mode == 'sphere':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
elif mode == 'cylinder':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
    z_iso = -depth_h/15000.

elif mode == 'rectangle':
    y_iso,z_iso = np.meshgrid(yt_ocean,depth_h)
    x_iso,z_iso = np.meshgrid(xt_ocean,depth_h) 
    #y_iso = y_iso*3
    z_iso =-z_iso/500.
#    z_iso = -np.tile(depth_h,(depth_h.shape[0],y_iso.shape[0]))/100.
m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalar,colormap=cmappu,vmin =150,vmax=210,opacity=1)
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

# Load in surface to plot
cruise='33RR20090320'
fn = cruise+'_oxygen.npz'
d=np.load(fn)
xt_ocean=d['lon']
yt_ocean=d['lat']
depth_h = d['depth']
scalar = d['oxygen'] #2D field to color on mesh surface

phi_iso, theta_iso = np.meshgrid(((yt_ocean*np.pi*2)/360.)+np.pi/2.,(xt_ocean*np.pi*2)/360.)

# Create variable dimensions
if mode == 'sphere':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
elif mode == 'cylinder':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
    z_iso = -depth_h/15000.

elif mode == 'rectangle':
    y_iso,z_iso = np.meshgrid(yt_ocean,depth_h)
    x_iso,z_iso = np.meshgrid(xt_ocean,depth_h) 
    #y_iso = y_iso*3
    z_iso =-z_iso/500.
#    z_iso = -np.tile(depth_h,(depth_h.shape[0],y_iso.shape[0]))/100.
m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalar,colormap=cmappu,vmin =150,vmax=210,opacity=1)
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]


#lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
#lut[:,0:3]=256*np.asarray(cmap_PuYl.colors)
#m.module_manager.scalar_lut_manager.lut.table = lut

c2 = mlab.scalarbar(m, nb_labels = 7, orientation = 'vertical',label_fmt='%.0f')
c2.scalar_bar_representation.proportional_resize=True
c2.scalar_bar_representation.position = [0.05, 0.05]
c2.scalar_bar_representation.position2 = [0.08, 0.9]
#c2.label_text_property.font_family = 'courier'
#c2.label_text_property.font_size = 10

#add trajectory
#d=np.load('/Users/vtamsitt/Dropbox/SO/tamsitt/Figures/SOSE_tracer/sample_trajectories.npz')
#lons_plot=d['lons_plot'][:,4]
#lats_plot=d['lats_plot'][:,4]
#depth_h=d['depths_plot'][:,4]
#phi_iso, theta_iso = np.meshgrid(((lats_plot*np.pi*2)/360.)+np.pi/2.,(lons_plot*np.pi*2)/360.)
#
## Create variable dimensions
#if mode == 'sphere':
#    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
#    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
#    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
#elif mode == 'cylinder':
#    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
#    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
#    z_iso = -depth_h/15000.
#
#elif mode == 'rectangle':
#    y_iso,z_iso = np.meshgrid(lats_plot,depth_h)
#    x_iso,z_iso = np.meshgrid(lons_plot,depth_h) 
#    #y_iso = y_iso*3
#    z_iso =-z_iso/200.
##print lons_plot
#print lats_plot
#print depth_h.shape
#mlab.plot3d(x_iso,y_iso,z_iso, -z_iso*200.,opacity=0.7,tube_radius=0.4,tube_sides=12,colormap='YlGnBu',vmin=500,vmax=3000)
#
#add trajectory
d=np.load('Indian_trajectories.npz')
ind=1 #which traj to plot?
lons_plot=d['xp'][ind,:]
lats_plot=d['yp'][ind,:]
depth_h=d['zp'][ind,:]
phi_iso, theta_iso = np.meshgrid(((lats_plot*np.pi*2)/360.)+np.pi/2.,(lons_plot*np.pi*2)/360.)

# Create variable dimensions
if mode == 'sphere':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/30000.)
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -depth_h/30000.)
    z_iso = np.cos(phi_iso) * (1 -depth_h/30000.)
elif mode == 'cylinder':
    x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
    y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
    z_iso = -depth_h/15000.

elif mode == 'rectangle':
    y_iso,z_iso = np.meshgrid(lats_plot,depth_h)
    x_iso,z_iso = np.meshgrid(lons_plot,depth_h) 
    #y_iso = y_iso*3
    z_iso =-z_iso/500.
cutoff=1900 #cutoff index for trajectory
mlab.plot3d(lons_plot[:cutoff],lats_plot[:cutoff],-depth_h[:cutoff]/500., depth_h[:cutoff],opacity=0.7,tube_radius=0.1,tube_sides=15,colormap='YlOrRd',vmin=1000,vmax=3000)


mlab.view(azimuth = azimuth, elevation = elevation, distance = distance)
#IPython.embed()
mlab.show()

