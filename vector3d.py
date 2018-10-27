#import dependent modules
import sys
import numpy as np
from mayavi import mlab



def vector3d(mode,xdata,ydata,zdata,udata,vdata,wdata,fig=None,zscale=500.,opacity=1.0,quiver_mode='2darrow', quiver_scale=1, quiver_spacing=8., set_view=None):
    """
    fig: integer or string, optional. Figure key will plot data on corresponding mlab figure, if it exists, or create a new one
    mode: string; coordinate system of 3D projection. Options are 'rectangle' (default), 'spherical' or 'cylindrical'
    xdata: 1D numpy array; longitude values for data array
    ydata: 1D numpy array; latitude values for data array
    zdata: 1D numpy array; depth values for data array
    scalardata: 2D numpy array, optional; 2D scalar field to plot colors on surface
    vmin: float, optional; colorbar minimum for data
    vmax: float, optional; colorbar maximum for data
    data_cmap: string, optional; colormap for data surface, default is blue-red
    data_alpha: float or int, optional; opacity for data surface from 0 to 1, default is 1
    set_view: array_like, optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is 
    """
        
    
    #make figure
    if fig is None: 
        mlab.figure(size = (1024,768),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
        mlab.clf()
    else:
        mlab.figure(figure=fig)
    
    #do coordinate transformation
    if xdata is not None and ydata is not None and zdata is not None:
        #TODO add an error message if not all data fields are provided
        #prep data grid
        phi_iso, theta_iso = np.meshgrid(((ydata*np.pi*2)/360.)+np.pi/2.,(xdata*np.pi*2)/360.)

        if mode is 'sphere':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -depth_h/zscale)
            y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -zdata/zscale)
            z_iso = np.cos(phi_iso) * (1 -zdata/zscale)
        elif mode is 'cylinder':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
            y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
            z_iso = zdata/zscale

        elif mode is 'rectangle':
            y_iso,z_iso = np.meshgrid(ydata,zdata)
            x_iso,z_iso = np.meshgrid(xdata,zdata)
            z_iso =-z_iso/zscale
    else:
        #raise error if all three fields are not provided
        print 'ERROR: not all data fields are provided. Must provide 1D data x, y and z data points'  
    
    #do quiver plot 
    mlab.quiver3d(x_iso, y_iso, z_iso, udata, vdata, wdata, mode=quiver_mode,opacity=opacity,scale_factor=quiver_scale,mask_points=quiver_spacing)   
 
    #optional: change mayavi camera settings
    if set_view is None:
        mlab.view(distance = 'auto')
    else:
        mlab.view(azimuth = set_view[0], elevation = set_view[1], distance = set_view[2], focalpoint = set_view[3])


    return mlab.gcf()

