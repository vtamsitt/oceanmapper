#import dependent modules
import sys
import numpy as np
#from mayavi import mlab
#mlab.init_notebook() #allows code to run within jupyter notebook

def topography3d(mode,topo_x=None,topo_y=None,topo_z=None,topo_limits=None,zscale=500.,topo_vmin=None,topo_vmax=None,topo_cmap='bone',topo_cmap_reverse=False,land_constant=False,land_color=(0.7,0.7,0.7),set_view=None):
    """
    mode: string; coordinate system of 3D projection. Options are 'rectangle' (default), 'sphere' or 'cylinder'
    topo: array_like, optional; input topography file, default is etopo 30 
##TODO: need to define sign of topography 
    topo_limits: array_like, optional; longitude and latitude limits for 3d topography plot [lon min, lon max, lat min, lat max], longitudes range -180 to 180, latitude -90 to 90, default is entire globe
    zscale: scalar, optional; change vertical scaling for plotting, such that the vertical axis is scaled as topo_z/zscale (assumes topo_z units are m); default zscale is 500
    topo_cmap: string, optional; set colormap for topography, default is bone 
    topo_cmap_reverse: string, optional; reverse topography colormap, default is false
    land_constant: string optional; if True, land is set to one colour, default is False
    land_color: color, optional; RGB triplet specifying land colour, defauly is gret
    set_view: array_like, optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is None 
    """
        
    #load topo data
    if topo_x is not None and topo_y is not None and topo_z is not None:
        xraw = topo_x
        yraw = topo_y
        zraw = topo_z

    else:
        tfile = np.load('etopo1_30min.npz')
        xraw = tfile['x']
        yraw = tfile['y']
        zraw = np.swapaxes(tfile['z'][:,:],0,1)
    

    #create coordinate variables
    phi = (yraw[:]*np.pi*2)/360.+np.pi/2.
    theta = (xraw[:]*np.pi*2)/360.
    c = zraw
    theta=np.append(theta,theta[0])
    c = np.concatenate((c,np.expand_dims(c[0,:],axis=0)),axis=0)

    if topo_vmin is None:
        tvmin = 0
    else:
        tvmin = topo_vmin
    if topo_vmax is None:
        tvmax = 7000
    else:
        tvmax = topo_vmax
   
    if topo_limits is not None:
        phi_1 = topo_limits[2]
        phi_2 = topo_limits[3]
        theta_1 = topo_limits[0]
        theta_2 = topo_limits[1]

        phi_ind1 = np.argmin(np.abs(yraw-phi_1))
        phi_ind2 = np.argmin(np.abs(yraw-phi_2))
        theta_ind1 = np.argmin(np.abs(xraw-theta_1))
        theta_ind2 = np.argmin(np.abs(xraw-theta_2))

        #restrict topo extent
        phi=phi[phi_ind1:phi_ind2]
        theta=theta[theta_ind1:theta_ind2]
        c = c[theta_ind1:theta_ind2:,phi_ind1:phi_ind2]


def topo_surface3d(mode,xdata=None,ydata=None,zdata=None,scalardata=None,vmin=None,vmax=None,data_cmap='blue-red',data_alpha=1,topo_x=None,topo_y=None,topo_z=None,topo_limits=None,zscale=500.,topo_vmin=None,topo_vmax=None,topo_cmap='bone',topo_cmap_reverse=False,land_constant=False,land_color=(0.7,0.7,0.7),set_view=None):
    '''
    mode = (string) coordinate system of 3D projection. Options are 'rectangle' (default), 'sphere' or 'cylinder'
    xdata = optional; (1D numpy array) longitude values for data array
    ydata = optional; (1D numpy array) latitude values for data array
    zdata = optional; (1D numpy array) depth values for data array
    scalardata = optional; (2D numpy array) scalar field to plot colors on surface
    vmin = (float) colorbar minimum for data
    vmax = (float) colorbar maximum for data
    data_cmap = colormap for data surface, default is blue-red
    data_alpha = (float or int) opacity for data surface from 0 to 1, default is 1
    topo = optional; input topography file, default is etopo 30 
    topo_limits = optional; longitude and latitude limits for 3d topography plot [lon_min, lon_max, lat_min, lat_max], longitudes range -180 to 180, latitude -90 to 90, default is entire globe
    zscale: scalar, optional; change vertical scaling for plotting, such that the vertical axis is scaled as topo_z/zscale (assumes topo_z units are m); default zscale is 500 
    topo_cmap = optional; default is bone 
    topo_cmap_reverse = optional; reverse topography colormap, default is false
    set_view = optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is 
    '''
    #TODO expand/clean descriptions
        
    #load topo data
    if topo_x is not None and topo_y is not None and topo_z is not None:
        xraw = topo_x
        yraw = topo_y
        zraw = topo_z

    else:    
        tfile = np.load('etopo1_30min.npz')
        xraw = tfile['x']
        yraw = tfile['y']
        zraw = np.swapaxes(tfile['z'][:,:],0,1)
    
    
    phi = (yraw[:]*np.pi*2)/360.+np.pi/2.
    theta = (xraw[:]*np.pi*2)/360.
    c = zraw
    theta=np.append(theta,theta[0])
    c = np.concatenate((c,np.expand_dims(c[0,:],axis=0)),axis=0)

    if topo_limits is not None:
        phi_1 = topo_limits[2]
        phi_2 = topo_limits[3]
        theta_1 = topo_limits[0]
        theta_2 = topo_limits[1]

        phi_ind1 = np.argmin(np.abs(yraw-phi_1))
        phi_ind2 = np.argmin(np.abs(yraw-phi_2))
        theta_ind1 = np.argmin(np.abs(xraw-theta_1))
        theta_ind2 = np.argmin(np.abs(xraw-theta_2))

        #restrict topo extent
        phi=phi[phi_ind1:phi_ind2]
        theta=theta[theta_ind1:theta_ind2]
        c = c[theta_ind1:theta_ind2:,phi_ind1:phi_ind2]
    phi, theta = np.meshgrid(phi,theta)
   


    if topo_vmin is None:
        tvmin = 0
    else:
        tvmin = topo_vmin
    if topo_vmax is None:
        tvmax = 7000
    else:
        tvmax = topo_vmax
    
    #make figure
    #mlab.figure(size = (1024,768),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
    #mlab.clf()
    # Plot Bathymetry mesh
    if mode is 'sphere':
        x = np.sin(phi) * np.cos(theta[::-1]) * (1 + c/zscale)
        y = np.sin(phi) * np.sin(theta[::-1]) * (1 + c/zscale)
        z = np.cos(phi) * (1 + c/zscale)
    
    elif mode is 'cylinder':
        x = np.sin(phi) * np.cos(theta[::-1])
        y = np.sin(phi) * np.sin(theta[::-1])
        z = c/zscale
    
    elif mode is 'rectangle':
        if topo_limits is not None:
            y, x = np.meshgrid(yraw[phi_ind1:phi_ind2],xraw[theta_ind1:theta_ind2])
            z = c/zscale
        else:
            y, x = np.meshgrid(yraw,xraw)
            z = c/zscale
    else:
        print('mode is not valid. Must be \'sphere\',\'cylinder\', or \'rectangle\'')    
    #make bathymetry mesh
    m = mlab.mesh(x, y, z, scalars = -c, colormap=topo_cmap,vmin=tvmin,vmax=tvmax)
    
    #optional: reverse bathymetry colormap
    if topo_cmap_reverse is True:
        lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
        ilut=lut[::-1]
        m.module_manager.scalar_lut_manager.lut.table = ilut

    #optional: plot constant color on land
    if land_constant is True:
        sl = mlab.mesh(x, y, z,mask = c<0,color =land_color)


    #optional: plot data surface
    if xdata is not None and ydata is not None and zdata is not None:
        #TODO add an error message if not all data fields are provided
        #prep data grid
        phi_iso, theta_iso = np.meshgrid(((ydata*np.pi*2)/360.)+np.pi/2.,(xdata*np.pi*2)/360.)
 
        if mode is 'sphere':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -zdata/zscale)
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
        
        if scalardata is not None:
            m1 = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalardata,colormap=data_cmap,vmin =vmin,vmax=vmax,opacity=data_alpha)
            m1.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]
        else:
            m1 = mlab.mesh(x_iso, y_iso, z_iso,vmin =vmin,vmax=vmax,opacity=data_alpha)
             
    #optional: change mayavi camera settings
    if set_view is None:
        mlab.view(distance = 'auto')
    else:
        mlab.view(azimuth = set_view[0], elevation = set_view[1], distance = set_view[2], focalpoint = set_view[3])


    return m


def surface3d(mode,xdata,ydata,zdata,fig=None,scalardata=None,vmin=None,vmax=None,data_cmap='blue-red',data_color=(0.5,0.5,0.5),data_alpha=1,zscale=500.,set_view=None):
    """
    fig: integer or string, optional. Figure key will plot data on corresponding mlab figure, if it exists, or create a new one
    mode: string; coordinate system of 3D projection. Options are 'rectangle' (default), 'sphere' or 'cylinder'
    xdata: 1D numpy array; longitude values for data array
    ydata: 1D numpy array; latitude values for data array
    zdata: 1D numpy array; depth values for data array
    scalardata: 2D numpy array, optional; 2D scalar field to plot colors on surface
    vmin: float, optional; colorbar minimum for data
    vmax: float, optional; colorbar maximum for data
    zscale: scalar, optional; change vertical scaling for plotting, such that the vertical axis is scaled as topo_z/zscale (assumes topo_z units are m); default zscale is 500 
    data_cmap: string, optional; colormap for data surface, default is blue-red
    data_color: triplet of floats ranging from 0 to 1, optional; sets color of surface, overrides colormap when scalardata is not included
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
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -zdata/zscale)
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
    
    #map data surface
    if scalardata is not None:
        m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalardata,colormap=data_cmap,vmin =vmin,vmax=vmax,opacity=data_alpha)
        m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

    else:
        m = mlab.mesh(x_iso, y_iso, z_iso,color=data_color,vmin =vmin,vmax=vmax,opacity=data_alpha)
        m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

def vector3d(mode,xdata,ydata,zdata,udata,vdata,wdata,scalardata=None,fig=None,zscale=500.,vector_color=(0,0,0),vector_cmap=None,alpha=1.0,vector_mode='2darrow', scale=1, spacing=8., set_view=None):
    """
    fig: integer or string, optional. Figure key will plot data on corresponding mlab figure, if it exists, or create a new one
    mode: string; coordinate system of 3D projection. Options are 'rectangle' (default), 'sphere' or 'cylinder'
    xdata: 1D array; longitude values for data array
    ydata: 1D array; latitude values for data array
    zdata: 1D array; depth values for data array
    udata: 2D or 3D array; u vector component
    vdata: 2D or 3D array; v vector component
    wdata: 2D or 3D array; w vector component
    zscale: scalar, optional; change vertical scaling for plotting, such that the vertical axis is scaled as topo_z/zscale (assumes topo_z units are m); default zscale is 500 
    vector_mode: string, optional; style of vector plot
    color: colormap or rgb triplet,optional; color of quiver plot default is black (0,0,0). 
    alpha: float or int, optional; opacity for data surface from 0 to 1, default is 1
    scale: float or int, optional; scaling for length of vectors, default is 1. 
    spacing: int, optional; If supplied, only one out of 'spacing' data points is displayed. This option is useful to reduce the number of points displayed on large datasets Must be an integer (int or long) or None
    set_view: array_like, optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is 
    """
        
    
    #make figure
    if fig is None: 
        mlab.figure(size = (1024,768),bgcolor = (1,1,1))
        mlab.clf()
    else:
        mlab.figure(figure=fig,bgcolor = (1,1,1))
    
    #do coordinate transformation
    if xdata is not None and ydata is not None and zdata is not None:
        #TODO add an error message if not all data fields are provided
        #prep data grid
        phi_iso, theta_iso = np.meshgrid(((ydata*np.pi*2)/360.)+np.pi/2.,(xdata*np.pi*2)/360.)

        if mode is 'sphere':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -zdata/zscale)
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
    if scalardata is not None:
        mlab.quiver3d(x_iso, y_iso, z_iso, udata, vdata, wdata, scalars=scalardata, scale_mode=None,colormap=vector_cmap,mode=vector_mode,opacity=alpha,scale_factor=scale,mask_points=spacing)   
    elif vector_cmap is not None:
        mlab.quiver3d(x_iso, y_iso, z_iso, udata, vdata, wdata, colormap=vector_cmap,mode=vector_mode,opacity=alpha,scale_factor=scale,mask_points=spacing)   
    else:
        mlab.quiver3d(x_iso, y_iso, z_iso, udata, vdata, wdata, color=vector_color,mode=vector_mode,opacity=alpha,scale_factor=scale,mask_points=spacing)   
 
    #optional: change mayavi camera settings


def trajectory3d(mode,xdata,ydata,zdata,fig=None,scalardata=None,vmin=None,vmax=None,color=(0,0,0),data_cmap=None,data_alpha=1,zscale=500.,tube_radius=0.01,tube_sides=15,set_view=None):
    """
    fig: integer or string, optional. Figure key will plot data on corresponding mlab figure, if it exists, or create a new one
    mode: string; coordinate system of 3D projection. Options are 'rectangle' (default), 'sphere' or 'cylinder'
    xdata: 1D array; longitude values for data array
    ydata: 1D array; latitude values for data array
    zdata: 1D array; depth values for data array
    scalardata: 1D array, optional; 1D scalar field to plot colors along trajectoy
    zscale: scalar, optional; change vertical scaling for plotting, such that the vertical axis is scaled as topo_z/zscale (assumes topo_z units are m); default zscale is 500 
    vmin: float, optional; colorbar minimum for data
    vmax: float, optional; colorbar maximum for data
    data_cmap: string, optional; colormap for data surface, default is blue-red
    data_alpha: float or int, optional; opacity for data surface from 0 to 1, default is 1
    tube_radius: float, optional; radius of tube
    tube_sides: int, optional; number of sides of tube
    set_view: array_like, optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is 
    """
        
    
    #make figure
    if fig is None: 
        mlab.figure(size = (1024,768),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
        mlab.clf()
    else:
        mlab.figure(figure=fig,bgcolor = (1,1,1))
    
    #do coordinate transformation

    if xdata is not None and ydata is not None and zdata is not None:
        phi_iso, theta_iso = np.meshgrid(((ydata*np.pi*2)/360.)+np.pi/2.,(xdata*np.pi*2)/360.)

        # Create variable dimensions
        if mode == 'sphere':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -zdata/zscale)
            y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -zdata/zscale)
            z_iso = np.cos(phi_iso) * (1 -zdata/zscale)
        elif mode == 'cylinder':
            x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
            y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
            z_iso = -zdata/zscale
        elif mode == 'rectangle':
            x_iso = xdata
            y_iso = ydata
            z_iso =-zdata/zscale
    else:
        #raise error if all three fields are not provided
        print 'ERROR: not all data fields are provided. Must provide 1D data x, y and z data points'  
    

    #map data surface
    if scalardata is not None:
        mlab.plot3d(x_iso,y_iso,z_iso, scalardata,opacity=data_alpha,tube_radius=tube_radius,tube_sides=tube_sides,color=color,vmin=vmin,vmax=vmax)
    
    else:
        mlab.plot3d(x_iso,y_iso,z_iso, opacity=data_alpha,tube_radius=tube_radius,tube_sides=tube_sides,color=color,vmin=vmin,vmax=vmax)

