


def map3d_surface(mode,xdata,ydata,zdata,scalardata,vmin,vmax,data_cmap='blue-red',data_alpha=1,topo=None,topo_limits=None,zscale=500.,topo_vmin=None,topo_vmax=None,topo_cmap='bone',topo_cmap_reverse=False,land_constant=False,land_color=(0.7,0.7,0.7),set_view=None):
    """
    mode = (string) coordinate system of 3D projection. Options are 'rectangle' (default), 'spherical' or 'cylindrical'
    xdata = (1D numpy array) longitude values for data array
    ydata = (1D numpy array) latitude values for data array
    zdata = (1D numpy array) depth values for data array
    scalardata = (2D numpy array) scalar field to plot colors on surface
    vmin = (float) colorbar minimum for data
    vmax = (float) colorbar maximum for data
    data_cmap = colormap for data surface, default is blue-red
    data_alpha = (float or int) opacity for data surface from 0 to 1, default is 1
    topo = optional; input topography file, default is etopo 30 
    topo_limits = optional; longitude and latitude limits for 3d topography plot [lon min, lon max, lat min, lat max], longitudes range -180 to 180, latitude -90 to 90, default is entire globe
    zscale = optional; change vertical scaling for plotting, default is 500
    topo_cmap = optional; default is bone 
    topo_cmap_reverse = optional; reverse topography colormap, default is false
    set_view = optional; set the mayavi camera angle with input [azimuth, elevation, distance, focal point], default is 
    """
    
    #import dependent modules
    import sys
    import numpy as np
    from mayavi import mlab

    #load topo data
    data = np.load('etopo1_30min.npz')
    xraw = data['x']
    yraw = data['y']
    zraw = np.swapaxes(data['z'][:,:],0,1)
    zraw[zraw>0]=0.
    phi = (yraw[:]*np.pi*2)/360.+np.pi/2.
    theta = (xraw[:]*np.pi*2)/360.
    c = zraw
    theta=np.append(theta,theta[0])
    c = np.concatenate((c,np.expand_dims(c[0,:],axis=0)),axis=0)

    if topo_limits is None:
        #set default limits (phi=latitude, theta=longitude)
        phi_1 = -90
        phi_2 = 90
        theta_1 = -180
        theta_2 = 180
    
    else:
        phi_1 = topo_limits[2]
        phi_2 = topo_limits[3]
        theta_1 = topo_limits[0]
        theta_2 = topo_limits[1]
    
    #find indices of topo limits
    phi_deg = data['y']
    theta_deg = data['x']
    phi_ind1 = np.argmin(np.abs(phi_deg-phi_1))
    phi_ind2 = np.argmin(np.abs(phi_deg-phi_2))
    theta_ind1 = np.argmin(np.abs(theta_deg-theta_1))
    theta_ind2 = np.argmin(np.abs(theta_deg-theta_2))

    #restrict topo extent
    phi=phi[phi_ind1:phi_ind2]
    theta=theta[theta_ind1:theta_ind2]
    c = c[theta_ind1:theta_ind2:,phi_ind1:phi_ind2]
    phi, theta = np.meshgrid(phi,theta)
   

    #prep data grid
    phi_iso, theta_iso = np.meshgrid(((ydata*np.pi*2)/360.)+np.pi/2.,(xdata*np.pi*2)/360.)
 
    if mode is 'sphere':
        x = np.sin(phi) * np.cos(theta[::-1]) * (1 + c/zscale)
        y = np.sin(phi) * np.sin(theta[::-1]) * (1 + c/zscale)
        z = np.cos(phi) * (1 + c/30000.)
    
        x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1]) * (1 -zdata/zscale)
        y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1]) * (1 -zdata/zscale)
        z_iso = np.cos(phi_iso) * (1 -zdata/zscale)
    
    elif mode is 'cylinder':
        x = np.sin(phi) * np.cos(theta[::-1])
        y = np.sin(phi) * np.sin(theta[::-1])
        z = c/zscale

        x_iso = np.sin(phi_iso) * np.cos(theta_iso[::-1])
        y_iso = np.sin(phi_iso) * np.sin(theta_iso[::-1])
        z_iso = -zdata/zscale
    
    elif mode is 'rectangle':
        y, x = np.meshgrid(phi_deg[phi_ind1:phi_ind2],theta_deg[theta_ind1:theta_ind2])
        z = c/zscale
    
        y_iso,z_iso = np.meshgrid(ydata,zdata)
        x_iso,z_iso = np.meshgrid(xdata,zdata) 
        z_iso =-z_iso/zscale

    if topo_vmin is None:
        tvmin = -6000
    else:
        tvmin = topo_vmin
    if topo_vmax is None:
        tvmax = 0
    else:
        tvmax = topo_vmax
    
    #make figure
    mlab.figure(size = (1024,768),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
    mlab.clf()
    # Plot Bathymetry mesh
    m = mlab.mesh(x, y, z, colormap=topo_cmap,vmin=tvmin/zscale,vmax=tvmax/zscale)
    
    #optional: reverse bathymetry colormap
    if topo_cmap_reverse is True:
        lut = m.module_manager.scalar_lut_manager.lut.table.to_array() 
        ilut=lut[::-1]
        m.module_manager.scalar_lut_manager.lut.table = ilut

    #optional: plot constant color on land
    if land_constant is True:
        sl = mlab.mesh(x, y, z,mask = z<0,color =land_color)


    #plot data surface
    m = mlab.mesh(x_iso, y_iso, z_iso,scalars=scalardata,colormap=data_cmap,vmin =vmin,vmax=vmax,opacity=data_alpha)
    m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]
    
    #optional: change mayavi camera settings
    if set_view is None:
        mlab.view(distance = 'auto')
    else:
        mlab.view(azimuth = set_view[0], elevation = set_view[1], distance = set_view[2], focalpoint = set_view[3])
    mlab.show()

    return

