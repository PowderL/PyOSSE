""" Functions to calculate GEOS-Chem model horizontal grids 

    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.10.28
    
    Functions
    ==============================================
    1. get_model_resolution: get model resolution from sizes of longitude and latitude grid.  
    2. get_model_lat: calculate model grid latitudes

    3. get_model_lon: calculate model grid longitudes

    4. get_model_lat_edge: calculate edges of model latitude grid.  
    

    5. get_model_lon_edge: calculate edges of model longitude grid.

    6. get_wnd_lons: get window longitude range

    7. get_wnd_lats: get window latitude range
    
    8. get_area: calculate areas for grid box
    
""" 

import numpy as  npy
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm


def get_model_resolution(nx, ny):
    """
    Get_model_resolution: get model resolution from sizes of longitude and latitude grid
    
    Inputs:
    -------------------------------------------------
    1. nx:<int>: size of longitude
    2. ny:<int>: size of latitude
    
    Returns:
    --------------------------------------
    1. model_res:<str>: model resolution
    
    
    """

    model_res=""
    if (nx==72 and ny== 46):
        model_res='4x5'
        
    elif (nx==144 and ny== 91):
        model_res='2x25'
    elif (nx==288 and ny== 181):
        model_res='1x125'
    elif (nx==360 and ny== 181):
        model_res='1x1'
    elif (nx==540 and ny== 361):
        model_res='0.5x0.666'
    else:
        msg='Not a valid GEOS model grid' 
        msm.show_err_msg(msg)
        model_res=""
    
    
    return model_res

def get_model_lat(model_res=None,nlat=0):
    """calculate model grid latitudes 
    Inputs:
        -------------------------------------
        1. model_res: <string>: model_resolution
        2. nlat:     <integer>: number of latitudes
        
        Returns:
        1. lats: <array, (ny,)>: model grid centre latitudes
       
    """

    if (model_res==None):
        ny=nlat
    
    elif (model_res=='4x5'):
        ny=46
    elif (model_res=='2x25'):
        ny=91
    elif (model_res=='1x125'):
        ny=181
    elif (model_res=='1x1'):
        ny=181
    elif (model_res=='0.5x0.666'):
        ny=361
    elif (model_res=='05x0666'):
        ny=361
    else:
        ny=nlat

    
    dy=180.0/(ny-1)
    fmid  = 0.5*(ny + 1 )
    
    lat=npy.zeros(ny, float)

    #  print *, 'assign value'
    
    for j in range(0, ny):
        lat[j]=dy*(j+1-fmid)
        
    # re-locate start/end points
        
    if (model_res=='4x5'):
        lat[0] = -89.0
        lat[-1] = 89.0
    
    elif (model_res=='2x25'):

        lat[0] = -89.5
        lat[-1] = 89.5
    elif (model_res=='1x125'):
        lat[0] = -89.75
        lat[-1] = 89.75
        
    elif (model_res=='1x1'):
        lat[0] = -89.75
        lat[-1] = 89.75
    elif (model_res=='0.5x0.666'):
        lat[0] = -89.875
        lat[-1] = 89.875
    elif (model_res=='05x0666'):
        lat[0] = -89.875
        lat[-1] = 89.875
    else:
        pass
    
    return lat

def get_model_lon(model_res=None,nlon=0):
    """calculate model grid longitudes 
    Inputs:
        -------------------------------------
        1. model_res: <string>: model_resolution
        2. nlon:     <integer>: number of longitudes
        
        Returns:
        1. lon: <array, (nx,)>: model grid center longitudes
       
    """

    if (model_res==None):
        nx=nlon
    
    elif (model_res=='4x5'):
        nx=72
    elif (model_res=='2x25'):
        nx=144
    elif (model_res=='1x125'):
        nx=288
    elif (model_res=='1x1'):
        nx=360
    elif (model_res=='0.5x0.666'):
        nx=540
    elif (model_res=='05x0666'):
        nx=540
    else:
        nx=nlon
        
    dx=360.0/(nx)
    lon=npy.zeros(nx, float)
    lon[0]=-180.0
    
    
    for j in range(1, nx):
        lon[j]=lon[j-1]+dx
    
    
    return lon    
                                            

    

def get_model_lat_edge(model_res=None,nlat=0):
    """calculate edges of model latitude grid.
    
    Inputs:
    -------------------------------------------------
    1. model_res: <string>: model_resolution
    2. nlat:     <integer>: number of latitudes
    
    Returns:
    ---------------------------------------------------
    1. lat: <array, (ny+1,)>: model grid edge latitudes
    
    """
    
    #S1 determine size
    
    if (model_res==None):
        ny=nlat
    elif (model_res=='4x5'):
        ny=46
    elif (model_res=='2x25'):
        ny=91
    elif (model_res=='1x125'):
        ny=181
    elif (model_res=='1x1'):
        ny=181
    elif (model_res=='0.5x0.666'):
        ny= 361
    elif (model_res=='05x0666'):
        ny= 361
    else:
        ny=nlat
    
    dy=180.0/(ny-1)
    fmid  = 0.5*(ny + 2 )
    
    lat=npy.zeros(ny+1, float)

    #S2 assign values 
    
    for j in range(0, ny):
        lat[j]=dy*(j+1-fmid)
        
    #S3 re-locate polar points
        
    lat[0]=-90.0
    lat[-1]=90.0
    
    
    return lat

def get_model_lon_edge(model_res=None,nlon=0):

    """calculate model grid longitudes 
    Inputs:
    ----------------------------------------------------
    1. model_res: <string>: model_resolution
    2. nlon:     <integer>: number of longitudes
    
    Returns:
    ------------------------------------------------
    1. lon: <array, (nx+1,)>: model grid edge longitudes
       
    """
    #S1 determine size
    
    if (model_res==None):
        nx=nlon
    elif (model_res=='4x5'):
        nx=72
    elif (model_res=='2x25'):
        nx=144
    elif (model_res=='1x125'):
        nx=288
    elif (model_res=='1x1'):
        nx=360
    elif (model_res=='0.5x0.666'):
        nx=540
    elif (model_res=='05x0666'):
        nx=540
    
    else:
        nx=nlon
    
    #S2 assign values
    
    dx=360.0/(nx)
    lon=npy.zeros(nx+1, float)
    lon[0]=-180.0-0.5*dx
    
    for j in range(1, nx+1):
        lon[j]=lon[j-1]+dx
    
    
    return lon    


def get_wnd_lons(glb_lons, ix, ifirst,  centre180):

    """get window longitude range 
    
    Inputs:
    -------------------------------------
    1. glb_lons: <array, (nx,)>: model grid longitudes
    2. ix:     <integer>: window x (longitude) size
    3. ifirst:  <integer>: window shift
    4. center180:<integer>: 0==model grid starting from -180.0
    
    Returns:
    -----------------------------------------------
    1. lons: <array, (ix)>: window longitudes
       
    """
    #S1 get the slice
    lons=glb_lons[ifirst-1:ifirst-1+ix]
    #S2 if required, re-arrange it so that '180 is at center'
    
    if (centre180==1):
        rlons=npy.array(lons)
        idx_lst=npy.where(rlons<0.0)
        nidx=npy.size(idx_lst)
        if (nidx>0):
            rlons[0:ix-nidx]=lons[nidx:]
            rlons[nidx:]=lons[0:ix-nidx]+360.0
            
        lons=rlons
    return lons


def get_wnd_lats(glb_lats, jx, jfirst):
    """get window latitude range 

    Inputs:
    -------------------------------------
    1. glb_lats: <array, (ny,)>: model grid latitudes
    2. jx:     <integer>: window y (latitude) size
    3. jfirst:  <integer>: window y shift
    
    Returns:
    ----------------------------------------------
    1. lats: <array, (jx)>: window latitudes
       
    """
    #S1 get the slice
    
    lats=glb_lats[jfirst-1:jfirst-1+jx]
    return lats

def get_area(lon_edge, lat_edge):
    
    """
    calculate areas for grid box (in m2)
    
    Inputs:
    --------------------------------------
    1. lon_edge:<array, (nx+1,)>: edges of the X (longitude) grid boxes
    2. lat_edge:<array, (ny+1,)>: edges of the Y (latitude) grid boxes
    
    Returns:
    ------------------------------
    1. area:<array, (nx, ny)>: area of grid box in m^2
      
    """
    
    deg2rad=npy.pi/180.0
    # S1:  get grid egdes 
    rlon=lon_edge
    rlat=lat_edge
      
        
    nx=npy.size(rlon)
    ny=npy.size(rlat)
    
      
    # S2:  dx (longitude) in RAD
    
    dx=npy.zeros(nx-1, float)
    dx=deg2rad*dx
    

    for ix in range(nx-1):
        # degree to radius
        dx[ix]=rlon[ix+1]-rlon[ix]
          
    dx=deg2rad*dx
    

      
    # S3:  R^2*sin(y); y=latitude
      
    a2=gc.re*gc.re*npy.sin(rlat*deg2rad) 
    
    # S4:  area=R2*dx*d(sin(y))
        
    area=npy.zeros([nx-1, ny-1], float)
    
    for ilat in range(ny-1):
        area[:, ilat]=(a2[ilat+1]-a2[ilat])*dx
          
    return area
