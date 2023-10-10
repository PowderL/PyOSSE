""" 
 Functions and interfaces to calculate model pressures
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2012.09.12
 

    Functions
    =================================================
    1. get_mod_pres:  get pressure at model levels from given surface pressure
    2. get_mod_pres_edge: get pressure at model level edges from given surface pressure
    3. get_geos5_ap_bp:   get pressure coefficients for GEOS5 grid 
    4. get_geos4_ap_bp:   get pressure coefficients for GEOS4 grid 
    
    """

import sigma_pres_mod as pm
import numpy as npy

#<<< FUNCTIONS >>> 

def get_mod_pres(ps, ap, bp):
    """ get pressures at model levels from given surface pressure
    
    Inputs:
    -------------------------------------
    1. ps:<array, (nx,) or (nx, ny)>: surface pressure in hPa
    2. ap:<array, (nz,)>:  surface pressure coefficients:a 
    3. bp:<array, (nz,)>:  surface pressure coefficients: b
    
    Returns:
    -------------------------------------------------------
    1. pres<array, <nx, (ny), nz>: pressure at model level center 
    

    """

    dims=npy.shape(ps)
    ndim=len(dims)
    print npy.shape(ap)
    print npy.shape(bp)
    print 'ps', npy.shape(ps)
    
    

    pres, pres_edge=pm.get_model_pres(ps, ap, bp)
    if (ndim==0): # ps is a single value 
        mod_pres=pres[0,0,:]
    elif (ndim==1): # ps is a 1d-array
        mod_pres=pres[:,0,:]
    elif (ndim==2): # ps is a 2d-array
        mod_pres=pres
        
    else:
        print 'Warning! Dimensions of PS are too many', dims
        mod_pres=None

    return mod_pres

def get_mod_pres_edge(ps, ap, bp):
    """ get pressure at model level edges from given surface pressure
    
    Inputs:
    -------------------------------------
    1. ps:<array, (nx,) or (nx, ny)>: surface pressure in hPa
    2. ap:<array, (nz,)>:  surface pressure coefficients:a 
    3. bp:<array, (nz,)>:  surface pressure coefficients: b
    
    Returns:
    1. pres<array, <nx, (ny), nz+1>: pressure at model level edge 
    
    
    """

    dims=npy.shape(ps)
    ndim=len(dims)
    pres, pres_edge=pm.get_model_pres(ps, ap, bp)
    if (ndim==0): # ps is a single value 
        mod_pres=pres_edge[0,0,:]
    elif (ndim==1): # ps is a 1d-array
        mod_pres=pres_edge[:,0,:]
    elif (ndim==2): # ps is a 2d-array
        mod_pres=pres_edge
    
    else:
        print 'Warning! Dimensions of PS are too many', dims
        mod_pres=None

    return mod_pres

        
    

    
def get_geos5_ap_bp(reduced_grid=1):
    """
    get pressure coefficients for GEOS5 grid 
    
    Inputs
    -------------------------------------
    1. reduced_grid:<integer>: 1==use reduced grid 
    
    Returns:
    --------------------------------------------
    1.ap, bp: <array, (nz,)>: ap and bp coefficient 
    
    """
    
    nlvl, ap, bp=pm.get_gc_ap_bp(5, reduced_grid)
    
    ap=ap[0:nlvl]
    bp=bp[0:nlvl]
    return ap, bp


def get_geos4_ap_bp(reduced_grid=1):
    """
    get pressure coefficients for GEOS4 grid 
    
    Inputs
    -------------------------------------
    1. reduced_grid:<integer>: 1==use reduced grid 
    
    Returns:
    -------------------------------------
    1.ap, bp: <array, (nz,)>: ap and bp coefficient 
    
    """
   
    nlvl, ap, bp=pm.get_gc_ap_bp(4, reduced_grid)
    ap=ap[0:nlvl]
    bp=bp[0:nlvl]
    return ap, bp



#<<< TEST >>>


if (__name__=='__main__'):
    
    #1 get ap, bp

    print '1. Get ap & bp for GEOS5 grid' 
    ap, bp=get_geos5_ap_bp(reduced_grid=1)
    
    print 'Shape of ap and bp', npy.shape(ap), npy.shape(bp)

    #2 set surface pressure as 1000 hPa and 850 hPa for two points
    
    print  '2. Set surface pressure as 1000 hPa and 850 hPa for two points'
    
    ps=npy.zeros(2, float)
    ps[0]=1000.0
    ps[1]=850.0
    
    #3 read in model pressures at these two points:


    print '3. Read in model level center pressures at these two points'

    mod_pres=get_mod_pres(ps, ap, bp)

    print 'Shape:', npy.shape(mod_pres)
    print 'Pressure at level 0:', mod_pres[:,0]
    print 'Pressure at  top:', mod_pres[:,-1]

    
    print '4. Read in model level edge pressures at these two points'
    
    mod_pres=get_mod_pres_edge(ps, ap, bp)
    
    print 'Shape:', npy.shape(mod_pres)
    print 'Pressure at level 0:', mod_pres[:,0]
    print 'Pressure at  top:', mod_pres[:,-1]
    
