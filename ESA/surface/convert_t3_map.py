"""
Functions to divide the regions into smaller areas 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.06

    Functions:
    ==============================================================
    1. save_t3_bf_to_netcdf: Convert T3 unit fluxes and regional maps into one netcdf file 
     
"""


import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.geo_constant as gc
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ctm_grid_2d as grd_2d
import ESA.util.read_data as rtd
import ESA.util.flux_regrid as fxgrd
    

import numpy as npy

###

def save_t3_map_to_netcdf(   flnm_map,\
                                 ncflnm,\
                                 dlon_map=1.0, \
                                 dlat_map=1.0\
                                 ):
    
    

    """
    
    Convert T3 regional maps into one netcdf file 
    
    Inputs:
    -----------------------------
    1. flnm_flux:<str>: file name for T3 flux perturbations
    2. flnm_map:<str>: file name for T3 region maps (reserved for future)
    
    3. dlon_map:<int>: longitude resolution for map 
    4. dlat_map:<int>: latitude resolution for map 
    
                         
    Notes:
    -------------------------------------------------------
    1. map is at  (1.0x1.0) to a higher resolution 
    
    """
    
    # S1: read map 

    # #T1:  setup the grid 
    
    nlat_std=int(180/dlat_map)
    nlon_std=int(360/dlon_map)
    
    std_lon, std_lat, std_map=fxgrd.read_std_map(nlon_std, nlat_std)

    lon_info=ovar.io_var_cl('longitude', 'f', ['longitude'], std_lon, \
                                varattr=None)
    
    lat_info=ovar.io_var_cl('latitude', 'f', ['latitude'], std_lat, \
                               varattr=None)
    
       

    map_info=ovar.io_var_cl('map', 'f', ['longitude','latitude'], std_map, \
                                varattr={"units":"1-23", "long_name":"region map", \
                                             "standard_name":"region_map"})
    
    
    
    ncfio.ncf_save_var(ncflnm,  [map_info], \
                           [lon_info, lat_info], \
                           create_new=True)
    

    
if (__name__=='__main__'):
    flnm_map='t3reg'
    ncflnm='t3_map.nc'
    
    save_t3_map_to_netcdf(flnm_map,\
                             ncflnm)
    
    
