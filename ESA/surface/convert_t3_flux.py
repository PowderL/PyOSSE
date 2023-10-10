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


#<Parameter> 

# #:notes: its unit is in kgC/m2 instead of kgCO2/m2

XNUMOL_CO2=gc.An/(1.0e-3*gc.mc)
t3_unit_factor=XNUMOL_CO2
t3_unit_factor=1.0e-4*t3_unit_factor # from kgC/m2 to molec/cm2

###

def save_t3_bf_to_netcdf(flnm_flux, \
                             flnm_map,\
                             ncflnm,\
                             dlon_flux=0.5,\
                             dlat_flux=0.5,\
                             dlon_map=1.0,\
                             dlat_map=1.0,\
                             flux_unit_convert=1.0,\
                             flux_unit='kgC/m2/s'\
                             ):
    
    

    """
    
    Convert T3 unit fluxes and regional maps into one netcdf file 
    
    Inputs:
    -----------------------------
    1. flnm_flux:<str>: file name for T3 flux perturbations
    2. flnm_map:<str>: file name for T3 region maps (reserved for future)
    
    3. ncflnm: <str>: the output file names 
    4. dlon_flux:<int>: longitude resolution for flux data 
    5. dlat_flux:<int>: latitude resolution for flux data
    6. dlon_map:<int>: longitude resolution for map 
    7. dlat_map:<int>: latitude resolution for map 
    
                         
    Notes:
    -------------------------------------------------------
    1. map will be regridded from a lower (1.0x1.0) to a higher resolution 
    
    """
    
    # S1: read T3 flux basis functions (at 0.5x0.5 degree by default)

    # #T1:  setup the grid 
    
    dlon=dlon_flux
    dlat=dlat_flux
    
    grd_lat=npy.arange(-90.0+dlat/2.0, 90.0-dlat/2.0+0.1, dlat)
    grd_lon=npy.arange(-180.0+dlon/2.0, 180.0-dlon/2.0+0.1, dlon)

    ctm_grd=grd_2d.ctm_grid_2d(grd_lon, grd_lat)
    
    nlat=ctm_grd.nlat
    nlon=ctm_grd.nlon
    
    
    print 'size (nlon, nlat):', nlon, nlat 
    
    
    # READ emission climatology from TRANSCOM file 
    
    ## ff90 fossil emission in 1990
    ## ff95 fossil emission in 1990
    ## nep  net biospheric fluxes
    ## ocean oceanic CO2 fluxes
    ## sf6  SF6 fluxes. 
    ## landunit: unified regional land GPP
    ## oceanunit: unified oceanic fluxes. 
    
    
    ff90,ff95,nep,ocean,sf6,landunit,oceanunit = rtd.read_transcom_map(flnm_flux)
    
    print 'shape landunit', npy.shape(landunit)
    
    
    # S2: read in map  

    
    nlat_std=int(180/dlat_map)
    nlon_std=int(360/dlon_map)
    
    std_lon, std_lat, std_map=fxgrd.read_std_map(nlon_std, nlat_std)

    # S3: regrid map to 0.5x0.5
    
    fine_map=fxgrd.refine_map(std_lon, std_lat,\
                                  grd_lon, grd_lat, \
                                  std_map)
    
    
    # S4: area 
    fine_area=ctm_grd.get_area()
    

    # S5: construct flux basis functions for all regions. 

    # ##c: get annual mean of the oceanunit 

    oceanunit=npy.mean(abs(oceanunit), axis=3)
    
    # ##c: get the number of land and number of ocean regions
    
    nland=npy.size(landunit[0,0,:])
    nocean=npy.size(oceanunit[0,0,:])
    

    t3_units=npy.zeros([nlon, nlat, nland+nocean+1], float)
    t3_units[:,:,1:nland+1]=landunit[:,:,0:nland]
    t3_units[:,:,nland+1:nland+nocean+1]=oceanunit[:,:,0:nocean]

   

    t3_units[:,:,nland+1:nland+nocean+1]=oceanunit[:,:,0:nocean]
    fine_map=fine_map+1
    nreg=max(fine_map.flat)
    all_flux=npy.zeros(npy.shape(fine_map), float)
    # S6: unit conversion 
    
    t3_units=flux_unit_convert*t3_units
    
    
    
    for ireg in range(nreg):
        reg_mask=npy.where(fine_map==ireg+1, 1, 0)
        all_flux=all_flux+t3_units[:,:,ireg]*reg_mask
    
    
    # S7: outputs 
    print npy.shape(grd_lon)
    
    
    nlayer=nreg
    

    xlayer=npy.arange(nlayer)
    pid=xlayer+1

    layer_info=ovar.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                  varattr=None)

    
    lon_info=ovar.io_var_cl('longitude', 'f', ['longitude'], grd_lon, \
                                varattr=None)
    
    lat_info=ovar.io_var_cl('latitude', 'f', ['latitude'], grd_lat, \
                               varattr=None)
    
       

    map_info=ovar.io_var_cl('map', 'f', ['longitude','latitude'], fine_map, \
                                varattr={"units":"1-23", "long_name":"region map", \
                                             "standard_name":"region_map"})
    
    
    # T3 flux in each layer (region)
    
    
    t3flux_info=ovar.io_var_cl('t3flux', 'f', ['longitude','latitude', 'layer'], t3_units, \
                                   varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                                "standard_name":"co2_surface_flux"})
    
    # combined T3 flux 
    
    flux_info=ovar.io_var_cl('flux', 'f', ['longitude','latitude'], all_flux, \
                                   varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                                "standard_name":"co2_surface_flux"})
    
    
    
    
    area_info=ovar.io_var_cl('area', 'f', ['longitude','latitude'], fine_area, \
                                 varattr={"units":"m2", "long_name":"grid area", \
                                              "standard_name":"area"})

    
    area_info=ovar.io_var_cl('area', 'f', ['longitude','latitude'], fine_area, \
                                 varattr={"units":"m2", "long_name":"grid area", \
                                              "standard_name":"area"})
    # parent ID
    
    pid_info=ovar.io_var_cl('pid', 'i', ['layer'], pid, \
                                varattr=None)
    
    
    ncfio.ncf_save_var(ncflnm,  [flux_info, t3flux_info, map_info, area_info, pid_info], \
                           [lon_info, lat_info, layer_info], \
                           create_new=True)
    

    
if (__name__=='__main__'):
    flnm_flux="t3_flux_input.dat"
    flnm_map='t3reg'
    ncflnm='reg_flux_t3_05x05.nc'
    
    save_t3_bf_to_netcdf(flnm_flux, flnm_map,\
                             ncflnm)
    
    
