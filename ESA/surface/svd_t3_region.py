"""
    Functions for doing composition of regional Basis functions (i.e, error)
    

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.24
    
    
    Functions:
    ===========================================
    svd_t3_region: do svd with T3 regions
    
    
"""

import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.gen_plots as gpl
import svd_region_flux as svd_flux_m
import ESA.util.flux_regrid as fgrd
import pylab as plb
import numpy as npy

def_svd_t3_flux_flnm='reg_23_ml.4x5.nc'

# names to read multiple layer flux map
def_svd_t3_varnames=['longitude', 'latitude', 'map', 'flux', 'area', 'pid']

# land/ocean  error correlation length 

def_svd_land_cor_len=800 # km 
def_svd_ocean_cor_len=2000 
# number of land regions 

def_svd_nland=12

def svd_t3_regions(flux_flnm, out_ncflnm,\
                       threshold, \
                       scale, \
                       dist_threshold, \
                       maxsel, \
                       nland=def_svd_nland, \
                       land_cor_length=def_svd_land_cor_len, \
                       ocean_cor_length=def_svd_ocean_cor_len, \
                       **keywords):
    
    """
    
    SVD error map over region to several  sub-regions defined by (nsb_lon, nsb_lat)
    
    Inputs:
   
    --------------------------------------
    1. grd_lon:<array, (nlon)>: longitude
    2. grd_lat:<array, (nlat)>: latitude
    3. reg_map:<array, (nlon, nlat)>: ratio  map defining one region. 
    4. reg_err:<array, (nlon, nlat)>: err (or BF) defined over given region
    5. cor_length <float>: length of the correlation 
    6. threshold:<float>: any grid boxes below this threshold will be treated as one single region 
    7. scale:<float>: scaling factor applied to this region 
    8. dist_threshold:<float>: the max distance/cor_length where correlation will be considered 
    9. maxsel:<int>: the maximum number to be keep
    
    Returns:
    --------------------------------------------
    1. flux:<array, (nlon, nlat, neig)>: new regions masked with sub-region number 
    2. new_map_lst:<list, t:array>: the layer map with each layer for each sub-regions

    """
     
    # S1: read in T3 flux and map 
    
    varnames=def_svd_t3_varnames
    
    if ('varnames' in keywords):
        varnames=keywords['varnames']

        
        
    lon, lat, reg_map, flux, area, pid= ncfio.ncf_read(flux_flnm, varnames)
    
    # #c: check flux unit 
    
    flux_vname=varnames[3]
    flux_unit=ncfio.ncf_get_attr(flux_flnm, flux_vname, 'units')
    
    nx, ny, nreg=npy.shape(flux)
    
    sel_reg_lst=range(1, nreg+1)
    
    # regions to be done svd
    
    if ('sel_reg_lst' in keywords):
        sel_reg_lst=sel_reg_lst['sel_reg_lst']

    all_svd_flux=[]
    all_svd_map=[]
    all_pid=[]
    
    parent_map=list()
    parent_flux=list()
    parent_id=list()
    
    for isel in sel_reg_lst:
            
        sel_reg_map=reg_map[:,:,isel-1]
        sel_reg_flux=flux[:,:,isel-1]
        if ('/cm2' in flux_unit):
            sel_reg_flux=1.0e4*sel_reg_flux*area
        elif ('/m2' in flux_unit):
            sel_reg_flux=sel_reg_flux*area
        
        sel_pid=pid[isel-1]
        
        parent_map=parent_map+[sel_reg_map]
        parent_flux=parent_flux+[sel_reg_flux]
        parent_id=parent_id+[sel_pid]
        
                
        if (isel==1):
            # don't do SVD for low emission region
            
            
            nflux=1
            
            all_svd_flux=all_svd_flux+[sel_reg_flux]
            all_svd_map=all_svd_map+[sel_reg_map]
            all_pid=all_pid+[sel_pid]

        
        else:
            
            if (isel>nland):
                # ocean region
                
                cor_length=ocean_cor_length
            else:
                # land region
                cor_length=land_cor_length
                
                
            out_dict=svd_flux_m.do_region_err_svd(lon, lat, \
                                                      sel_reg_map, sel_reg_flux, \
                                                      cor_length, \
                                                      threshold=threshold, \
                                                      scale=scale, \
                                                      dist_threshold=dist_threshold, \
                                                      maxsel=maxsel)
            
            
            
                
            svd_flux=out_dict['flux']
            svd_map=out_dict['map']
            nflux=len(svd_flux)
            
            
            # #c: add element to list
            
            all_svd_flux=all_svd_flux+svd_flux
            
            all_svd_map=all_svd_map+svd_map
            
            all_pid=all_pid+nflux*[sel_pid]
            
    
    
    # S4: convert to array
    
    nlayer=len(all_svd_map)
    
    # elements 

    all_svd_flux=npy.array(all_svd_flux)
    all_svd_flux=npy.transpose(all_svd_flux, [1,2,0])
    
    all_svd_map=npy.array(all_svd_map)
    all_svd_map=npy.transpose(all_svd_map, [1,2,0])
    print all_pid
    
    all_pid=npy.array(all_pid)
    
    # parents 
    
    parent_map=npy.array(parent_map)
    parent_map=npy.transpose(parent_map, [1,2,0])
    
    parent_flux=npy.array(parent_flux)
    parent_flux=npy.transpose(parent_flux, [1,2,0])
    
    parent_id=npy.array(parent_id)
    


    # S5: save to netcdf file 

    nparent_layer=npy.size(parent_id)

    
    xlayer=npy.arange(nlayer)
    
    xparent_layer=npy.arange(nparent_layer)
    
    layer_info=ovar.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                  varattr=None)

    parent_layer_info=ovar.io_var_cl('parent_layer', 'i', ['parent_layer'], xparent_layer, \
                                         varattr=None)
    
    
    
    lon_info=ovar.io_var_cl('longitude', 'f', ['longitude'], lon, \
                                varattr=None)
    
    lat_info=ovar.io_var_cl('latitude', 'f', ['latitude'], lat, \
                               varattr=None)
    
    map_info=ovar.io_var_cl('map', 'f', ['longitude','latitude', 'layer'],  all_svd_map, \
                                varattr={"units":"1", "long_name":"region map", \
                                             "standard_name":"region_map"})
    print 'all_pid', all_pid
    
    
    pid_info=ovar.io_var_cl('pid', 'i', ['layer'],  all_pid, \
                                varattr={"units":"1", "long_name":"parent ID", \
                                             "standard_name":"pid"})
   
    parent_id_info=ovar.io_var_cl('parent_id', 'i', ['parent_layer'],  parent_id, \
                                      varattr={"units":"1", "long_name":"parent ID", \
                                                   "standard_name":"pid"})
   
    if ('/cm2' in flux_unit):
        flux_unit=flux_unit.replace('/cm2', '')
        
    elif ('/m2' in flux_unit):
        flux_unit=flux_unit.replace('/m2', '')
    
    
    flux_info=ovar.io_var_cl('flux', 'f', ['longitude','latitude', 'layer'], all_svd_flux, \
                                 varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                              "standard_name":"co2_surface_flux"})
    
    
    area_info=ovar.io_var_cl('area', 'f', ['longitude','latitude'], area, \
                                 varattr={"units":"m2", "long_name":"grid area", \
                                              "standard_name":"area"})
    

    parent_flux_info=ovar.io_var_cl('parent_flux', 'f', ['longitude','latitude', 'parent_layer'], parent_flux, \
                                 varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                              "standard_name":"co2_surface_flux"})
    
    parent_map_info=ovar.io_var_cl('parent_map', 'f', ['longitude','latitude', 'parent_layer'], parent_map, \
                                       varattr={"units":"1", "long_name":"region map", \
                                             "standard_name":"region_map"})

    
    
    ncfio.ncf_save_var(out_ncflnm,  [map_info, flux_info, area_info, \
                                         pid_info, parent_flux_info, parent_map_info, parent_id_info], \
                           [lon_info, lat_info, layer_info, parent_layer_info], \
                           create_new=True)   
    
    

#<<< TEST >>> 

if (__name__=='__main__'):
    
    out_ncflnm='reg_flux_svd.nc'
    flux_flnm='reg_23_ml.4x5.nc'
    threshold=0.3
    scale=1.0
    dist_threshold=4.0
    maxsel=30
    
    svd_t3_regions(flux_flnm, out_ncflnm,\
                       threshold, \
                       scale, \
                       dist_threshold, \
                       maxsel, \
                       nland=def_svd_nland, \
                       land_cor_length=def_svd_land_cor_len, \
                       ocean_cor_length=def_svd_ocean_cor_len)
    
    
