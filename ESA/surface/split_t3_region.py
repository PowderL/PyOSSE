"""
   Wrapper for functions in divide_region.py to 
   generate 144 global regions from T3 standard regions

   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.01.12
    
   
   
   Functions:
   ==============================================================
   1. get_reg_op_lst: generate a list for action (split) for each T3 region 
   2. split_t3_region: split t3 region into  sub-regions. 
   
   
   """

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.geo_constant as gc
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ctm_grid_2d as grd_2d
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio

import  numpy as npy
import ESA.util.gen_plots as gpl
import pylab as plb
import divide_region as dreg


#<<< Parameters >>>

def_t3_varnames=['longitude', 'latitude', 'map', 'flux', 'area', 'pid']
def_t3_mapflnm='reg_flux_t3_05x05.nc'

def_t3_reg_id_lst=range(1, 24)
def_t3_reg_type_lst=[2]+11*[0]+11*[1]  # snow, land, ocean 

# sub-regions number: 1x1 for snow, 3x3 for land, 2x2 for ocean 

def_t3_nsub_lon_lst=[1]+11*[3]+11*[2]  
def_t3_nsub_lat_lst=[1]+11*[3]+11*[2]



#<<< Functions >>>

def get_reg_op_lst(reg_id_lst, nsub_lon_lst, nsub_lat_lst, \
                       split_op=dreg.divide_region_to_subregion):
    """
    construct reg_op_lst for dividing regions 

    Inputs:
    ---------------------------------------------------
    1. reg_id_lst:<list, t:int>:  region id for parent regions
    2. nsub_lon_lst:<list, t:int>: number of subregions along longitudes 
    3. nsub_lat_lst: <list, t:int>: number of subregions along latitudes
    
    
    returns:
    -------------------------------------------
    1. reg_op_lst:<list, t:reg_op_cl>: class for operations on region maps 
    

    """

    
    reg_op_lst=[]
    nreg=len(reg_id_lst)
    for ireg in range(nreg):
        reg_op=dreg.reg_op_cl(reg_id_lst[ireg],\
                                  nsub_lon_lst[ireg],\
                                  nsub_lat_lst[ireg], split_op)
        
        reg_op_lst.append(reg_op)
    
    return reg_op_lst
        


def split_t3_region(out_ncflnm, do_debug=False, **keywords):
    """
    
    split T3 regions into sub-regions
    
    Inputs:
    1. out_ncflnm:<str>: name for the name map file 
    2. do_debug:<T/F>: debug 
    3. keywords:<dict>: extra inputs 
    
    ---Reserved keywords
    --->1.map_flnm:<str>: original file 
    --->2.varnames:<list, t:str>: names for lon,lat and map in original map file 
    --->3. reg_id_lst:<list, t:id>: ids for regions to be splited 
    --->4. nsub_lon_lst:<list, t:int>: number of subregions along longitudes 
    --->5.  nsub_lat_lst: <list, t:int>: number of subregions along latitudes
    
    
    """
    
    if ('map_flnm' in keywords):
        map_flnm=keywords['map_flnm']
    else:
        map_flnm=def_t3_mapflnm
  
    if ('varnames' in keywords):
        varnames=keywords['varnames']
    else:
        varnames=def_t3_varnames
    
    # #c: read in flux
    
    lon, lat, reg_map, flux, area, pid= ncfio.ncf_read(map_flnm, varnames)
    
    # #c: check flux unit 
    
    flux_vname=varnames[3]
    flux_unit=ncfio.ncf_get_attr(map_flnm, flux_vname, 'units')
    

    if ('reg_id_lst' in keywords):
        reg_id_lst=keywords['reg_id_lst']

    else:
        reg_id_lst=def_t3_reg_id_lst

        
    if ('nsub_lon_lst' in keywords):
        nsub_lon_lst=keywords['nsub_lon_lst']

    else:
        nsub_lon_lst=def_t3_nsub_lon_lst
        
    
    if ('nsub_lat_lst' in keywords):
        nsub_lat_lst=keywords['nsub_lat_lst']

    else:
        nsub_lat_lst=def_t3_nsub_lat_lst

        
    reg_op_lst=get_reg_op_lst(reg_id_lst, nsub_lon_lst, nsub_lat_lst)
    new_map, pid_lst= dreg.divide_map(reg_map, lon, lat, reg_op_lst)
    
    if (do_debug):
        gpl.plot_map(new_map, lon, lat, use_pcolor=1)
        plb.show()
    

    tarea=npy.sum(area)
    tflux=npy.sum(flux*area)
    mean_flux=tflux/tarea
    # #c: reset low emission region to a small values 
    flux=npy.where(new_map==1, 1.e-9*mean_flux, flux)
    
    nlayer=len(pid_lst)
    xlayer=npy.arange(nlayer)
    
    pid_lst=npy.array(pid_lst)

    
    layer_info=ovar.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                  varattr=None)

    
    
    lon_info=ovar.io_var_cl('longitude', 'f', ['longitude'], lon, \
                                varattr=None)
    
    lat_info=ovar.io_var_cl('latitude', 'f', ['latitude'], lat, \
                               varattr=None)
    
    map_info=ovar.io_var_cl('map', 'f', ['longitude','latitude'],  new_map, \
                               varattr={"units":"1", "long_name":"region map", \
                                            "standard_name":"region_map"})
    
    
    pid_info=ovar.io_var_cl('pid', 'i', ['layer'],  pid_lst, \
                                varattr={"units":"1", "long_name":"parent ID", \
                                             "standard_name":"pid"})
    
    
    
    

    
    
    flux_info=ovar.io_var_cl('flux', 'f', ['longitude','latitude'], flux, \
                                 varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                              "standard_name":"co2_surface_flux"})
    
    
    area_info=ovar.io_var_cl('area', 'f', ['longitude','latitude'], area, \
                                 varattr={"units":"m2", "long_name":"grid area", \
                                              "standard_name":"area"})

    ncfio.ncf_save_var(out_ncflnm,  [map_info, flux_info, area_info, pid_info], \
                           [lon_info, lat_info, layer_info], \
                           create_new=True)   
    

#<<< TEST >>> 

if (__name__=='__main__'):
    
    out_ncflnm='reg_flux_144_05x05.nc'
    split_t3_region(out_ncflnm, do_debug=True)
    
    
