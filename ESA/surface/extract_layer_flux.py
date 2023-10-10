"""
   Functions to regrid and store regional fluxes (basis functions) 
   into multiple layer maps at coarse resolution 
   
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.01.21
   
   Classes:
   =============================================================
   1. region_flux_cl: class for mult-layer regional fluxes 
   
   
   Functions:
   ==============================================================
   1.  extract_region_flux: extract regional fluxes from 0.5x0.5 to 
   multi-layer maps 
   
"""

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.geo_constant as gc
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ctm_grid_2d as grd_2d
import ESA.atmosphere.gc_grid_2d as gcgrd_2d

import ESA.util.gen_plots as gpl
import divide_region as dreg
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.flux_regrid as fgrd

import  numpy as npy
import pylab as plb

#<<< parameters >>>
# ##c@ names of variables to be read from flux file 

def_flux_varnames=['longitude', 'latitude', 'map', 'flux', 'area', 'pid']

# ##c: file name for single-layer flux 

def_fluxflnm='reg_flux_144_05x05.nc'

# ##c: region number (not in use)

def_nreg=144

# ##c: region number (not in use)

def_reg_id_lst=range(1, def_nreg+1)

#<<< Classes >>> 

class region_flux_cl:
    """
    class for mult-layer regional fluxes 
    
    """
    def __init__(self, filename, fopen, \
                     fread, \
                     fclose, \
                     fio_keywords={},\
                     **keywords):
        
        self.lon=None
        self.lat=None
        self.reg_map=None
        self.reg_flux=None
        self.pid=None
        
        self.filename=filename
        self.fopen=fopen
        self.fread=fread
        self.fclose=fcose
        self.fwrite=fwrite


def extract_region_flux(flux_flnm, \
                            out_ncflnm, \
                            ctm_grd, \
                            do_debug=False, **keywords):
    """
    Extract regional fluxes to a multi-layer map 
    
    Inputs:
    --------------------------------------------------------
    1. flux_flnm:<str>: name for single-layer regional map and flux
    2. out_ncflnm:<str>: name for the name map file 
    3. ctm_grd:<ctm_grd_2d>: class for output model grid
    4. do_debug:<T/F>: debug 
    5. keywords:<dict>: extra inputs 
    
    ---Reserved keywords
    --->varnames:<list, t:str>: names for lon,lat and map in original map file 
    --->reg_id_lst:<list, t:id>: ids for regions to be splited 
    --->nsub_lon_lst:<list, t:int>: number of subregions along longitudes 
    --->nsub_lat_lst: <list, t:int>: number of subregions along latitudes
    
        
    """
    
    
    # S1: read in flux from file 
    
    if ('varnames' in keywords):
        varnames=keywords['varnames']
    else:
        varnames=def_flux_varnames
    
    
    lon, lat, reg_map, flux, area, pid= ncfio.ncf_read(flux_flnm, varnames)
    
    # #T1: convert from per area to per box 
    # #c:from m2 to cm2
    
    flux_vname=varnames[3]
    flux_unit=ncfio.ncf_get_attr(flux_flnm, flux_vname, 'units')

    print 'max flux', max(flux.flat)
    
    if ('/cm2' in flux_unit):
        flux=1.0e4*area*flux

    elif ('/m2' in flux_unit):
        flux=area*flux


    print 'max flux after', max(flux.flat)
    
    
    nlayer=int(npy.max(reg_map.flat))
    
    # #T2: setup output grid 
 
    
    grd_lon=ctm_grd.get_lon()
    grd_lat=ctm_grd.get_lat()
    
    grd_area=ctm_grd.get_area()
    # print npy.shape(grd_area)
    
    # S2: extract regional map and flux to multiple-layer data
    
    # print npy.max(reg_map), npy.min(reg_map)
    # gpl.plot_map(reg_map, lon, lat, use_pcolor=1)
    # flux_0=npy.where(reg_map==1, flux, 0.0)
    # print 'max flux_0', npy.max(flux_0.flat)
    
    # plb.show()
    
    
    grd_map, grd_flux= fgrd.extract_flux_to_multi_layer(lon, lat, \
                                                            grd_lon,grd_lat,\
                                                            reg_map, nlayer, flux)
    
    print npy.shape(grd_map)
    
    # S3: convert flux from per grid box to per m2

    
    usd_idx=npy.where(grd_area>0)
    un_usd_idx=npy.where(grd_area==0)
    if ('/cm2' in flux_unit):
        area_factor=1.0e-4
    else:
        area_factor=1.0
        
    for ilayer in range(nlayer):
        tx=grd_flux[:,:,ilayer]
        if ('/cm2' in flux_unit):
            tx[usd_idx]=tx[usd_idx]/grd_area[usd_idx]
        elif ('/m2' in flux_unit):
            tx[usd_idx]=tx[usd_idx]/grd_area[usd_idx]
        
        tx[un_usd_idx]=0.0
        tx=area_factor*tx
        
        grd_flux[:,:,ilayer]=tx[:,:]
    
        
    
    
    if (do_debug):
        plb.subplot(2,1,1)
        gpl.plot_map(grd_map[:,:,9], grd_lon, grd_lat, use_pcolor=1)
        
        plb.subplot(2,1,2)
        gpl.plot_map(grd_map[:,:,10], grd_lon, grd_lat, use_pcolor=1)
        
        plb.figure(2)
        
        plb.subplot(2,1,1)
        gpl.plot_map(grd_flux[:,:,9], grd_lon, grd_lat, use_pcolor=1)
        
        plb.subplot(2,1,2)
        gpl.plot_map(grd_flux[:,:,10], grd_lon, grd_lat, use_pcolor=1)
        
        
        plb.show()
    
    

    
    # S4: save data
        
    
    # ##c: layer
        
    
    xlayer=npy.arange(nlayer)
    layer_info=ovar.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                  varattr=None)

    # ##c: lon
        
    lon_info=ovar.io_var_cl('longitude', 'f', ['longitude'], grd_lon, \
                                varattr=None)
    # ##c: lat

    lat_info=ovar.io_var_cl('latitude', 'f', ['latitude'], grd_lat, \
                               varattr=None)
    # ##c: map
    map_info=ovar.io_var_cl('map', 'f', ['longitude','latitude', 'layer'],  grd_map, \
                                varattr={"units":"1", "long_name":"region map", \
                                             "standard_name":"region_map"})
    
    # ##c: flux
    
    flux_info=ovar.io_var_cl('flux', 'f', ['longitude','latitude', 'layer'], grd_flux, \
                                 varattr={"units":flux_unit, "long_name":"CO2 surface flux", \
                                              "standard_name":"co2_surface_flux"})
    
    # ##c: pid
    
    pid_info=ovar.io_var_cl('pid', 'i', ['layer'],  pid, \
                                varattr={"units":"1", "long_name":"parent ID", \
                                             "standard_name":"pid"})
    
    
    # ##c: area
    area_info=ovar.io_var_cl('area', 'f', ['longitude','latitude'], grd_area, \
                                 varattr={"units":"m2", "long_name":"grid area", \
                                              "standard_name":"area"})
    
    ncfio.ncf_save_var(out_ncflnm,  [map_info, flux_info, area_info, pid_info], \
                           [lon_info, lat_info, layer_info], \
                           create_new=True)   
    

# <<< TEST >>> 

if (__name__=='__main__'):
    # build up a 4x5
    ctm_grd=gcgrd_2d.gc_grid_2d(0, 0, mod_res='4x5')
    flux_flnm='reg_flux_144_05x05.nc'
    out_ncflnm='reg_144_ml.4x5.nc'
    
    extract_region_flux(flux_flnm, out_ncflnm, ctm_grd, \
                            do_debug=True)
    
