"""
   Test for constructing emission ensembles for 144 regions from a single-layer flux 
mape for 23 T3 regions. 

 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.01.22
   
   
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

import split_t3_region as t3dvd
import extract_layer_flux as mpflux
import gen_ensemble_flux as gems
import flux_ensemble_def as cfg


# S1: split regions from 23 to 144

print '>1. split regions from 23 to 144'

t3_mapflnm='t3_reg_flux_05x05.nc'
t3_varnames=['longitude', 'latitude', 'map', 'flux', 'area']
t3_reg_id_lst=range(1, 24)
t3_reg_type_lst=[2]+11*[0]+11*[1]  # snow, land, ocean 

# sub-regions number: 1x1 for snow, 3x3 for land, 2x2 for ocean 

t3_nsub_lon_lst=[1]+11*[3]+11*[2]  
t3_nsub_lat_lst=[1]+11*[3]+11*[2]
flux_ncflnm='reg_144_05x05.nc'



t3dvd.split_t3_region(flux_ncflnm, map_flnm=t3_mapflnm, \
                          varnames=t3_varnames,\
                          reg_id_lst=t3_reg_id_lst,\
                          nsub_lon_lst=t3_nsub_lon_lst,\
                          nsub_lat_lst=t3_nsub_lat_lst)




# S2: extract the flux to mult-layer map and regrid them to GEOS-Chem 4x5 grid

print '>2. extract the flux to mult-layer map and regrid them to 4x5'

ctm_grd=gcgrd_2d.gc_grid_2d(0, 0, mod_res='4x5')
ml_ncflnm='reg_144_ml.4x5.nc'

mpflux.extract_region_flux(flux_ncflnm, ml_ncflnm, ctm_grd, \
                               do_debug=False)

# S3: sample multi-layer flux maps to generate ensemble emissions. 

print '>3. sample basis functions to construct emissions ensemble'

# ##c: the configuration have been defined by cfg

cfg.bf_flnm=ml_ncflnm

gems.gen_ensemble_flux_by_reg(cfg)



