""" Functions for write file to GEOS-Chem BPCH2 file 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.12
    History: v0.95, 2013.02.01
    

    Functions:
    ==========================================================
    1. convert_nc_ens_flux_bpch: read into ensemble flux and save them to bpch2 file. 
    
    """
import ESA.surface.write_flux_bpch as wfbp 
import ESA.util.otool_var_io as ovio
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.time_module as tm
import numpy as npy



def convert_nc_ens_flux_bpch(ncflnm, tau0, tau1,\
                                 varname_lst=['lon', 'lat', 'map', 'flux'], \
                                 varname_dict={'lon':'lon', 'lat':'lat', 'regmap':'map','reg_flux':'flux'},\
                                 bpch_flux_flnm='region_flux_ml', \
                                 bpch_map_flnm='region_map_ml', \
                                 save_split_necdf=True,\
                                 ):
    
    
    """
    read into ensemble flux and save them to bpch2 file. 
    
    Inputs:
    ----------------------------------------------------------
    1. ncflnm:<str>: netcdf file name 
    2. tau0, tau1:<float>: start and end time 
    3. varname_lst:<list, t:str>: list of variables to be read in 
    4. varname_dict:<dict>: translate name used in codes to name from netcdf file 
    5. bpch_flux_flnm:<str>: name for output flux 
    6. bpch_map_ml:<str>: name for output map 
    7. save_split_necdf:<T/F>: If ture, netcdf conuterparts will be  
    
    Notes:
    
    
    """

    # S1: read netcdf ensemble flux file 

    nc_varname_lst=list()
    nvar=len(varname_lst)
    
    for ivar in range(nvar):
        varname=varname_lst[ivar]
        if varname in varname_dict:
            varname=varname_dict[varname]
        
        nc_varname_lst.append(varname)
    
        

    var_lst=ncfio.ncf_read(ncflnm, nc_varname_lst)
    nvar=len(var_lst)

    # S2: variable 

    for ivar in range(nvar):
        varname=varname_lst[ivar]
        
        if (varname=='lon'):
            lon=var_lst[ivar]
        elif (varname=='lat'):
            lat=var_lst[ivar]
        elif (varname=='map'):
            new_map=var_lst[ivar]

        elif (varname=='flux'):
            new_flux=var_lst[ivar]
        
            
    # S3: save map 
    ndims=npy.shape(new_map)
    if (len(ndims)==3):
        mlayer=ndims[2]
    else:
        mlayer=1
        
        
    wfbp.write_ens_map_ml(bpch_map_flnm, \
                              new_map,\
                              lon, lat, mlayer,\
                              category='LANDMAP')
    
    
    # S4: save map 
    ndims=npy.shape(new_flux)
    if (len(ndims)==3):
        mlayer=ndims[2]
    else:
        mlayer=1
    

    wfbp.write_ens_flux_ml(bpch_flux_flnm, \
                               new_flux,\
                               lon, lat, mlayer,\
                               tau0, tau1,\
                               1,\
                               category='CO2_FLUX',\
                               unit="mole/s")
    

    # S5: save netcdf counterpart 
    
    if (save_split_necdf):
        
        lat_info=ovio.io_var_cl('lat', 'f', ['lat'], lat, \
                                    varattr=None)
        lon_info=ovio.io_var_cl('lon', 'f', ['lon'], lon, \
                                    varattr=None)
        
        xlayer=npy.arange(mlayer)
        layer_info=ovio.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                      varattr=None)
        
        flux_info=ovio.io_var_cl('flux', 'f', ['lon', 'lat', 'layer'], new_flux, \
                                     varattr={"units":"moleC/s/cm2", "long_name":"CO2 surface flux", \
                                                  "standard_name":"co2_surface_flux"})
        
        map_info=ovio.io_var_cl('map', 'f', ['lon', 'lat', 'layer'], new_map, \
                                      varattr={"units":"ratio", "long_name":"region map ratio", \
                                                   "standard_name":"region_map_ratio"})

        out_ncflnm='sp_'+bpch_flux_flnm+'.nc'
        ncfio.ncf_save_var(out_ncflnm,  [flux_info], \
                               [lon_info, lat_info, layer_info], \
                               create_new=True)   

        
        out_ncflnm=bpch_map_flnm+'.nc'
        
        ncfio.ncf_save_var(out_ncflnm,  [map_info], \
                               [lon_info, lat_info, layer_info], \
                               create_new=True)   
    

        
if (__name__=='__main__'):
    """ generate the files """
    
    import ensemble_flux_io as ens_fio
    yyyy_st, mm_st, dd_st=2010, 1, 1
    time_step=32
    nstep=12
    # S1: check the time step 
    ens_doy_lst, ens_yyyy_lst=ens_fio.define_ens_step(yyyy_st, mm_st, dd_st, time_step, nstep)
    
    # S2: write co2_emissions for each time step to bpch2 

    for istep in range(0, nstep+1):

        # #T: construct file name for each step
        
        cur_doy=ens_doy_lst[istep]
        cur_yyyy=ens_yyyy_lst[istep]
        next_doy=ens_doy_lst[istep+1]
        next_yyyy=ens_yyyy_lst[istep+1]
        sdoy='%4.4dD%3.3d' % (cur_yyyy, cur_doy)
        # #c: original flux file name 
        ncflnm='CO2_EMISSION'+"."+sdoy+'.nc'
        # #c: output file name 
        bp_flux_flnm='CO2_EMISSION'+"."+sdoy
        bp_map_flnm='CO2_MAP'+"."+sdoy
        
        yyyy, imm, idd=tm.doy_to_time_array(cur_doy, cur_yyyy)
        tau0     = tm.get_tau(yyyy, imm, 1)
        tau0=tau0/3600.0
        yyyy, imm, idd=tm.doy_to_time_array(next_doy, next_yyyy)
        tau1     = tm.get_tau(yyyy, imm, 1)
        tau1=tau1/3600.0
        
        # #T2: convert data to bpch2
        print 'origin file:',  ncflnm 
        print 'bpch flux file:',  bp_flux_flnm
        print 'bpch map file:',   bp_map_flnm
        
        
        convert_nc_ens_flux_bpch(ncflnm, tau0, tau1, \
                                     bpch_flux_flnm=bp_flux_flnm, \
                                     bpch_map_flnm=bp_map_flnm)
        


    # loop istep end
    

    
        

        
    
                                                     

    

        
