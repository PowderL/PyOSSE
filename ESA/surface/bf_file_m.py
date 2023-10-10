""" Functions to read flux perturbations (Basis functions)

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2012.09.12
    this module define axis classes
    
    

    Functions:
    ===============================================
    
    1.  open_bf_file: setup file access class
    2.  read_bf_file: read file to variable list 
    3.  close_bf_file: close file 
    4.  get_bf: get value at given location 
    
    """
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_var_io as ovio
import ESA.util.otool_ncfile_io as nfio

import ESA.util.time_module as tm
import ESA.util.process_nf_array as pnf

import ESA.util.gp_grid_m as grid
import ESA.util.gp_axis_m as axis_m


# parameters for bf file 

bf_filename_format='bf.'+'XYYYYX'+'XMMX'+'.nc'

bf_lon_nm='lon'  # longitude name used in bf file  
bf_lat_nm='lat'   # latitude name used in bf file  
bf_layer_nm='layer'    # time name used in bf file 
bf_map_nm='map'   # total column bf 
bf_flux_nm='flux'
bf_area_nm='area'
bf_pid_nm='pid'

# variables to be read in

bf_varname_lst=[bf_lon_nm,\
                    bf_lat_nm,\
                    bf_layer_nm,\
                    bf_map_nm,\
                    bf_flux_nm, \
                    bf_area_nm,\
                    bf_pid_nm]

#<<< FUNCTIONS >>> 

def open_bf_file(datapath, flnm, \
                     yyyy, mm, dd ,\
                     varname_lst=bf_varname_lst,\
                     var_ncname_dict={},\
                     read_to_array=False,\
                     fio_keywords={},\
                     **keywords):
    
    
    
    """ open bf netcdf file 
    Inputs:
    ---------------------------------------
    1. datapath:<str>: file path
    2. flnm:<str>: file name as format template of basis function file 
    3. yyyy:<int>: year
    4. mm:<int>: month 
    5. dd:<int>: day 
    6. varname_lst: list of the (default) varname to be read in from netcdf file 
    7. var_ncname_dict:<dict>: translations from var names to names in netcdf file 
    8. read_to_array:<T/F>: if Ture, a copy of data, instead of the method,  will be read.
    9. keywords:<dict>: attributes 
    
    Returns
    -----------------------------------------
    1. fdesc: <ncfile_desc_cl>: file access to netcdf file 
    
    """
    
    # S1: construct IO variables to be read
    
    
    var_lst=[]  # list of the var to be read 
    
    for varname in varname_lst:
        if (varname in var_ncname_dict):
            nc_vname=var_ncname_dict[varname]
        else:
            nc_vname=varname
        #
        new_var= ovio.io_var_cl(nc_vname,'f',  [nc_vname], None)
        var_lst.append(new_var)
        

    # S2: construct fdesc
    
    fdesc=nfio.ncfile_desc_cl(var_lst=var_lst, \
                                  varname_dict=var_ncname_dict,\
                                  read_to_array=read_to_array,fio_keywords=keywords)
    
    
    
    # S4: set file path 
    
    fdesc.set_file_path(datapath)

    # S5:  set file name format 
    
    fdesc.set_filename_format(flnm)
    
    # S6: construct full file name 
    doy=tm.day_of_year(yyyy, mm, dd)
  
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy, \
                                           XMMX=mm, \
                                           XDDX=dd, XDOYX=doy)
    
    # #t: save the filename 
    

    
    fdesc.set_attr('filename', full_flnm)
    
    
    return fdesc

                   
    
def read_bf_file(fdesc, **keywords):
    
    """
    Read data from netcdf file
    
    Inputs:
    ---------------------------------------
    1. fdesc:<ncfile_desc_cl>: netcdf file access 
    """
    
    # S1: get file name 
    
    full_flnm=fdesc.get_attr('filename')
    
    # S2: read the required variable 
    
    fdesc.read_var_from_file(full_flnm)
    
    # S3: change time to tau 
    
    lon=fdesc.get_lon(bf_lon_nm)
    lat=fdesc.get_lat(bf_lat_nm)
    
    lat=npy.array(lat[:])
    lon=npy.array(lon[:])
    
    
    
    ax_lon=axis_m.gp_axis_cl(bf_lon_nm, lon)
    ax_lat=axis_m.gp_axis_cl(bf_lat_nm, lat)
    
    ax_lst=[ax_lon, ax_lat]
    
    
    # sort the list of axis according to definition in the file
    # 

    # ax_lst---> sorted axis_list 
    # ax_ord --> location of [lon, lat, time] in the new axis_list

    # 

    ctm_grd=grid.gp_grid_cl(ax_lst)
    
    fdesc.set_grid(ctm_grd)
    
                       
    return fdesc

def close_bf_file(fdesc):
    """
    1. fdesc:<ncfile_desc_cl>: netcdf file access 
    
    """
    fdesc.close_file()
    del fdesc
    return 


def get_bf(fdesc,ilayer, **keywords):
                
    """
    
    get bf data for given or default variable list 
    
    Inputs: 
    ----------------------------------------------
    1.fdesc:<ncfile_desc_cl>: file access
    2. layer:<int>: ilayer 
    3. keywords:<dict>: extra inputs
    ---reserved keys:
    ---varname_lst:<list, t:str>: name of variables.
    ---flux_nm:<str>:name of flux
    ---map_nm:<str>: name of map
    ---area_nm:<str>: name of area
    ---lon_nm:<str>: name of lon
    ---lat_nm:<str>: name of lat
        
    
    Returns:
    --------------------------------------------
    
    1. var_lst:<list>: list of variable read (see Note 1)
    

    Notes:
    ---------------------------------------
    1. if varname_lst is not set in keywords. A default list of variable will be returned 
    This list includes:
    ---1. lon:<array, (nlon)>: longitude
    ---2. lat:<array, (nlat)>: latitude
    ---3. sel_flux:<array, (nlon, nlat)>: basis functions 
    ---4. sel_map:<array, (nlon, nlat)>: map define the region  
    ---5. sel_area:<array, (nlon, nlat)>: areas for the region 
    
    
    """
    
    # S1: check whether the name list is given 
    if ('varname_lst' in keywords):
        varname_lst=keywords['varname_lst']
        var_lst=list()
        for varname in varname_lst:
            var=fdesc.get_data(varname)
            var_lst.append(npy.array(var))
        
        return var_lst
    
        
    # S2: if varname_lst is not given, a default list will be return
    
    # #c: flux 

    flux_nm=bf_flux_nm
    if ('flux_nm' in keywords):
        flux_nm=keywords['flux_nm']
    
    reg_flux=fdesc.get_data(flux_nm)
    
    # #c: map
    map_nm=bf_map_nm
    if ('map_nm' in keywords):
        map_nm=keywords['map_nm']
    
    reg_map=fdesc.get_data(map_nm)
    
    # #c: area
    
    area_nm=bf_area_nm
    
    if ('area_nm' in keywords):
        area_nm=keywords['area_nm']
    
    reg_area=fdesc.get_data(area_nm)
    
    # lon 
    lon_nm=bf_lon_nm
    if ('lon_nm' in keywords):
        lon_nm=keywords['lon_nm']
    

    lon=fdesc.get_data(lon_nm)

    # lat
    lat_nm=bf_lat_nm
    if ('lat_nm' in keywords):
        lat_nm=keywords['lat_nm']
    

    lat=fdesc.get_data(lat_nm)
    
    reg_flux=npy.array(reg_flux)
    reg_map=npy.array(reg_map)
    reg_area=npy.array(reg_area)
    
    lon=npy.array(lon)
    lat=npy.array(lat)

    dims=npy.shape(reg_flux)
    
    if (len(dims)==3):
        # it is layer map
        sel_map=reg_map[:,:,layer]
        sel_flux=reg_flux[:,:,layer]
    else:
        # map starting from 1 
        sel_map=npy.where(reg_map==layer+1, 1, 0)
        sel_flux=reg_flux*sel_map
        
    sel_area=reg_area
    return lon, lat, sel_flux, sel_map, sel_area
    
    

 


if (__name__=='__main__'):
    import ESA.util.gen_plots as gpl
    import pylab as plb

    datapath='./'
    flnm='reg_144_05x05.nc'
    fdesc= open_bf_file(datapath, flnm, \
                            2009, 1, 1, var_ncname_dict={bf_lon_nm:'longitude', \
                                                             bf_lat_nm:'latitude',\
                                                             bf_map_nm:'map',\
                                                             bf_flux_nm:'flux'})
    
    
    
    fdesc=read_bf_file(fdesc)
    
    reg_flux, reg_map, reg_area=get_bf(fdesc, 3)
    
    close_bf_file(fdesc)
    
    ctm_grd=fdesc.grd
    
    ax_lon=ctm_grd[bf_lon_nm]
    ax_lat=ctm_grd[bf_lat_nm]
    lon=ax_lon[:]
    lat=ax_lat[:]
    
    # plb.subplot(2,1,1)
    gpl.plot_map(reg_flux, lon, lat, use_pcolor=1)
    plb.show()
    
    
