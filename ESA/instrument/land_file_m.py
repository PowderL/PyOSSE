"""
   Functions for reading LAI percentage file 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.26
   
   Functions
   =================================================
   
   # file access

   1. open_lai_file:  open MODIS LAI file
   2. close_lai_file: close MODIS LAI file 
   3. get_lai_file: get MODIS LAI  map 
   
   """

import ESA.util.gp_axis_m as axis_m
import ESA.util.gp_grid_m as grid

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.otool_var_io as ovio
import ESA.util.geo_constant as gc
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.otool_ncfile_io as ncfio
import numpy as npy

lai_lon_nm='lon'
lai_lat_nm='lat'
lai_time_nm='time'
lai_gp_nm='lai'

lai_varname_lst=[lai_lon_nm,\
                     lai_lat_nm,\
                     lai_time_nm,\
                     lai_gp_nm]
# variable name to name in lai file 

lai_varname_dict={lai_lon_nm:'lon',\
                      lai_lat_nm:'lat',\
                      lai_time_nm:'time',\
                      lai_gp_nm:'MODIS'}



lai_filename_format='MODIS.LAIv.V5.generic.025x025.'+'XYYYYX'+'.nc'



#<<< FUNCTIONS >>> 
def get_lai_file_ref(yyyy, mm, dd, **keywords):
    """ check file reference
    for lai filec, the reference number is year
        
    Inputs:
    ----------------------------------------------------
    1. yyyy:<array, (ntime,)>: year
    2. mm: <array, (ntime,)>: month
    3. dd:<array, (ntime,)>: dd
    4. keywords:<dict>: extra information
    
    Returns:
    ------------------------------------------------------
    1. ref:<int>: reference number 
    --- ref can be [2008, 2009] for year. 
    
    """
    ref=yyyy
    return ref
    

        

def open_lai_file(datapath, flnm, \
                     yyyy, mm, dd ,\
                     varname_lst=lai_varname_lst,\
                     var_ncname_dict=lai_varname_dict,\
                     read_to_array=False,\
                     fio_keywords={},\
                     **keywords):
    
    
    
    """ open bf netcdf file 
    Inputs:
    ---------------------------------------
    1. datapath:<str>: file path
    2. flnm:<str>: file name as format template of lai file 
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
    
    ref=get_lai_file_ref(yyyy, mm, dd, **keywords)
    
    fdesc=ncfio.ncfile_desc_cl(ref=ref, var_lst=var_lst, \
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


                   
    
def read_lai_file(fdesc, **keywords):
    
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
    
    lon=fdesc.get_lon(lai_lon_nm)
    lat=fdesc.get_lat(lai_lat_nm)
    
    lat=npy.array(lat[:])
    lon=npy.array(lon[:])
    
    
    
    ax_lon=axis_m.gp_axis_cl(lai_lon_nm, lon)
    ax_lat=axis_m.gp_axis_cl(lai_lat_nm, lat)
    
    ax_lst=[ax_lon, ax_lat]
    
    
    # sort the list of axis according to definition in the file
    # 

    # ax_lst---> sorted axis_list 
    # ax_ord --> location of [lon, lat, time] in the new axis_list

    # 

    ctm_grd=grid.gp_grid_cl(ax_lst)
    
    
    
    fdesc.set_grid(ctm_grd)
    
                       
    return fdesc

def close_lai_file(fdesc):
    """
    1. fdesc:<ncfile_desc_cl>: netcdf file access 
    
    """
    fdesc.close_file()
    del fdesc
    return 



def add_lai_data(fdesc, ref, **keywords):
    """ add the cloud data defined by ref to the fdesc
    
    Inputs:
    =====================================
    1. fdesc:<grdfile_io_cl>: class for file access
    2. ref:<str/numeric>: reference for data sets read from file 
    3. keywords:<dict>: extra parameters

    Returns:
    =====================================
    1. idx:<int>: index of the data set with reference =ref. 
    

    """
    
    # S1 construct the name 
    
    yyyy=ref
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy)
    fdesc.set_filename(full_flnm)
    
    
    # S2 read the data
    
    fdesc=read_lai_file(fdesc, **keywords)
    
    
    return fdesc
 



def get_lai_data(fdesc, olon, olat, \
                     oyyyy, omm, odd,\
                     **keywords):
    
    """
    read cloud data at given point 
    Inputs:
    ---------------------------------------
    1. fdesc:<grdfile_io_cl>: file access 
    2. olon: <array, (nob,)>: longitude
    3. olat: <array, (nob,)>: latitude
    4. oyyyy:<int>: year 
    5. omm: <int>: month 
    6. odd:<int>: day (for future use)
    
    8. keywords:<dict>: extra inputs

    ---Reserved Entries
    ---common_ref:<str/numeric>: the reference shared by observations. 

    Returns:
    ================================================
    lai:<array, (nob, )>: lai at the (olon, olat, otime)
    
    """
    
    ref=get_lai_file_ref(oyyyy, odd, omm)

    if (ref<>fdesc.ref):
        fdesc=add_lai_data(fdesc,ref, **keywords)
    
        
    nob=npy.size(olon)
    # construct ref and ref_set (a set for different references for obs)
    # 
    
    grd=fdesc.get_grid()
    
    gdata=fdesc.get_data('lai')
    # print npy.shape(gdata)
    
    gdata=gdata[omm-1, :,:]
    gdata=npy.transpose(gdata)
    
    
    ax_lon=grd['lon']
    ax_lat=grd['lat']
    
    plon=ax_lon.get_closest_point(olon, fdesc.mask_val)
    plat=ax_lat.get_closest_point(olat, fdesc.mask_val)
    lai=gdata[plon, plat]

    
    return lai

    


    
if (__name__=='__main__'):
    import ESA.util.gen_plots as gpl
    import pylab as plb
    
    datapath='/home/lfeng/local_disk/otool_data/LAI/'
    flnm='MODIS.LAIv.V5.generic.025x025.2008.nc'
    fdesc= open_lai_file(datapath, flnm, \
                             2009, 1, 1)
    
    
    fdesc=read_lai_file(fdesc)

    olon=npy.arange(-60, 60, 10)
    olat=npy.arange(-60, 60, 10)
    yyyy=2009
    mm=1
    dd=2
    
    lai=get_lai_data(fdesc, olon, olat, \
                         yyyy, mm, dd)
    
    print npy.shape(lai)
    print lai
    
    close_lai_file(fdesc)
    
    
    


#<<< FUNCTIONS >>> 

