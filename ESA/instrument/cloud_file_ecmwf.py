""" Functions to read cloud file  
Functions:
===============================================

1.  open_cloud_file: setup file access class
2.  read_cloud_file: read file to variable list 
3.  close_cloud_file: close file 
4.  get_value: get value at given location 


Parameters:
=========================================================
1. cld_filename_format='cld_ecmwf.'+'XYYYYX'+'XMMX'+'.nc'
2. cf_lon_nm='longitude'  # longitude name used in cloud file  
3. cf_lat_nm='latitude'   # latitude name used in cloud file  
4. cf_time_nm='time'      # time name used in cloud file 
5. cf_var_nm='tcc'        # total column cloud 

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


# configuration of the cloud file 

cld_filename_format='cld_ecmwf.'+'XYYYYX'+'XMMX'+'.nc'

cf_lon_nm='longitude'  # longitude name used in cloud file  
cf_lat_nm='latitude'   # latitude name used in cloud file  
cf_time_nm='time'      # time name used in cloud file 
cf_var_nm='tcc'        # total column cloud 

cld_def_yyyy=2003
cld_def_mm=1
cld_def_dd=1




def open_cloud_file(datapath, flnm, yyyy, mm, dd ,\
                        read_to_array=False, \
                        use_current_time=False):
    
                    
    """ open cloud netcdf file 
    Inputs:
    ---------------------------------------
    1. datapath:<str>: file path
    2. flnm:<str>: file name (not used)
    3. yyyy:<int>: year
    4. mm:<int>: month 
    5. dd:<int>: day 
    6. read_to_array:<T/F>: if Ture, a copy of data, instead of the method,  will be read.
    7.  use_current_time:<T/F>: 
    ---if Ture, the input year (year) will be used
    ---if False, the default year (cld_def_year) will be used
    
    Returns
    -----------------------------------------
    1. fdesc: <ncfile_desc_cl>: file access to netcdf file 
    
    """
    
    # S1 construct IO variables to be read
    
    
    ## lon
    vlon= ovio.io_var_cl(cf_lon_nm,'f',  ['longitude'], None)
    
    vlon.set_attr_rd_list(['units'])
    
    ## lat 
    vlat= ovio.io_var_cl(cf_lat_nm,'f',  ['latititude'], None)
    vlat.set_attr_rd_list(['units'])
    
    ## time 
    vtime=ovio.io_var_cl(cf_time_nm, 'f', ['time'], None)
    vtime.set_attr_rd_list(['units'])
    
    ## total cloud column  
    tcc=ovio.io_var_cl(cf_var_nm, 'f', [cf_time_nm, cf_lat_nm, cf_lon_nm], None)
    
    tcc.set_attr_rd_list(['units', 'scale_factor',  'add_offset', 'missing_value'])
    
    
    var_lst=[vlon, vlat, vtime, tcc]
    # S2 construct  file access 
    ## set read_to_array to false so that no data is read read in 
 
    fdesc=nfio.ncfile_desc_cl(var_lst=var_lst, read_to_array=read_to_array)

    # S3 set file  path 
    
    fdesc.set_file_path(datapath)

    # set file name format 
    fdesc.set_filename_format(cld_filename_format)
    
    # S4 construct full file name 
    
    if (use_current_time):
        # use_current_time
        full_flnm=fdesc.construct_filename(XYYYYX=yyyy, XMMX=mm, XDDX=dd)
        tau_shift=0.0
    else:
        full_flnm=fdesc.construct_filename(XYYYYX=cld_def_yyyy, XMMX=mm, XDDX=dd)
        
        # difference between current year and year saved the data 
        tau_shift=tm.get_tau(yyyy, 1,1)-tm.get_tau(cld_def_yyyy, 1,1)

    

    fdesc.set_attr('filename', full_flnm)
    fdesc.set_attr('tau_shift', tau_shift)
    
    
    return fdesc

                   
    
def read_cloud_file(fdesc):
    
    """
    Read data from netcdf file
    
    Inputs:
    ---------------------------------------
    1. fdesc:<ncfile_desc_cl>: netcdf file access 
    """
    # S1 get file name 
    
    full_flnm=fdesc.get_attr('filename')
    
    # S2 read the required variable 

    
    fdesc.read_var_from_file(full_flnm)
    
    # S3 change time to tau 

    t=fdesc.get_var('time')
    units=t.get_attr('units')
    
    pos=units.index('-')
    st_date=units[pos-4:]
    
    yyyymmdd=st_date[0:10]
    
    ## T1 staring time in time unit 
    
    yyyy0, mm0, dd0=yyyymmdd.split('-')
    yyyy0, mm0, dd0=int(yyyy0), int(mm0),int(dd0)
    
    st_tau=tm.get_tau(yyyy0, mm0, dd0)
    
    tau_shift=fdesc.get_attr('tau_shift')
    
    ## T2  convert time to tau and shift the current year
    
    t.var=npy.array(t.var)
    t.var=t.var+st_tau+tau_shift

    
    
    
    # S4 setup grid and conversion function for tcc  
    
    tcc=fdesc[cf_var_nm]

    
    ## T1 conversion data 
    
    scale_factor=tcc.get_attr('scale_factor')
    add_offset=tcc.get_attr('add_offset')
    missing_value=tcc.get_attr('missing_value')
    
    tcc.set_scale(scale_factor[0])
    tcc.set_offset(add_offset[0])
    tcc.set_mask(missing_value[0])
    
    ## T2  setup the grid 
    
    lon=fdesc.get_lon(cf_lon_nm)
    lat=fdesc.get_lat(cf_lat_nm)
    tau=fdesc.get_time(cf_time_nm)
    
    tau=npy.array(tau)
    lat=npy.array(lat)
        
    if (lat[0]>lat[-1]):
        ### latitude from [90 to -90]

        lat=-lat  
        tcc.set_attr('lat_reverse', 1)
    else:
        tcc.set_attr('lat_reverse', 0)
    
    lon=npy.array(lon)
    
    if (lon[-1]>180.0):  
        ### longitude in  (0-360) 
        tcc.set_attr('lon_shift', 1)
    else:
        tcc.set_attr('lon_shift', 0)
        
    
    ndim=len(tcc.dimname)
    
    ax_tau=axis_m.gp_axis_cl('time', tau)
    ax_lon=axis_m.gp_axis_cl('lon', lon)
    ax_lat=axis_m.gp_axis_cl('lat', lat)

    
    # sort the list of axis according to definition in the file
    # 

    # ax_lst---> sorted axis_list 
    # ax_ord --> location of [lon, lat, time] in the new axis_list

    # 

    ax_lst=3*[0] #
    ax_ord=3*[0]
    ip=0
    
   
    for dname in tcc.dimname: 
        dname=dname.strip()
        
        
        if (dname==cf_lon_nm):
            ax_lst[ip]=ax_lon
            ax_ord[0]=ip
            ip=ip+1
            
        elif (dname==cf_lat_nm):
            
            ax_lst[ip]=ax_lat
            ax_ord[1]=ip
            ip=ip+1
            
        elif (dname==cf_time_nm):
            
            ax_lst[ip]=ax_tau
            ax_ord[2]=ip
            ip=ip+1
    
    
    # set the grid 
    
    
    ctm_grd=grid.gp_grid_cl(ax_lst)
    
    tcc.set_grid(ctm_grd)
    
    # save the order

    tcc.set_attr('ax_order', ax_ord)
    
                   
    return fdesc

def close_cloud_file(fdesc):
    # close the file 
    
    fdesc.close_file()
    del fdesc
    

def get_value(fdesc, olon, olat, otau):
    """
    
    get cloud data at given coordinate 

    Inputs: 
    ----------------------------------------------
    1.fdesc:<<ncfile_desc_cl>: file access
    2. olon:<float>: longitude
    3. olat:<float>: latitude
    4. otau:<float>: time (in hours since 1985.01.01)
    
    
    Returns:
    --------------------------------------------
    1. pcld:<float>: cloud cover percentage 
    
    
    """
    
    tcc=fdesc[cf_var_nm]
    
    # S1 lontitude shift from (-180 to 180) to (0, 360) when necessary 

    
    lon_shift=tcc.get_attr('lon_shift')
    
    if (lon_shift==1):
        if (olon<0.):
            olon=360.0+olon
        
    # S2 latitude
            
    # reverse lat when necessary 
    lat_reverse=tcc.get_attr('lat_reverse')
    
    if (lat_reverse==1):
            olat=-olat
    
    # reorder the coordinate according the axis list 
    # bydefault coordinate is in [lon, lat, time], which will be shifted 
    # what saved in the grid 

    ## index from (lon, lat, time) to the current order of axis 
    
    ax_ord=tcc.get_attr('ax_order')
    

    p=[olon, olat, otau]
    p2=list(p)
    
    # sort the data in right order
    
    p2[ax_ord[0]]=p[0]
    p2[ax_ord[1]]=p[1]
    p2[ax_ord[2]]=p[2]
    
    # S3 get interoplation poistion and weight from grid 
    
    grd_pts_lst, wgt_lst=tcc.grd.allocate_point(p2)
        
    
    # S4 do interpolation by looping over each corner (boundary) points. 
    
    pcld=0.0  # cloud coverage 
    wgt=0.0  # total wgt 
    
    npt=len(wgt_lst)  # number of bound points 
    
    
    for icp in range(npt):
    
        xp=grd_pts_lst[icp]
        
        ### change from tuple to i,j, k as required by non-numpy array
        
        
        i,j,k=xp
        i, j,k =int(i), int(j), int(k)
        
        ### post process data
        
        val=tcc.post_process(tcc.var[i,j,k], 1.0)
            
            
        pcld=pcld+wgt_lst[icp]*val
        wgt=wgt+wgt_lst[icp]
            

    
    pcld=pcld/wgt
    
            
    return pcld

    

if (__name__=='__main__'):

    datapath='/home/lfeng/local_disk_2/otool_data/ecmwf/'
    
    flnm='cld_ecmwf.200301.nc'
    fdesc= open_cloud_file(datapath, flnm, 2003, 1, 1)
    
    fdesc=read_cloud_file(fdesc)
    t=fdesc['time']
    print t.var[0:10]
    
    tcc=fdesc['tcc']

    print tcc.var[3, 4, 5]
    data=tcc.post_process(tcc.var[3,4,5],1)

    tau=tm.get_tau(2003, 1,2,12)
    
    olon, olat, otau=50, 40, tau
    
    pt=get_value(fdesc, olon, olat, otau)
    print pt
    
    
    close_cloud_file(fdesc)

