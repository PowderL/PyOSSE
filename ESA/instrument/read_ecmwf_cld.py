"""
   Functions for reading  ECMWF cloud percentage file 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.26
   
   Functions
   =================================================
   
   # file access

   1. open_ecmwf_cld:  open ECMWF cloud file
   2. close_ecmwf_cld: close ECMWD cloud file 
   3. get_ecmwf_cld_map: get ecmwf cloud map 
   
   """

import numpy as npy
import Scientific.IO.NetCDF as snf
import ESA.util.gen_plots as gpl
import pylab as plt

import datetime

def open_ecmwf_cld(datapath, flnm, yyyy, mm, dd, hh):
    """ open cloud file ecmwf file 
    
    Inputs:
    =============================================
    1. datapath:<str>: path
    2. flnm: <str>: file name
    3. yyyy, mm, dd, hh:<int>: time 
    
    Returns:
    ===========================================
    1. cld_fl:<obj>: cloud file handle
    
    """
    
    yyyy, mm, dd=tm.doy_to_time_array(doy, yyyy)
    syyyymm=r'%4.4d%2.2d' % (yyyy, mm)
    cld_flnm=datapath+'/'+flnm+"."+syyyymm+".nc"
    cld_fl=snf.NetCDFFile(cld_flnm, "r")
    
    return cld_fl

def close_ecmwf_cld(fl):
    """
    close file 
    
    Inputs:
    =============================================
    1. fl:<obj>: file handle
    
    
    """
    snf.close(fl)
    
    

def get_ecmwf_cld_map(fl, yyyy, mm, dd, hh, tcc=None):
    """

    get ecmwf cloud map 
    
    Inputs:
    --------------------------------------------
    1. fl :<handle>: file handle 
    2. yyyy, mm, dd:<integer>: year, month and day 
    3. hh: <integer>: hour
    4. tcc:<variable, (nt, nlat, nlon)>: time series of the total cloud ratio, 
    

    Returns:
    ----------------------------------------------
    1. lon:<array, nlon>: longitude grid 
    2. lat:<array, nlon>: latitude grid 
    3. sel_tcc:<array, (nlon, nlat)>: cloud map at current time  
    4. tcc:<variable, (nt, nlat, nlon)>:  time series of cloud map 
        
    """
    
    # longitude 
    lon=fl.variables['longitude']
    lon=npy.array(lon)

    # latitude 
    lat=fl.variables['latitude']
    lat=npy.array(lat)

    # time  
    t=fl.variables['time']
    units=getattr(t, 'units')
    t=npy.array(t)
    t=t.astype(float)
    
    ##  convert current time to hours since starting time

  
    pos=units.index('-')
    st_date=units[pos-4:]
    
    yyyymmdd=st_date[0:10]
    ## staring time 
    yyyy0, mm0, dd0=yyyymmdd.split('-')
    yyyy0, mm0, dd0=int(yyyy0), int(mm0),int(dd0)

    st_date=datetime.date(yyyy0, mm0, dd0)
    st_hour=plt.date2num(st_date)
    
    ##  get current time relative to starting time 
    

    sel_date=datetime.date(yyyy, mm, dd)
    sel_hour=plt.date2num(sel_date)
    
    
    sel_hour=24*(sel_hour-st_hour)+hh
    
    
    
        
    # find time index 
    
    dt=abs(t-sel_hour)
    t_idx=npy.argmin(dt)
    t_idx=npy.squeeze(t_idx)
    t_idx=int(t_idx)
    
    # load tcc if necessary 
    
    if (tcc==None):
        tcc=fl.variables['tcc']
        
    # scaling factor 
    scale_factor=getattr(tcc, 'scale_factor')

    # offset 

    add_offset=getattr(tcc, 'add_offset')
    
   
    # tcc=npy.array(tcc)
    
    # tcc at the current time 
    sel_tcc=tcc[t_idx,:,:]
    
    sel_tcc=npy.array(sel_tcc)
    sel_tcc=sel_tcc.astype(float)
    
    # re-calculate tcc

    usd_idx=npy.where(sel_tcc<>-32767)
    sel_tcc[usd_idx]=scale_factor*sel_tcc[usd_idx]+add_offset
    sel_tcc=npy.where(sel_tcc==-32767, -999.0, sel_tcc)
    sel_tcc=npy.transpose(sel_tcc)

    # re-order lat

    lat=lat[::-1]
    sel_tcc=sel_tcc[:, ::-1]
    
    shft_idx=npy.size(lon)/2
    # lon shift
    lon=npy.where(lon>180, lon-360.0, lon)
    new_lon=npy.array(lon)
    new_lon[0:shft_idx]=lon[shft_idx:]
    new_lon[shft_idx:]=lon[0:shft_idx]
    lon=new_lon
    
    # tcc shift 
    
    new_sel_tcc=npy.array(sel_tcc)
    new_sel_tcc[0:shft_idx,:]=sel_tcc[shft_idx:, :]
    new_sel_tcc[shft_idx:,:]=sel_tcc[0:shft_idx, :]
    sel_tcc=new_sel_tcc
    
    
    return lon, lat, sel_tcc, tcc

if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/ecmwf/'
    
    flnm=datapath+'cld_ecmwf.200301.nc'
    fl=snf.NetCDFFile(flnm)
    tcc=None
    
    lon, lat, sel_tcc, tcc=get_ecmwf_cld_map(fl, 2003, 1,1,13, tcc)
    
    print npy.shape(sel_tcc)
    print npy.shape(tcc)
    
    lon, lat, sel_tcc2, tcc=get_ecmwf_cld_map(fl, 2003, 1,4,13, tcc)
    fl.close()
    
    gpl.plot_map(sel_tcc2, rlon=lon, rlat=lat, use_pcolor=1, vmax=1.0, vmin=0.0, dv=0.1)
    plt.title('ECMWF opertional Cloud coverage 2003.01.04 12:00:00')
    plt.savefig('ecmwf_cld.png')
    
    plt.show()
    
    
             
