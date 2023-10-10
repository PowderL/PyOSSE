""" functions to reach cloud file  


"""
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_var_io as ovio
import ESA.util.otool_ncfile_io as nfio

import ESA.util.time_module as tm

import ESA.util.process_nf_array as pnf


cld_filename_format='cld_ecmwf.'+'XYYYYX'+'XMMX'+'.nc'
def open_cloud_file(datapath, flnm, yyyy, mm, dd, hh):
    """ open cloud file ecmwf file 
    
    """
    
    vlon= ovio.io_var_cl('longitude','f',  ['longitude'], None)
    vlon.set_attr_rd_list(['units'])
    
    vlat= ovio.io_var_cl('latitude','f',  ['latititude'], None)
    vlat.set_attr_rd_list(['units'])
    
    vtime=ovio.io_var_cl('time', 'f', ['time'], None)
    vtime.set_attr_rd_list(['units'])
    
    tcc=ovio.io_var_cl('tcc', 'f', ['longitude', 'latititude', 'time'], None)
    tcc.set_attr_rd_list(['units', 'scale_factor',  'add_offset', 'missing_value'])
    
    
    var_lst=[vlon, vlat, vtime, tcc]
    
    fdesc=nfio.ncfile_desc_cl(var_lst=var_lst)
    # set path 
    fdesc.set_file_path(datapath)
    # set file name format 
    fdesc.set_filename_format(cld_filename_format)
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy, XMMX=mm, XDDX=dd)
    fdesc.set_attr('filename', full_flnm)
    return fdesc

                   
    
def read_cloud_file(fdesc):
    
    full_flnm=fdesc.get_attr('filename')
    fdesc.read_var_from_file(full_flnm)
    
    # change time to tau 

    t=fdesc.get_var('time')
    units=t.get_attr('units')
    
    pos=units.index('-')
    st_date=units[pos-4:]
    
    yyyymmdd=st_date[0:10]
    ## staring time 
    yyyy0, mm0, dd0=yyyymmdd.split('-')
    yyyy0, mm0, dd0=int(yyyy0), int(mm0),int(dd0)
    
    st_tau=tm.get_tau(yyyy0, mm0, dd0)
    
    print t.var[0:10]
    
    t.var=t.var+st_tau
    
    print t.var[0:10]
    
    # change tcc to value 
    
    tcc=fdesc['tcc']

    scale_factor=tcc.get_attr('scale_factor')
    add_offset=tcc.get_attr('add_offset')
    print scale_factor
    print add_offset
    
    missing_value=tcc.get_attr('missing_value')
    
    missing_value= missing_value[0]
    scale_factor=scale_factor[0]
    add_offset=add_offset[0]
    ist=0
    iend=5
    ntime=npy.size(t.var)
    
    pnf.fill_nf_3d_int(tcc.var,missing_value,1)

    while ist<ntime:
        print ist, iend
        
        data=tcc.var[ist:iend, :, :]
    
        print data[3,3,3]

    
        data=pnf.process_nf_3d(data,scale_factor,add_offset,missing_value,1.0)
        tcc.var[ist:iend, :, :]=data
        print data[3,3,3]
        print tcc.var[ist+3, 3,3]
        
        ist=iend
        iend=ist+5
        if (iend>ntime):
            iend=ntime
    


    # usd_idx=npy.where(tcc.var<>missing_value)
    # tcc.var[usd_idx]=scale_factor*tcc.var[usd_idx]+add_offset
    # tcc.var=where(tcc.var==missing_value, fdesc.mask_val, tcc.var)
    
    
    return fdesc

def close_cloud_file(fdesc):
    
    del fdesc
    


if (__name__=='__main__'):

    datapath='/home/lfeng/local_disk_2/otool_data/ecmwf/'
    
    flnm='cld_ecmwf.200301.nc'
    fdesc= open_cloud_file(datapath, flnm, 2003, 1, 1, 13)

    fdesc=read_cloud_file(fdesc)
    t=fdesc['time']
    print t.var[0:10]
    
    
