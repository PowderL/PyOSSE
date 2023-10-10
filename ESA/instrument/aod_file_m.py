"""
  Functions for reading  AOD PDF file 

  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.06.28
  History: v0.95, 2013.02.24
  
  Functions
  =================================================
  1. get_aod_column_name:  check names of columns in AOD file.
  2. get_aod_surface_type: check names of surface types.
  3. aod_rescale: rescale AOD axis values.
  4. open_aod_file: set up class file_desc_cl for file access.
  5. read_aod_file: read AOD  data into class file_desc_cl.
  6. close_aod_file: delete class file_desc_cl.
  7. get_aod_file_ref: check file reference
  8. get_aod_data: get aod data at given locations
  9. add_aod_data: add AOD data to the data list 
  
    
"""

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_grdfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import pylab as plb

#===< PARAMETERS >======


## aod values 

def_aod_lvl=npy.arange(0.0, 21.0)
def_aod_lvl=0.05*def_aod_lvl


##2 columns in AOD file

def_aod_colno_lst=[0,1]+range(2, 23)  # column 23 and 24 will not be read


def_aod_colnm_lst=['lat', 'lon']  # first 2 columns. 


##3 replacements for special character used in text lines

def_aod_replaces={} 


##4  axises (i.e., dependences of aod)
 
aod_axis_name=['lon', 'lat']
## dictionary between varnames used in file and varnames used in the classes

aod_varname_dict={'lat':'lat', 'lon':'lon', 'aod':'aod'}

##5 'verical' levels

aod_zname='aod'
aod_zval=def_aod_lvl



##6 starting row number of data section  
aod_data_start=0

##7 starting row number of header section
aod_header_start=-1

##8 template for file name 

aod_filename_format='modis_misr_ccm_histo_XMONTHX.dat'


#===< FUNCTIONS >====



def get_aod_file_ref(yyyy, mm, dd, **keywords):
    """ check file reference
    for aod climatology, the reference number is 
    climatology months [1, 4, 7, 10]
        
    Inputs:
    ----------------------------------------------------
    1. yyyy:<array, (ntime,)>: year
    2. mm: <array, (ntime,)>: month
    3. dd:<array, (ntime,)>: dd
    4. keywords:<dict>: extra information
    
    Returns:
    ------------------------------------------------------
    1. ref:<int>: reference number 
    --- ref can be [1, 4, 7 or  10] for climatology month. 
    
    """

    if (oob.get_ot_type(mm)==oob.ot_array):
        ## if mm is an array
        
        ref=npy.zeros(npy.size(mm), int)
        
        idx=npy.where((mm<3) | (mm>11))  # winder
        ref[idx]=1
        
        idx=npy.where((mm>=3) & (mm<6))  # spring  
        ref[idx]=4

        idx=npy.where((mm>=6) & (mm<9))  # summer
        ref[idx]=7
        
        idx=npy.where((mm>=9) & (mm<12))  # fall 
        
        ref[idx]=10
        
        
    else:
        # single value 

        if (mm<3):  # winter 
            ref=1   
        elif (mm<6):  # spring
            ref=4
        elif (mm<9):  # summer
            ref=7
        elif (mm<12):  # fall
            ref=10
        else:
            ref=1
    
    return ref
    


def get_aod_column_name():
    
    """
    check column names defined in AK file 
        
    Returns
    ---------------------------------------
    1. fl_colnm_lst:<list, str>: name of columns in file 
    
    """

    fl_colnm_lst=def_aod_colnm_lst

    for iz in range(npy.size(aod_zval)):
        # append 'Z??' to column name list 
        
        szl=r'Z%d' % iz
        fl_colnm_lst=fl_colnm_lst+[szl]
    
    return fl_colnm_lst




def open_aod_file(datapath, flnm, \
                      yyyy, mm, dd, \
                      **keywords):
    
    """ setup fdesc for aerosol optical depth file
    Inputs:
    -----------------------------------------------------------
    1. flnm:<str>:   name of file (reserved for future use )
    2. datapath:<str>: file path
    3. viewtype:<str>: instrument type 
    4. viewmode:<str>: nadir or glint view
    5. yyyy: <int>: year
    6. doy: <int>: day of year 
    7. keywords:<dict>: extra for reading file 
    ---Reserved words
    ---fl_colnm_lst:<list, t:str>: names of all columns in the averaging kernel file 
    
    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access

    """

        
    # S1 set up a full column list 
    
    if ('fl_colnm_lst' in keywords):
        fl_colnm_lst=keywords['fl_colnm_lst']
    else:
        fl_colnm_lst=def_aod_colnm_lst  # +['C1']
        
        for iz in range(npy.size(aod_zval)):
            # append 'Z??' to column name list 
 
            szl=r'Z%d' % iz
            fl_colnm_lst=fl_colnm_lst+[szl]

        
    
    delim=' '
    
    if ('delim' in keywords):
        delim=keywords['delim']
        
    
    # S2 construct fdesc
    
    ref=get_aod_file_ref(yyyy, mm, dd, **keywords)
    
    if (ref==1):
        month='jan'
    
    elif (ref==4):
        month='apr'
    
    elif (ref==7):
        month='jul'
    
    elif (ref==10):
        month='oct'
    
    

    
    fdesc=ofo.grdfile_desc_cl(flnm, ref, aod_header_start, \
                                  aod_data_start, \
                                  delim=delim, \
                                  colnm_lst=fl_colnm_lst)
    
    # S3 set format 
    fdesc.set_filename_format(aod_filename_format)
    
    # S4 set path 
    fdesc.set_file_path(datapath)
    
    # S5 set full name 
    
    full_flnm=fdesc.construct_filename(XMONTHX=month)
    fdesc.set_filename(full_flnm, ref)
    
    return fdesc
 
    
def read_aod_file(fdesc, **keywords):
 
    
    """
    Read averaging kernel data into fdesc
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    2. keywords:<dict>: extra for file reading 
    ---Reserved keys

    --->ax_colnm_lst:<list, t:str>: column names for 'horizontal' (such as sza and aod) axises  
    --->zname: <str>: axis name for vertical (i.e, log10(pressure)) axis
    --->zval: <array, (nz,)>: vertical axis 
    --->'replaces':<dict>: words to be replaced when decoding line 
    
    Returns:
    -------------------------------------
    1. ndata:<int>: size of AK table 


    Notes:
    -------------------------------------
    1.fdesc.grd and fdesc.gdata will store 1) grid; and 2) gridded table 
    
    
    """
    # S1 load full filename

    flnm=fdesc.cur_flnm
    
    fl=open(flnm, 'r')
    lines=fl.readlines()
    fl.close()
   

    
    # S2 load format information 
    ## T1 columns for horizontal axis 
    
    if ('ax_colnm_lst' in keywords):
        ax_colnm_lst=keywords['ax_colnm_lst']
    
    else:
        ax_colnm_lst=aod_axis_name

    
    ## T2 vertical axis 
    
    if ('zname' in keywords):
        zname=keywords['zname']
    else:
        zname=aod_zname
    
    if ('zval' in keywords):
        zval=keywords['zval']
    else:
        zval=aod_zval

    
    ## T3 replacements for line decoding 
    
    if ('replaces' in keywords):
        replaces=keywords['replaces']
    else:
        replaces=def_aod_replaces
   
    
    ## T4  columns needs to be read 
    
    if ('fl_colno_lst' in keywords):
        fl_colno_lst=keywords['fl_colno_lst']
    else:
        fl_colno_lst=def_aod_colno_lst
    

    if ('replaces' in keywords):
        replaces=keywords['replaces']
    else:
        replaces=def_aod_replaces
        
    # S3 read data into maxtrix (fdesc.data)
   
    data=fdesc.read_matrix_from_lines(lines, \
                                          colno_lst=fl_colno_lst, **replaces)
    
    
    # S4 regrid data to fdesc.gdata

    # for testing 
    
    ax_colno_lst=fdesc.get_colno_list(ax_colnm_lst)
    
    
    grd, gdata=fdesc.convert_matrix_to_grid(data, ax_colnm_lst, \
                                                ax_colno_lst, \
                                                zname=zname, \
                                                zval=zval)
    # S4 rescale axis aod of fdesc.grd

    fdesc.append_data(fdesc.cur_flnm, data, grd, gdata, fdesc.cur_ref)
    
    
    
    return fdesc


def add_aod_data(fdesc, ref, **keywords):
    """ add AOD data defined by ref to the fdesc
    
    Inputs:
    =====================================
    1. fdesc:<grdfile_io_cl>: class for file access
    2. ref:<str/numeric>: reference for data sets read from file 
    3. keywords:<dict>: extra parameters
    
    --->ax_colnm_lst:<list, t:str>: column names for 'horizontal' (such as lon and lat) axises     
    --->zname: <str>: axis name for vertical (i.e, log10(pressure)) axis
    --->zval: <array, (nz,)>: vertical axis 
    --->'replaces':<dict>: words to be replaced when decoding line 
    
    Returns:
    =====================================
    1. idx:<int>: index of the data set with reference =ref. 
    

    """
    
    # S1 construct the name 
    
    if (ref==1):
        month='jan'
    elif (ref==4):
        month='apr'
    elif (ref==7):
        month='jul'
    elif (ref==10):
        month='oct'
        
    
    full_flnm=fdesc.construct_filename(XMONTHX=month)
    fdesc.set_filename(full_flnm, ref)
    # S2 read the data
    
    read_aod_file(fdesc, **keywords)
    
    # S3 re check the index
    
    idx=fdesc.get_index(ref)
    
    return idx

def get_aod_data(fdesc, olon, olat, \
                     oyyyy, omm, odd, osec, \
                     do_debug=False,\
                     **keywords):
    
    """
    read aod data at given point 
    Inputs:
    ---------------------------------------
    1. fdesc:<grdfile_io_cl>: file access 
    2. olon: <array, (nob,)>: longitude
    3. olat: <array, (nob,)>: latitude
    4. oyyyy:<array (nob)>: year (for future use)
    --- oyyyy can also be a single value 
    
    5. omm: <int/array (nob)>: month 
    --- omm can also be a single value 
    
    6. odd:<int/array (nob)>: day (for future use)
    7. osec:<float/array (nob)>: seconds   (for future use)
    
    8. keywords:<dict>: extra inputs

    ---Reserved Entries
    ---common_ref:<str/numeric>: the reference shared by observations. 
    --->ax_colnm_lst:<list, t:str>: column names for 'horizontal' (such as sza and aod) axises  
    --->zname: <str>: axis name for vertical (i.e, log10(pressure)) axis
    --->zval: <array, (nz,)>: vertical axis 
    --->'replaces':<dict>: words to be replaced when decoding line 


    Returns:
    ================================================
    1.zval:<array, (nob,nz)>: aod values
    2.aod_pb:<array, (nob,nz)>: aod PDF at the (olon, olat, otime)
    
    """
    
    nob=npy.size(olon)
    # construct ref and ref_set (a set for different references for obs)
    # 
    
    if ('common_ref' in keywords):
        # shared common reference for time (season)
        
        ref=keyword['common_ref']
        ref_set=[ref]
        ref=nob*[ref]
        ref=npy.array(ref)
        
    else:
        
        # get reference 
        
        ref=get_aod_file_ref(oyyyy,omm, odd, **keywords)
        
        if (npy.size(ref)==1):
            ## convert to array 
            # ref_set is for different ref only 
            
            ref_set=[ref]
            ref=nob*[ref]
            ref=npy.array(ref)
        else:
            ## find different refs needed
            
            ref_set=set(ref)
            ref_set=list(ref_set)
    

        
    
    # S2 get aod at (olon, olat, otime (ref))
    
    if (oob.get_ot_type(olon)==oob.ot_array):
        # if olon and olat are the size same array 
        aod_pb=None
        for iref in ref_set: # loop over different reference number 
            
            ## find idx for data in fdesc with ref=iref
        
            idx=fdesc.get_index(iref)
            if (idx==fdesc.mask_val):
                ## load data if it is not in it 
                idx=add_aod_data(fdesc,iref)
            
            ## get the data and grid at found idx
            
            grd=fdesc.grd_lst[idx]
                
            gdata=fdesc.gdata_lst[idx]
            ax_lon=grd['lon']
            ax_lat=grd['lat']

            grd=fdesc.grd_lst[idx]
            ax_z=grd[aod_zname]
            nz=ax_z.ax_np
            sel_z=ax_z[:]
            
            if (aod_pb==None):
                
                zval=npy.zeros([nob, nz], float)
                aod_pb=npy.zeros([nob, nz], float)
                aod_pb[:,:]=fdesc.mask_val
            
                
            
            
            ## slice observation for ref==iref 
            
            sel_oidx=npy.where(ref==iref)
            sel_olat=olat[sel_oidx]
            sel_olon=olon[sel_oidx]

            ## allocated the subset 

            plon=ax_lon.get_closest_point(sel_olon, fdesc.mask_val)
            plat=ax_lat.get_closest_point(sel_olat, fdesc.mask_val)
            
            pidx=range(npy.size(sel_olon))

            ## get aod pb at (plon, plat)
            
            sel_pb=gdata[plon[pidx], plat[pidx], :]
            
            
            ##  fill back to aod_pb 
            
            aod_pb[sel_oidx[:],:]=sel_pb[:,:]
            
            ## set zval
            
            zval[sel_oidx[:], :]=sel_z[npy.newaxis, :]
        #l#l ref_set
            
            if (do_debug):
            
                rlon=ax_lon[:]
                rlat=ax_lat[:]
                print sel_olon
                print rlon[plon]
                
                print sel_olat
                print rlat[plat]
                sel_aod_lvl=def_aod_lvl
                gdata_aod=gdata[:,:,:]*sel_aod_lvl[npy.newaxis, npy.newaxis, :]
                gdata_aod=npy.sum(gdata, axis=2)
                

                gpl.plot_map(gdata_aod, rlon, rlat, use_pcolor=1)
                plb.show()

            
            
        
    else:
        ## if olon and olat are single values
        
        iref=ref[0]
        ### find index for data set in fdesc 
        idx=fdesc.get_index(iref)

        
        if (idx==fdesc.mask_val):
            idx=add_aod_data(fdesc, iref)
        
        ### fetch the grid and data from fdesc with index=idx
        
        grd=fdesc.grd_lst[idx]
        gdata=fdesc.gdata_lst[idx]
        ax_lon=grd['lon']
        ax_lat=grd['lat']
        
        ax_z=grd[aod_zname]
        nz=ax_z.ax_np
        sel_z=ax_z[:]
        
        plon=ax_lon.get_closest_point(olon, fdesc.mask_val)
        plat=ax_lat.get_closest_point(olat, fdesc.mask_val)

        
        ### fill the cloud array 
        aod_pb=gdata[plon, plat,:]
        zval=npy.array(sel_z)

        if (do_debug):
            
            rlon=ax_lon[:]
            rlat=ax_lat[:]
            print sel_olon
            print rlon[plon]
                
            print sel_olat
            print rlat[plat]
            sel_aod_lvl=def_aod_lvl
            gdata_aod=gdata[:,:,:]*sel_aod_lvl[npy.newaxis, npy.newaxis, :]
            gdata_aod=npy.sum(gdata, axis=2)
            

            gpl.plot_map(gdata_aod, rlon, rlat, use_pcolor=1)
            plb.show()

    return zval, aod_pb


def close_aod_file(fdesc):
    """
    Close fdesc file 
  
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    Returns:
    1. None
    
    """
    del fdesc
    
    return None 


#===< TESTS >====

if (__name__=='__main__'):
    
    
    datapath='/home/lfeng/local_disk/otool_data/clim_dat/'
    viewtype='aqua'
    viewmode='nadir'
    surface='snow'
    dd=2
    mm=3
    yyyy=2006
    flnm=''
    fdesc=open_aod_file(datapath, flnm,  yyyy, mm, dd)
    


    
    fdesc=read_aod_file(fdesc)
    
    

    print 'nset', fdesc.nset
    
    grd=fdesc.grd_lst[0]
    print 'grd.dims',  grd.dims
    gdata=fdesc.gdata_lst[0]
    
    print 'gdata--shape', npy.shape(gdata)
    
    olon=npy.arange(-180, 180, 60.0)
    olat=npy.arange(-90, 90, 30.0)

    print 'values at omm=3'

    oyyyy=2006
    omm=3
    odd=20
    
    osec=3600.0
    
    
    zval, aod=get_aod_data(fdesc, olon, olat, \
                         oyyyy, omm, odd, osec)
    print 'aod[:,0] for 200603', aod[:,0]
    
    
    
    zval, aod=get_aod_data(fdesc, olon[2], olat[2], \
                               oyyyy, omm, odd, osec)
    
    
    print 'single aod[:] at olon[2], olat[2]', aod[:]
    

    print 'values at omm=6'
    

    oyyyy=2006
    omm=6
    odd=20
    
    osec=3600.0
    
    
    zval, aod=get_aod_data(fdesc, olon, olat, \
                         oyyyy, omm, odd, osec)
    print 'aod[:,0] for 200603', aod[:,0]
    
    
    
    zval, aod=get_aod_data(fdesc, olon[2], olat[2], \
                               oyyyy, omm, odd, osec)
    
    
    print 'single aod at olon[2], olat[2]', aod[:]

    zval, aod=get_aod_data(fdesc, olon[0], olat[0], \
                               oyyyy, omm, odd, osec)
    
    print 'single aod at olon[0], olat[0]', aod[:]
    
    

    
    print 'values at omm=[6,6, 3, 2, 9, 11]'
    
    oyyyy=6*[2006]
    omm=npy.array([6,6, 3, 2, 9, 11])
    odd=6*[20]
    osec=6*[3600.0]
    
    
    
    zval, aod=get_aod_data(fdesc, olon, olat, \
                               oyyyy, omm, odd, osec)
    
    print 'months:', omm
    print 'aod [:,0]',   aod[:,0]
    
    print '-------------------------------------'
    
    print '>> aod prof at [2,:]',   aod[2,:]
    print  ' '
    print '>> zval at [2,:]',   zval[2,:]
    
    print '-------------------------------------'

    print 'aod prof at [0,:]:',   aod[0,:]
    print  ' '
    print 'zval at [0,:]:',   zval[0,:]
    
    

    
    
