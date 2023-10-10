"""
   Functions for reading  Cloud PDF file 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.26
   
   Functions
   =================================================
   
   # file access

   1. get_cld_column_name:  check names of columns in AK file
   2. get_cld_surface_type: check names of surface types
   
   3. cld_rescale: rescale CLD axis 
   
   4. open_cloud_file: set up class file_desc_cl for file access
   5. read_cloud_file: read cloud   data into class file_desc_cl
   6. close_cloud_file: delete class file_desc_cl
   7. get_data: get cloud data at given locations
   
   # file (dataset) manage
   
   8. get_cloud_file_ref: check file reference
   9. add_cloud_data: add cloud data to the list 
   
   """

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_grdfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import pylab as plb



#===< PARAMETERS >======


##1 levels (i.e., percentage)

def_cld_lvl=npy.arange(0.0, 20.0)
def_cld_lvl=0.05*def_cld_lvl


##2 columns in CLOUD file to be read


def_cld_colno_lst=[0,1,2]
def_cld_colnm_lst=['lat', 'lon', 'val']
## translations from var name 
cld_varname_dict={'lat':'lat', 'lon':'lon', 'cloud':'val'}

##3 replacements for special character used in text lines

def_cld_replaces={'-9999.0':'-999.0'} #c Rescaling by a factor of 10 is needed 


##4  axises (i.e., dependences of ak)
 
cld_axis_name=['lon', 'lat']

##5 verical levels

cld_zname=None
cld_zval=None



##6 starting row number of data section  
cld_data_start=0

##7 starting row number of header section
cld_header_start=-1

##8 template for file name 

cld_filename_format='meancloud_frac_XMONTHX_landsea.dat'


#===< FUNCTIONS >====

def get_cld_file_ref(yyyy, mm, dd, **keywords):
    """ check file reference
    for cloud climatology, the reference number is climatology months [1, 4, 7, 10]
        
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
    

        
def get_cld_column_name():
    
    """
    check column names defined in cloud file 
    
    Returns
    ---------------------------------------
    1. fl_colnm_lst:<list, str>: name of columns in file 
    
    """

    fl_colnm_lst=def_cld_colnm_lst

    return fl_colnm_lst




def open_cloud_file(datapath, flnm,  \
                      yyyy, mm,dd, \
                      **keywords):
    
    """ setup fdesc for average kernal data
    
    Inputs:
    -----------------------------------------------------------
    1. flnm:<str>:   name of file (reserved for future use )
    2. datapath:<str>: file path
    3. yyyy: <int>: year
    4. mm:<int>: month
    5. dd: <int>: day 
    6:  keywords:<dict>: extra for reading file
    
    ---Reserved Words:
    --->'fl_colnm_lst': <list, t:str>: names of all 
       columns in the cloud file 
    
    --->'delim':<str>: separator in the columns 
    
    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access

    """

    
     # S1 set up a full column list 
    
    if ('fl_colnm_lst' in keywords):
        fl_colnm_lst=keywords['fl_colmn_lst']
    else:
        fl_colnm_lst=def_cld_colnm_lst 
        
    
    delim=' '
    
    if ('delim' in keywords):
        delim=keywords['delim']
    
    
    # S2 construct fdesc
    
    ref=get_cld_file_ref(yyyy, mm, dd, **keywords)
    
    if (ref==1):
        month='jan'
    
    elif (ref==4):
        month='apr'
    
    elif (ref==7):
        month='jul'
    
    elif (ref==10):
        month='oct'
    
        

    fdesc=ofo.grdfile_desc_cl(flnm, ref, 
                              cld_header_start, \
                                  cld_data_start, \
                                  delim=delim, \
                                  colnm_lst=fl_colnm_lst)
    
    
    # S3 set format 
    
    fdesc.set_filename_format(cld_filename_format)
    
    # S4 set path 
    fdesc.set_file_path(datapath)
    
    # S5 set full name and data register month 
    
    
    
        
    full_flnm=fdesc.construct_filename(XMONTHX=month)
    
    fdesc.set_filename(full_flnm, ref)
    
    return fdesc
 
    
def read_cloud_file(fdesc, **keywords):
    
    
    """
    Read cloud  data into fdesc
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<grdfile_desc_cl>: class for file access
    2. keywords:<dict>: extra for data reading. 
    
    ---Reserved Entries
    --->replaces: <dict, old:new>: texts should be replaced. 
    ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
    ---> zval: <array, (nz,)>: vertical axis 
    ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
    ---> fl_colno_lst:<lst, t:int>: columns to be read 
        
 
    Returns:
    -------------------------------------
    1. fdesc:<grdfile_desc_cl>: file access 
    
    Notes:
    -------------------------------------
    1.fdesc.grd_lst and fdesc.gdata_lst will store 1) grid; and 2) gridded data
    
    
    
    """
    # S1 load full filename
    flnm=fdesc.cur_flnm
    
    
    # S2 read data into maxtrix (fdesc.data)
    
    fl=open(flnm, 'r')
    lines=fl.readlines()
    fl.close()
    
    # S2 load format information 
    ## T1 columns for horizontal axis 
    
    if ('ax_colnm_lst' in keywords):
        ax_colnm_lst=keywords['ax_colnm_lst']
    
    else:
        ax_colnm_lst=cld_axis_name

    
         

    
    ## T2 vertical axis 
    
    if ('zname' in keywords):
        zname=keywords['zname']
    else:
        zname=cld_zname
    
    if ('zval' in keywords):
        zval=keywords['zval']
    else:
        zval=cld_zval

    
    ## T3 replacements for line decoding 
    
    if ('replaces' in keywords):
        replaces=keywords['replaces']
    else:
        replaces=def_cld_replaces
   

    ## columns needs to be read 
    
    if ('fl_colno_lst' in keywords):
        fl_colno_lst=keywords['fl_colno_lst']
    else:
        fl_colno_lst=def_cld_colno_lst
    

    # S3 read data into matrix (array)
    data=fdesc.read_matrix_from_lines(lines, \
                                          colno_lst=fl_colno_lst, \
                                          **replaces)
    
    
    # S4 regrid data to fdesc.gdata
    
    ## T1 get column number for the axis 
    
    ax_colno_lst=fdesc.get_colno_list(ax_colnm_lst)
    
    ## T2 convert array to grid 
    grd, gdata=fdesc.convert_matrix_to_grid(data, ax_colnm_lst, \
                                                ax_colno_lst, \
                                                zname=zname, \
                                                zval=zval)
    
    fdesc.append_data(fdesc.cur_flnm, data, grd, gdata, fdesc.cur_ref)
    
    
    return fdesc




def close_cloud_file(fdesc):
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


def add_cloud_data(fdesc, ref, **keywords):
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
    
    read_cloud_file(fdesc, **keywords)
    
    # S3 re check the index
    
    idx=fdesc.get_index(ref)
    
    return idx

def get_cloud_data(fdesc, olon, olat, \
                       oyyyy, omm, odd, osec, \
                       do_debug=False,\
                       **keywords):
    
    """
    read cloud data at given point 
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

    Returns:
    ================================================
    cld:<array, (nob, )>: cloud at the (olon, olat, otime)
    
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
        
        ref=get_cld_file_ref(oyyyy,omm, odd, **keywords)
        
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
    

        
    
    # S2 get cloud at (olon, olat, otime (ref))
    
    if (oob.get_ot_type(olon)==oob.ot_array):
        # if olon and olat are the size same array 
        cld=npy.zeros(nob, float)
        for iref in ref_set: # loop over different reference number 
            
            ## find idx for data in fdesc with ref=iref
        
            idx=fdesc.get_index(iref)
            if (idx==fdesc.mask_val):
                ## load data if it is not in it 
                idx=add_cloud_data(fdesc,iref)
            
            ## get the data and grid at found idx
            
            grd=fdesc.grd_lst[idx]
                
            gdata=fdesc.gdata_lst[idx]
            ax_lon=grd['lon']
            ax_lat=grd['lat']
            gdata=npy.where(gdata<0.0, 1.0, gdata)
        
            
            ## slice observation for ref==iref 
            
            
            sel_oidx=npy.where(ref==iref)
            
            sel_olat=olat[sel_oidx]
            sel_olon=olon[sel_oidx]

            ## allocated the subset 

            plon=ax_lon.get_closest_point(sel_olon, fdesc.mask_val)
            plat=ax_lat.get_closest_point(sel_olat, fdesc.mask_val)

            if (do_debug):
                rlon=ax_lon[:]
                rlat=ax_lat[:]
                print sel_olon
                print rlon[plon]
                
                print sel_olat
                print rlat[plat]
                
                gpl.plot_map(gdata, rlon, rlat, use_pcolor=1)
                plb.show()

            pidx=range(npy.size(sel_olon))

            ## find cloud 
            sel_cld=gdata[plon[pidx], plat[pidx]]
            
            ##  put back to cld 
            
            cld[sel_oidx[:]]=sel_cld[:]
            
        #l#l ref_set

        
    else:
        ## if olon and olat are single values
        
        iref=ref[0]
        ### find index for data set in fdesc 
        idx=fdesc.get_index(iref)

        
        if (idx==fdesc.mask_val):
            idx=add_cloud_data(fdesc, iref)
        
        ### fetch the grid and data from fdesc with index=idx
        
        grd=fdesc.grd_lst[idx]
        gdata=fdesc.gdata_lst[idx]
        ax_lon=grd['lon']
        ax_lat=grd['lat']
        gdata=npy.where(gdata<0.0, 1.0, gdata)
            
        if (do_debug):
            rlon=ax_lon[:]
            rlat=ax_lat[:]
            print sel_olon
            print rlon[plon]
                
            print sel_olat
            print rlat[plat]
                
            gpl.plot_map(gdata, rlon, rlat, use_pcolor=1)
            plb.show()

        
        plon=ax_lon.get_closest_point(olon, fdesc.mask_val)
        plat=ax_lat.get_closest_point(olat, fdesc.mask_val)
        
        

        ### fill the cloud array 
        cld=gdata[plon, plat]
    
    return cld
    
        

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
    fdesc=open_cloud_file(datapath, flnm,  yyyy, mm, dd)
    


    
    fdesc=read_cloud_file(fdesc)
    
    

    print 'nset', fdesc.nset
    
    grd=fdesc.grd_lst[0]
    print 'grd.dims',  grd.dims
    gdata=fdesc.gdata_lst[0]
    
    print 'gdata--shape', npy.shape(gdata)
    
    olon=npy.arange(-180, 180, 60.0)
    olat=npy.arange(-90, 90, 30.0)

    oyyyy=2006
    omm=3
    odd=20
    
    osec=3600.0
    
    
    cld=get_cloud_data(fdesc, olon, olat, \
                           oyyyy, omm, odd, osec)
    print 'cld for 200603', cld
    
    
    
    cld=get_cloud_data(fdesc, olon[2], olat[2], \
                           oyyyy, omm, odd, osec)
    
    
    print 'single cld', cld
    
    
    oyyyy=6*[2006]
    omm=npy.array([6,6, 3, 2, 9, 11])
    odd=6*[20]
    
    
    cld=get_cloud_data(fdesc, olon, olat, \
                           oyyyy, omm, odd, osec)
    
    print 'months:', omm
    print 'cloud',   cld
    
    
