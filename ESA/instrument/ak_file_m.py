"""
   Reading  averaging kernel (AK) file to grdfile_desc_fl. 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.24
      
   Functions
   =================================================
   
   # file (data) access
   
   1 . aod_rescale: rescale AOD axis 
   2.  open_ak_file: set up class file_desc_cl for file access
   3. read_ak_file: read AK  data into class file_desc_cl
   4. load_ak_file: add AK data into class file desc when necessary 
   5. get_ak: get ak values
   6. close_ak_file: delete class file_desc_cl
   7. get_ak_column_name:  check names of columns in AK file
   
   # file manage 

   8. get_ak_file_ref:     check ref number for AK files
   9. decode_ak_ref:       decode ref number into viewmode and surface type 
   
   10. get_ak_sur_type:     check names of surface types
   11. get_ak_surf_type_id: check id of surface types
   
   
   12. get_ak_viewmode:     check names of viewmode
   13. get_ak_viewmode_id:  check id of viewmode
   
   
   
   """


import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_grdfile_io as ofo
import ESA.util.time_module as tm


#===< PARAMETERS >======


## retrieval pressure levels in hPa 

def_ak_pres=npy.array([  1., 70.0, 100.0,\
                             200.0, 350.0, 450.0,\
                             550.0, 650.0, 750.0, \
                             850.0, 950.0, 1000.])
##1 surface type 

def_ak_surf_type=['soil', 'ocean', 'snow', 'desert', 'conifer']
def_ak_surf_type_dict={'soil':3, 'ocean':0, 'snow':2, 'desert':4, 'conifer':1}

## #c: convert glc2000 type to surface type

ak_lc_stype_lst=['conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'conifer',\
                     'soil',\
                     'soil',\
                     'soil',\
                     'desert',\
                     'ocean',\
                     'snow',\
                     'soil']



##2 view mode

def_ak_viewmode=['nadir', 'glint']
def_ak_viewmode_dict={'nadir':0, 'glint':1}

##3 columns in AK file

def_ak_colnm_lst=['sza', 'aod']

##4 replacements for special character used in text lines

def_ak_replaces={'od':'0.'} #c Rescaling by a factor of 10 is needed 

def_aod_scale=10.0 

##5  axises (i.e., dependences of ak)
 
ak_axis_name=['sza', 'aod']

##6 verical levels

ak_zname='lgp'
ak_zval=npy.log10(def_ak_pres)


##7 starting row number of data section  
ak_data_start=1

##8 starting row number of header section
ak_header_start=0

##9 template for file name 
ak_filename_format='aknorm_XSURFACEX.XEXTX'


#===< FUNCTIONS >====



def get_ak_file_ref(iviewmode,isurf_type, **keywords):
    """ check file reference for averaging kernels 
    
    Inputs:
    ----------------------------------------------------
    1. isurf_type:<int/array>: type of surface 
    ---type includes: 
    2. iviewmode: <int/array>: view model 
    
    3. keywords:<dict>: extra information
    
    Returns:
    ------------------------------------------------------
    1. ref:<int/array>: reference number 

    Notes:
    ---------------------------------
    1, ref_no=1000*iviewmode+isurf_type
    
    """
    
    ref=1000*iviewmode+isurf_type
    return ref
   

def get_ak_viewmode_id(viewmode, **keywords):
    
    """
    Get viewmode ID

    Inputs:
    ---------------------------------------------------------

    1. viewmode:<str>: for view mode 
    2. keywords:<dict> extra inputs

    ---reserved keywords: 
    ---> 'ak_viewmode_dict':<dict>: the view modes for documented averaging kernel. 
    
    
    Returns:
    ------------------------------------------------------------
    1. iviewmode:<Int>: the view model ID
    
    
    """
    
    # S1 loading the dictionary 
    
    if ('ak_viewmode_dict' in keywords):
        ak_viewmode_dict=keywords['ak_viewmode_dict']
    else:
        ak_viewmode_dict=def_ak_viewmode_dict
    
    # S2 find the viewmode ID from dictionary 
        
    
    iviewmode=ak_viewmode_dict[viewmode]
    
    

    return iviewmode


def get_ak_viewmode(iviewmode, **keywords):
    
    """
    Get viewmode.

    Inputs:
    ---------------------------------------------------------

    1. iviewmode:<int>: ID for view mode 
    2. keywords:<dict> extra inputs

    ---reserved keywords: 
    ---> 'ak_viewmode_dict':<dict>: the view modes for documented averaging kernel. 
    
    
    Returns:
    ------------------------------------------------------------
    1. viewmode:<str>: the matched view mode
    
    
    """
    
    # S1 loading the dictionary 
    
    if ('ak_viewmode_dict' in keywords):
        ak_viewmode_dict=keywords['ak_viewmode_dict']
    else:
        ak_viewmode_dict=def_ak_viewmode_dict
    
    # S2 find the match view mode 
    
    for viewmode in ak_viewmode_dict:
        id_viewmode=ak_viewmode_dict[viewmode]
        if (id_viewmode==iviewmode):
            return viewmode
    

    # S3 if no match found, an empty string will be returned


    return ""



def get_ak_surf_type_id(surf_type, **keywords):
    
    """
    Get surface type 
    

    Inputs:
    ---------------------------------------------------------

    1. surf_type:<str>: surface type 
    
    2. keywords:<dict> extra inputs

    ---reserved keywords: 
    ---> 'ak_surf_type_dict':<dict>: surface type for documented averaging kernel. 
    
    
    Returns:
    ------------------------------------------------------------
    1. isurf_type:<int>: the matched surface type ID
    
    
    """
    
    # S1 loading the dictionary 
    
    if ('ak_surf_type_dict' in keywords):
        ak_surf_type_dict=keywords['ak_surf_type_dict']
    else:
        ak_surf_type_dict=def_ak_surf_type_dict

    # S2 find ID for surface type  
    
    isurf_type=ak_surf_type_dict[surf_type]

    # S3 if no match found, an empty string will be returned


    return isurf_type




def get_ak_surf_type(isurf_type, **keywords):
    
    """
    Get surface type 
    

    Inputs:
    ---------------------------------------------------------

    1. iviewmode:<int>: ID for view mode 
    2. keywords:<dict> extra inputs

    ---reserved keywords: 
    ---> 'ak_surf_type_dict':<dict>: surface type for documented averaging kernel. 
    
    
    Returns:
    ------------------------------------------------------------
    1. surf_type:<str>: the matched surface_type
    
    
    """
    
    # S1 loading the dictionary 
    
    if ('ak_surf_type_dict' in keywords):
        ak_surf_type_dict=keywords['ak_surf_type_dict']
    else:
        ak_surf_type_dict=def_ak_surf_type_dict

    # S2 find the match surface type  
    
    for surf_type in ak_surf_type_dict:
        id_surface=ak_surf_type_dict[surf_type]
        if (id_surface==isurf_type):
            return surf_type

    # S3 if no match found, an empty string will be returned


    return ""





def decode_ak_ref(ref, **keywords):

    """ check file reference
    for averaging kernels 
    
    Inputs:
    ----------------------------------------------------
    1. ref:<int>: reference number 
    2. keywords:<dict>: additional inputs 
    ---reserved words
    --->'ak_viewmode_dict':<dict>: dictionary for view mode
    --->'ak_surf_type_dict': <dict>: dictionary for surface type
    
    

    Returns:
    ------------------------------------------------------
    1. surface_type:<str>: surface type 
    2. viewmode:<str>: view mode 
    
    """
    ivewemode=int(ref/1000)
    isurf_type=ref-1000*ivewemode
    
    viewmode= get_ak_viewmode(ivewemode, **keywords)
    surf_type=get_ak_surf_type(isurf_type, **keywords)
    
    return viewmode, surf_type

   


def get_ak_column_name():
    
    """
    check column names defined in AK file 
        
    Returns
    ---------------------------------------
    1. fl_colnm_lst:<list, str>: name of columns in file 
    
    """

    fl_colnm_lst=def_ak_colnm_lst

    for iz in range(npy.size(ak_zval)):
        # append 'Z??' to column name list 
        
        szl=r'Z%d' % iz
        fl_colnm_lst=fl_colnm_lst+[szl]
    
    return fl_colnm_lst



def aod_rescale(grd, aod_axnm='aod', scale_factor=def_aod_scale):
    """
    rescaling aod axis in a grid 
    
    Inputs:
    ----------------------------------------------
    1. grd:<gp_grid_cl>: grid 
    2. aod_axnm:<str>: AOD axis name
    3. scale_factor:<float>: scaling factor
    
    Returns:
    ------------------------------------------
    1. grd:<gp_grid_cl>: grid after axis aod is re-scaled
    
    
    """
    # S1 get axis 
    ax_aod=grd[aod_axnm]
    # S2 re-scaling 
    new_val=scale_factor*ax_aod[:]
    # S3 update axis values

    ax_aod.set_axis(new_val)
    # as a result, aod axis in class grd will be updated 
    
    return grd


def open_ak_file(flnm, datapath, \
                     viewtype, \
                     viewmode, surface, \
                     yyyy, mm, dd,\
                     **keywords):
    
    
    """ setup fdesc for average kernal data
    Inputs:
    -----------------------------------------------------------
    1. flnm:<str>:   name of file (reserved for future use )
    2. datapath:<str>: file path
    3. viewtype:<str>: instrument type 
    4. viewmode:<str>: nadir or glint view
    5. surface:<str>: surface type
    6. yyyy: <int>: year
    7. mm: <int>: month
    8. dd:<int>: day 
    
    9. keywords:<dict>: extra inputs for file reading 
    ----Reserved words
    --->'fl_colnm_lst':<list, t:str>: names of all columns in the
    --->'ak_viewmode_dict':<dict>: dictionary for view mode
    --->'ak_surf_type_dict': <dict>: dictionary for surface type
    --->'delim':<str>: splitor in averaging kernel (text) file
    
    
    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access

    """

        
    # S1 set up a full column list 
    
    if ('fl_colnm_lst' in keywords):
        fl_colnm_lst=keywords['fl_colmn_lst']
    else:
        fl_colnm_lst=def_ak_colnm_lst
        
        for iz in range(npy.size(ak_zval)):
            # append 'Z??' to column name list 
            szl=r'Z%d' % iz
            fl_colnm_lst=fl_colnm_lst+[szl]
    
    # S2 construct fdesc

    ## T1 get ref
    
    isurf_type=get_ak_surf_type_id(surface, **keywords)
    iviewmode=get_ak_viewmode_id(viewmode, **keywords)

    ref= get_ak_file_ref(iviewmode,isurf_type, **keywords)
    
    ## T2

    delim=' '
    
    if ('delim' in keywords):
        delim=keywords['delim']
    
        
    fdesc=ofo.grdfile_desc_cl(flnm, ref,\
                                  ak_header_start, \
                                  ak_data_start, \
                                  delim=delim, \
                                  colnm_lst=fl_colnm_lst)
    
                              
    # S3 set format 
    fdesc.set_filename_format(ak_filename_format)

    
    # S4 set path 
    fdesc.set_file_path(datapath)
    
    # S5 set full name 
    
    if (viewmode=='nadir'):
        ext='dat_1'
    else:
        ext='dat_lg_glint'
    
    full_flnm=fdesc.construct_filename(XSURFACEX=surface, \
                                           XVIEWMODEX=viewmode, \
                                           XEXTX=ext)
    fdesc.set_filename(full_flnm, ref)
    
    
    return fdesc


def add_ak_data(fdesc, ref, **keywords):

    """

    Inputs:
    =====================================
    1. fdesc:<grdfile_io_cl>: class for file access
    2. ref:<str/numeric>: reference for data sets read from file 
    3. keywords:<dict>: extra parameters
    ----Reserved words
    --->'ak_viewmode_dict':<dict>: dictionary for view mode
    --->'ak_surf_type_dict': <dict>: dictionary for surface type
    ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
    ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
    ---> zval: <array, (nz,)>: vertical axis 
    --->'replaces': <dict>: words to be replaced when decoding line 
    

    Returns:
    =====================================
    1. idx:<int>: index of the data set with reference =ref. 
    
    """
    # S1 check the view mode etc 
    
    viewmode, surface=decode_ak_ref(ref, **keywords)
    
    if (viewmode=='nadir'):
        ext='dat_1'
    else:
        ext='dat_lg_glint'
    
    full_flnm=fdesc.construct_filename(XSURFACEX=surface, \
                                           XVIEWMODEX=viewmode, \
                                           XEXTX=ext)
    fdesc.set_filename(full_flnm, ref)
    
    read_ak_file(fdesc, **keywords)
    idx=fdesc.get_index(ref)
    
    return idx
    
    
def read_ak_file(fdesc, \
                     **keywords):
    

    
    """
    Read averaging kernel data into fdesc
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<grdfile_desc_cl>: class for file access
    2. keywords:<dict>: extra for data reading. 
    
    ---Reserved Entries
    --->replaces: <dict, old:new>: texts should be replaced. 
    ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
    ---> zval: <array, (nz,)>: vertical axis 
    ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
    
    Returns:
    -------------------------------------
    1. fdesc:<grdfile_desc_cl>:class for file access  


    Notes:
    -------------------------------------
    1.fdesc.grd_lst and fdesc.gdata_lst will store 1) grid; and 2) gridded data
   
    
    
    """
    # S1 load full filename
    
    flnm=fdesc.cur_flnm
    fl=open(flnm, 'r')
    lines=fl.readlines()
    fl.close()
    
    
    # S2 read data into maxtrix (fdesc.data)
    

     ## T1 columns for horizontal axis 
    
    if ('ax_colnm_lst' in keywords):
        ax_colnm_lst=keywords['ax_colnm_lst']
    
    else:
        ax_colnm_lst=ak_axis_name
    
    
    ## T2 vertical axis 
    if ('zname' in keywords):
        zname=keywords['zname']
    else:
        zname=ak_zname
    
    if ('zval' in keywords):
        zval=keywords['zval']
    else:
        zval=ak_zval

    ## T3 replacements for line decoding 
    
    if ('replaces' in keywords):
        replaces=keywords['replaces']
    else:
        replaces=def_ak_replaces
        
    
    
    data=fdesc.read_matrix_from_lines(lines, **replaces)

    
    # S3 regrid data to fdesc.gdata
    
    ax_colno_lst=fdesc.get_colno_list(ax_colnm_lst)
        
    grd, gdata=fdesc.convert_matrix_to_grid(data, ax_colnm_lst, \
                                                ax_colno_lst, \
                                                zname=zname, \
                                                zval=zval)
    # S4 rescale axis aod of fdesc.grd
    
    grd=aod_rescale(grd)
    
    fdesc.append_data(fdesc.cur_flnm, data, grd, gdata, fdesc.cur_ref)

    
    return fdesc 




def close_ak_file(fdesc):
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


def get_ak_data(fdesc, osza, oaod, \
                    iviewmode, \
                    isurf_type,\
                    **keywords):
    
    """
    Read cloud data at given point 
    
    Inputs:
    ----------------------------------------------
    1. fdesc:<grdfile_io_cl>: file access 
    2. osza:<array, (nob,)>: solar zenith angles of observation
    3. oaod:<array, (nob,)>: aerosol optical depth 
    
    4. iviewmode:<int/array>: view mode ID (see Note 1)
    5. isurf_type:<int/array>:  surface type ID (see Note 1)
    6. keywords:<dict>: extra inputs
    
    ---Reserved Entries
    ---common_ref:<str/numeric>: the reference shared by observations. 
    --->'ak_viewmode_dict':<dict>: dictionary for view mode
    --->'ak_surf_type_dict': <dict>: dictionary for surface type
    ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
    ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
    ---> zval: <array, (nz,)>: vertical axis 
    --->'replaces': <dict>: words to be replaced when decoding line 
    
    Returns:
    ------------------------------------------------
    1. zval:<array, (nob, nz)>: vertical coordinate (logp)
    2. ak:<array, (nob, nz)>: averaging kernel at the (osza, oaod)
    
    

    Notes:
    ---------------------------------------------
    1. if iviewmode and isurf_type are either given as a single integer or  
    an array of the same size of osza
    
    """
   
    if (oob.get_ot_type(isurf_type)==oob.ot_list):
        isurf_type=npy.array(isurf_type)
        iviewmode=npy.array(iviewmode)
    
    nob=npy.size(osza)
    
    # construct ref and ref_set (a set for different references for obs)
    # 
    
    
    if ('common_ref' in keywords):
        # shared common reference
        
        ref=keyword['common_ref']
        ref_set=[ref]
        ref=nob*[ref]
        ref=npy.array(ref)
        
    else:
        
        # get reference 
        
        
        if (oob.get_ot_type(iviewmode)<>oob.ot_array):
            # single value 

            ref=get_ak_file_ref(iviewmode, isurf_type)
            ref_set=[ref]
            ref=nob*[ref]
            ref=npy.array(ref)
            
        else:
            ref=get_ak_file_ref(iviewmode, isurf_type)
            ## find different refs needed
            ref_set=set(ref)
            ref_set=list(ref_set)
    
    
        
    
    # S2 get cloud at (olon, olat, otime (ref))
    
    
    if (oob.get_ot_type(osza)==oob.ot_array):
        
        # if olon and olat are the size same array 
        
        ak=None
        
        
        for iref in ref_set: # loop over different reference number 
            
            ## find idx for data in fdesc with ref=iref
        
            idx=fdesc.get_index(iref)
            if (idx==fdesc.mask_val):
                ## load data if it is not in it 
                idx=add_ak_data(fdesc,iref, **keywords)
            
            ## get the data and grid at found idx
            
            grd=fdesc.grd_lst[idx]
            ax_z=grd[ak_zname]
            nz=ax_z.ax_np
            ## zvalue 
            
            sel_z=ax_z[:]
            
            
            if (ak==None):
                zval=npy.zeros([nob, nz], float)
                ak=npy.zeros([nob, nz], float)
            
            
            gdata=fdesc.gdata_lst[idx]
            
            ax_aod=grd['aod']
            ax_sza=grd['sza']
            
            ## slice observation for ref==iref 
            
            sel_oidx=npy.where(ref==iref)

            sel_aod=oaod[sel_oidx]
            sel_sza=osza[sel_oidx]

            ## allocated the subset 
            
            psza=ax_sza.get_closest_point(sel_sza, fdesc.mask_val)
            paod=ax_aod.get_closest_point(sel_aod, fdesc.mask_val)
            
            pidx=range(npy.size(sel_aod))
            
            ## find ak
            
            sel_ak=gdata[psza[pidx], paod[pidx], :]
            
            ##  put back to ak
            
            ak[sel_oidx[:],:]=sel_ak[:,:]
            
            ## set zval
            
            zval[sel_oidx[:], :]=sel_z[npy.newaxis, :]
            
            
        #l#l ref_set

        
    else:
        ## if osza and oaod are single values
        
        iref=ref[0]
        
        ### find index for data set in fdesc 
        
        idx=fdesc.get_index(iref)
        
        
        if (idx==fdesc.mask_val):
            idx=add_ak_data(fdesc, iref, **keywords)
        
        ### fetch the grid and data from fdesc with index=idx
        
        grd=fdesc.grd_lst[idx]
        
        ax_z=grd[ak_zname]
        
        nz=ax_z.ax_np
        sel_z=ax_z[:]
            
        gdata=fdesc.gdata_lst[idx]
        ax_sza=grd['sza']
        ax_aod=grd['aod']
        
        psza=ax_sza.get_closest_point(osza, fdesc.mask_val)
        paod=ax_aod.get_closest_point(oaod, fdesc.mask_val)
        
        ### fill the cloud array 
        
        ak=gdata[psza, paod, :]
        zval=npy.array(sel_z)
    
    return zval, ak




    


#===< TESTS >====

if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/oco/oco_ak/'
    viewtype='aqua'
    viewmode='nadir'
    surface='snow'
    
    mm=3
    dd=2
    yyyy=2006
    flnm=''
    fdesc=open_ak_file(flnm,  datapath, viewtype, \
                           viewmode, surface, yyyy, mm, dd)

    
    fdesc=read_ak_file(fdesc)
    
    grd=fdesc.grd_lst[0]
    print grd.dims
    
    gdata=fdesc.gdata_lst[0]
    
    print npy.shape(gdata)
    
    ax_z=grd[ak_zname]
    
    print 'lgp', ax_z[:]
    

    ax_aod=grd['aod']
     
    print 'aod', ax_aod[:]
    
    
    ax_sza=grd['sza']
    
    print 'sza', ax_sza[:]
    
    # sza, 100*aod
    rpos=[15, 0.25]
    prof=grd.get_profile(rpos, gdata, [0,1])
    
    print prof
    
    print 'try to get data from two observations '
    
    osza=npy.array([10.0,15.0])
    oaod=npy.array([0.2, 0.10])
    iviewmode=[1,0]
    isurf_type=[1,0]
    
    zval, ak=get_ak_data(fdesc, osza, oaod, \
                             iviewmode, \
                             isurf_type)
    

    
    print zval
    print ak


    
    print 'try to get data for one observation'
    
    osza=osza[1]
    oaod=[1]
    iviewmode=iviewmode[1]
    isurf_type=isurf_type[1]
    
    zval, ak=get_ak_data(fdesc, osza, oaod, \
                             iviewmode, \
                             isurf_type)
    
    print zval
    print ak
    

