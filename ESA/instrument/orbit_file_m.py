"""
    Functions to read orbit file 
    Authors: L. Feng, Edinburgh University
    
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.19
    

    Functions
    =================================================
    1. get_orbit_column_name: check names of columns in orbit file 
    2. open_orbit_file: set up class file_desc_cl for file access
    3. read_orbit_file: read orbit data into class file_desc_cl
    4. close_orbit_file: delete class file_desc_cl
    
    5. get_orbit_data: get orbit at given date
    
    """


import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_txtfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm



#===< PARAMETERS >===

##1 columns in the text file 
orb_column_lst=['lat', 'lon', 'number', 'date', 'sza',  'effective']

### dictionary for translateing variable names to those stored in the file

orb_varname_dict={'lon':'lon', 'lat':'lat', 'sza':'sza', \
                      'cnt':'number', 'time':'date'}

##2 starting row number of data section  
orb_data_start=4

##3 starting row number of header section
orb_header_start=3


##4 template for file name 

orb_filename_format='XVIEWTYPEX_num_XVIEWMODEX_XDOYX.dat'

#===< FUNCTIONS >====

def get_orbit_column_name():
    
    """get default column name in the orbit file 
    
    Returns
    ---------------------------------------
    1. orb_colnm_lst:<list, str>: default name of columns in file 
    """
    
    return orb_column_lst


def open_orbit_file(flnm, datapath, \
                        viewtype, viewmode, \
                        yyyy, \
                        mm, \
                        dd,\
                        **fio_keywords):
    
    """ open file to read orbit data 
    
    Inputs:
    -----------------------------------------------------------
    1. flnm:<str>:   name of file (reserved for future use )
    2. datapath:<str>: file path
    3. viewtype:<str>: instrument type 
    4. viewmode:<str>: nadir or glint view
    5. yyyy: <int>: year
    6. mm: <int>: month
    7. dd:<int>: day of year
    8. fio_keywords:<dict>: extra inputs. 
    
     --Reserved words:
     --->fl_colnm_lst: list of all column names in the file 
     --->sel_colnm_lst: list of the names for column to be read
     
     --->delim:<str>: dividing word
     --->header_start:<int>: header start line
     --->data_start:<int>: data start line 

    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    
    """

    # construct file_desc_cl
    
    ## which columns includes the file 

    if ('fl_colnm_lst' in fio_keywords):
        
        fl_colnm_lst=fio_keywords['fl_colnm_lst']
    else:
        ###c using the default one 
        fl_colnm_lst=orb_column_lst
    
    
    # #c: delim for dividing items 
        
    if ('delim' in fio_keywords):
        delim=fio_keywords['delim']
    else:
        delim=' '
    

    # #c: start line and end line
    
    
    header_start=orb_header_start
    data_start=orb_data_start
    
    if ('header_start' in fio_keywords):
        header_start=fio_keywords['header_start']
    

    if ('data_start' in fio_keywords):
        data_start=fio_keywords['data_start']
    
    # #c: create fdesc

    
    fdesc=ofo.file_desc_cl(header_start, data_start, \
                               delim=delim, colnm_lst=fl_colnm_lst)
    
    # set file format 

    fdesc.set_filename_format(orb_filename_format)

    # set file path 
    fdesc.set_file_path(datapath)

    # construct name 
    
    doy=tm.day_of_year(yyyy, mm, dd)
    sdoy=r'%d' %doy
    
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy, XMMX=mm, XDDX=dd, \
                                           XVIEWTYPEX=viewtype, \
                                           XVIEWMODEX=viewmode, XDOYX=sdoy)
    
   
    # set file name 
    
    fdesc.set_attr('filename', full_flnm)
    
    
    return fdesc

def read_orbit_file(fdesc, **fio_keywords):
    
    """
    Read orbit data into fdesc.detable
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    2. fio_keywords:<dict>: extra inputs for file read write
    
     --Reserved words:
     --->fl_colnm_lst: list of all column names in the file 
     --->sel_colnm_lst: list of the names for column to be read
        
    
    
    

    Returns:
    -------------------------------------
    1. ndata:<int>: length of fdesc.dtable 


    Notes:
    -------------------------------------
    1. data will be stored as recarray by fdesc.dtable
        
    """
    
    if ('sel_colnm_lst' in fio_keywords):
        sel_colnm_lst=keywords[sel_colnm_lst]
    else:
        # use default 
        sel_colnm_lst=orb_column_lst 
    
    flnm=fdesc.get_attr('filename')
    sel_colno_lst=fdesc.get_colno_list(sel_colnm_lst)
    
    fdesc.read_table_from_file(flnm, sel_colnm_lst, sel_colno_lst)
    numb=fdesc.dtable['number']
    usd_idx=npy.where(numb>0)
    fdesc.dtable=fdesc.dtable[usd_idx]
    return len(fdesc.dtable)


def close_orbit_file(fdesc):
    
    """ Close fdesc file 
  
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    Returns:
    1. None
    
    """
    
    del fdesc
    
    return None


def get_orbit_data(fdesc,  **fio_keywords):
    """
    get orbit data

    Inputs:
    ---------------------------------------------------------------
    1. fdesc: <fdesc_cl>: class of desc
    2. fio_keywords :<dict>: extra parameters for file access 
    
    Returns:
    -------------------------------------
    1. fdesc.dtable:<recordarray>: table read from orbit file 

    """
    
    return fdesc.dtable



#===< TESTS >====

if (__name__=='__main__'):
    datapath='/home/lfeng/local_disk/otool_data/oco/aqua_0.25x0.25/'
    viewtype='aqua'
    viewmode='nadir'
    doy=205
    yyyy=2006
    flnm=''
    fdesc=open_orbit_file(flnm,  datapath, viewtype, viewmode, yyyy, doy)
    nlines=read_orbit_file(fdesc)
    
    print nlines
    lon=fdesc.get_lon()
    print lon[0:10]

    lat=fdesc.get_lat()
    print lat[0:10]
    
    sza=fdesc.get_column_data('sza')

    close_orbit_file(fdesc)
    
    print sza[0:10]
    
    
    
