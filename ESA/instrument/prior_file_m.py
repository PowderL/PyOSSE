"""
    
     Functions to access apriori data, which is assumed to be used for retrieval 
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.21
    History: v0.95, 2013.02.21
    
    
    Functions
    =================================================
    
    1. get_prior_column_name: get column names in prior files 
    2. open_apr_file: set up class file_desc_cl for file access
    3. read_apr_file: read apr data into class file_desc_cl
    4. close_apr_file: delete class file_desc_cl
    5. get_apr_data: get apriori data at given location and date
    
    """

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_txtfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.vertical_interp_m as vitpl_m



#====< PARAMETERS >====
##1 columns in the text file 
apr_column_lst=['pres', 'vmr']

### dictionary for translating variable names to those stored in the file

apr_varname_dict={'pres':'pres', 'vmr':'vmr'}

apr_filename_format='clim_XGPNAMEX'

##2 starting row number of data section  
apr_data_start=0

##3 starting row number of header section
apr_header_start=0


#<<< FUNCTIONS >>> 

def get_prior_column_name():
    
    """get column names in prior files 
    
    Returns
    ---------------------------------------
    1. apr_colnm_lst:<list, str>: name of columns in file 
    """
    
    return apr_column_lst


def open_prior_file(flnm, datapath, \
                        gpname,\
                        yyyy, \
                        mm, \
                        dd,\
                        **fio_keywords):
    
    """ open file to read prior data 

    Inputs:
    -----------------------------------------------------------
    1. flnm:<str>:   name of file (reserved for future use )
    2. datapath:<str>: file path
    3. gpname:<str>: name of the tracer
    4. yyyy: <int>: year
    5. mm: <int>: month
    6. dd:<int>: day 
    7. fio_keywords:<dict>: extra inputs. 
    
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

    # S1: construct file_desc_cl
    
    # #c: names for columns in the text file 
    
    if ('fl_colnm_lst' in fio_keywords):
        
        fl_colnm_lst=fio_keywords['fl_colnm_lst']
    else:
        ###c using the default one 
        fl_colnm_lst=apr_column_lst
    
    # #c: delim for dividing items 
        
    if ('delim' in fio_keywords):
        delim=fio_keywords['delim']
    else:
        delim=' '
    
    # #c: start line and end line
    
    
    header_start=apr_header_start
    data_start=apr_data_start
    
    if ('header_start' in fio_keywords):
        header_start=fio_keywords['header_start']
    

    if ('data_start' in fio_keywords):
        data_start=fio_keywords['data_start']
    
    # #c: create fdesc
    
    fdesc=ofo.file_desc_cl(header_start, \
                               data_start, \
                               delim=delim, \
                               colnm_lst=fl_colnm_lst)
    
    # S2: set file format 

    fdesc.set_filename_format(apr_filename_format)

    # S3: set file path 
    fdesc.set_file_path(datapath)

    # S4: construct name 
    
    doy=tm.day_of_year(yyyy, mm, dd)
    
    
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy, \
                                           XMMX=mm, XDDX=dd, \
                                           XDOYX=doy, XGPNAMEX=gpname)
    
   
    # set file name 
    
    fdesc.set_attr('filename', full_flnm)
    
    
    return fdesc

def read_prior_file(fdesc, **fio_keywords):
    
    """
    Read prior data into fdesc.detable
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    2. yyyy, mm, dd:<int>: year month day
    3. fio_keywords:<dict>: extra inputs for file read write
    
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
        sel_colnm_lst=apr_column_lst 
    
    
    flnm=fdesc.get_attr('filename')
    sel_colno_lst=fdesc.get_colno_list(sel_colnm_lst)
    fdesc.read_table_from_file(flnm, sel_colnm_lst, sel_colno_lst)
    

    return fdesc



def close_prior_file(fdesc):
    
    """ Close fdesc file 
  
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    Returns:
    1. None
    
    """
    
    del fdesc
    
    return None


def get_prior_data(fdesc, **fio_keywords):
    
    """
    
    get prior data at observation locations
    
    
    Inputs:
    ---------------------------------------------------------------
    1. fdesc:<fdesc_cl>: class for file access
    2. fio_keywords :<dict>: extra parameters for file access 
    
    Returns:
    -------------------------------------
    1. dtable:<recarray>: origin prior data  set 
    
    """
    
    return fdesc.dtable





#===< TESTS >====

if (__name__=='__main__'):
    datapath='/home/lfeng/local_disk/otool_data/clim_dat/'
    doy=205
    yyyy=2006
    flnm=''
    gpname='CO2'
    mm=3
    dd=6
    
    fdesc=open_prior_file(flnm,  datapath, gpname, yyyy, mm,dd)
    nlines=read_prior_file(fdesc)
    close_prior_file(fdesc)
    
    pres=fdesc.get_column_data('pres')
    co2=fdesc.get_column_data('vmr')
    
    print pres
    print co2
    
    
    
    
    



