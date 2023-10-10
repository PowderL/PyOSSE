"""
Functions for creating and  reading  configuration for ensemble runs

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15

Functions
=================================================

1. open_enr_desc_file: set up class file_desc_cl for file access
2. read_enr_desc_file: read enr_desc   data into class file_desc_cl
3. close_enr_desc_file: delete class file_desc_cl
4. get_enr_desc_table: get enr_desc_table_lst at given time
5. write_enr_desc_colum_heads: write out the column names etc
6. write_enr_desc_colum_val: write out the column values
7. open_enr_desc_file_write: open one new run description file 
8. close_new_enr_desc_file: close a new run description file


"""

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_descfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio

#SECTION (configure outputs)
#========< PARAMETERS >======

# S1:  Columns in ensemble run description text file

enr_desc_column_lst=['step', 'mem_st',  'mem_end',  \
                        'year_st',  'year_end',  \
                        'doy_st',  'doy_end',  \
                        'flnm', 'output_name', \
                        'index_method']


enr_desc_coltype_lst=['i4', 'i4',  'i4',  \
                         'i4',  'i4',  \
                         'i4',  'i4',  \
                         'S40', 'S40', \
                         'i4']

# S2: dictionary for tranlating column names to names used to access data by the code 


enr_desc_colnm_dict={'step':'step', \
                        'mem_st':'mem_st', \
                        'mem_end': 'mem_end',\
                        'year_st':'year_st', \
                        'year_end':'year_end',  \
                        'doy_st':'doy_st',  \
                        'doy_end':'doy_end', \
                        'flnm':'flnm_ems', \
                        'output_name':'flnm_mod'}




def_delim=','
def_enr_header_start=0
def_enr_data_start='000'  # if a string is given, it will be used to find data section beginning with this string 


enr_datapath='/XYYYYX/'  
enr_desc_flnm_format='enr_desc_XYYYYX.dat'


        
def get_enr_desc_file_ref(yyyy, mm, dd, **keywords):
    """ check file reference
    
    for enremble configuration. By default 
    the reference number is year 
    
    Inputs:
    ----------------------------------------------------
    1. yyyy:<array, (ntime,)>: year
    2. mm: <array, (ntime,)>: month
    3. dd:<array, (ntime,)>: dd
    4. keywords:<dict>: extra information
    
    Returns:
    ------------------------------------------------------
    1. ref:<int>: reference number 
    --- by default, the reference number is year
        
    """

    if (oob.get_ot_type(yyyy)==oob.ot_array):
        ref=npy.array(yyyy)
        return ref
    else:
        ref=yyyy
    return ref 

def open_enr_desc_write(flnm, \
                            datapath, \
                            yyyy, \
                            mm, \
                            dd,\
                            **fio_keywords):
    full_flnm=ovio.construct_filename(datapath, flnm, \
                                          XYYYYX=yyyy, XMMX=mm, XDDX=dd)
    fl=open(full_flnm, 'w')
    return fl


def open_enr_desc_file(flnm, \
                           datapath, \
                           yyyy, \
                           mm, \
                           dd,\
                           **fio_keywords):
    
    """ construct fdesc to read ensemble run configurations data
    
    
    Inputs:
    ==============================================
    1. flnm:<str>: format for filename
    2. datapath:<str>: file path
    3. yyyy:<array, (ntime,)>: year
    4. mm: <array, (ntime,)>: month
    5. dd:<array, (ntime,)>: dd
    6. fio_keywords:<dict>: extra inputs 
    ---Reserved keywords:
    --->fl_colnm_lst: <list>:  list of all column names in the file 
    --->sel_colnm_lst: <list>: list of the names for column to be read
    --->colnm_dict:<dict>: dictionary for translating column names
    --->delim: <str>: separator. 
    --->data_start:<int>: starting line number of the data section 
    --->header_start:<int>: line number of header section
    
    
    Returns:
    =============================================
    1. fdesc:<file_desc_cl>: class for file access
    
    """

    # S1 open the file for entry table

    
    
    # #T1:  which columns includes the file 

    if ('fl_colnm_lst' in fio_keywords):
        fl_colnm_lst=fio_keywords['fl_colnm_lst']
    else:
        ###c using the default one 
        fl_colnm_lst=enr_desc_column_lst
    
        
    # #T2: separator for the 

    if ('delim' in fio_keywords):
        delim=fio_keywords['delim']
    else:
        # #c:  using the default one 
        
        delim=def_delim
    
    
    # #T3:  data start 
    
    if ('data_start' in fio_keywords):
        data_start=fio_keywords['data_start']
    else:
        # #c: using the default one 
        
        data_start=def_enr_data_start

    # #T4: head start

        
    if ('header_start' in fio_keywords):
        header_start=fio_keywords['header_start']
    else:
        # #c: using the default one 
        
        header_start=def_enr_header_start
    
    # #T5: name translation 
    
    colnm_dict=enr_desc_colnm_dict
    
    if ('colnm_dict' in fio_keywords):
        colnm_dict=keywords['column_dict']
        
    # S3: construct fdesc
        
    ref=get_enr_desc_file_ref(yyyy, mm, dd, **fio_keywords)
    
    fdesc=ofo.descfile_desc_cl(flnm, ref,\
                                   header_start, \
                                   data_start, \
                                   delim=delim, \
                                   colnm_lst=fl_colnm_lst, \
                                   colnm_dict=colnm_dict)
    
    # S4: set format 
    
    if (len(flnm)>0):
        flnm_format=flnm
    else:
        flnm_format=enr_desc_flnm_format
        
    fdesc.set_filename_format(flnm_format)
    
    # S5:  set path 
    
    fdesc.set_file_path(datapath)
    
    # S6: construct name 
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy)
    
    
    
    # S7: set file name 
    
    fdesc.set_filename(full_flnm, ref)
    
    
    return fdesc


def add_enr_desc_data(fdesc, ref, **keywords):
    
    """ add table for ensemble configurations  to data_lst in fdesc
    
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
    
    
    full_flnm=fdesc.construct_filename(XYYYYX=ref)
    fdesc.set_filename(full_flnm, ref)
    
    # S2 read the data
    
    fdesc=read_enr_desc_file(fdesc, **keywords)
    
    # S3 re check the index
    
    idx=fdesc.get_index(ref)
    
    return idx

def read_enr_desc_file(fdesc, **fio_keywords):
    
    """
    Read ensemble configurations as table into fdesc.data_lst
        
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    2. fio_keywords:<dict>: extra inputs for file read write
    ---Reserved words:
    --->fl_colnm_lst: list of all column names in the file 
    --->sel_colnm_lst: list of the names for column to be read
    --->colnm_dict:<dict>: dictionary for translating column names
    
    Returns:
    -------------------------------------
    1. ndata:<int>: length of fdesc.dtable 

    
    Notes:
    -------------------------------------
    1. data will be stored as recarray by fdesc.dtable
       

    """
    
    # read into lines
    
    flnm=fdesc.cur_flnm
    print flnm
    
    fl=open(flnm, 'r')
    lines=fl.readlines()
    fl.close()

        
    
    if ('sel_colnm_lst' in fio_keywords):
        sel_colnm_lst=keywords['sel_colnm_lst']
    else:
        # use default 
        sel_colnm_lst=enr_desc_column_lst 

    
    if ('sel_coltype_lst' in fio_keywords):
    
        sel_coltype_lst=keywords['sel_coltype_lst']
    else:
        # use default 
        sel_coltype_lst=enr_desc_coltype_lst 
    
    
    sel_colno_lst=fdesc.get_colno_list(sel_colnm_lst)
    
    # S3 read data into table 
    
    ## change the name of columns in the table 
    new_colnm_lst=fdesc.translate_column_names(sel_colnm_lst)
    
    
    ## read data 

    data=fdesc.read_table_from_lines(lines, \
                                         new_colnm_lst, \
                                         sel_colno_lst, sel_coltype_lst)
    ## append the data
    

    head=fdesc.read_header_from_lines(lines)
    
    fdesc.append_data(fdesc.cur_flnm, data, fdesc.cur_ref, head)
    
    return fdesc



def close_enr_desc_file(fdesc):
    
    """ Close fdesc file 
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    Returns:
    ---------------------------------------------
    1. None
    
    """
    
    del fdesc
    
    return None



def get_enr_desc_table(fdesc,  oyyyy, omm, odd,\
                          **keywords):
    """
    read run description  data table at a given point (such year)
    
    Inputs:
    ---------------------------------------
    1. fdesc:<grdfile_io_cl>: file access 
    2. oyyyy:<int>: year
    3. omm: <int>: month 
    4. odd:<int>: day (for future use)
    
    6. keywords:<dict>: extra inputs
    ---Reserved Entries
    
    Returns:
    
    ================================================
    1. fdesc.data_lst:<list, t:recarray>:class for runs covering given date
    
    """
    
    # S1:  get reference 
    
    iref=get_enr_desc_file_ref(oyyyy,omm, odd, **keywords)
    
    # S2: find index for data set in fdesc 
    
    idx=fdesc.get_index(iref)
    
    if (idx==fdesc.mask_val):
        # S3: if not in the list, add the new data
        
        
        idx=add_enr_desc_data(fdesc, iref, **keywords)
        
    
    desc_table=fdesc.data_lst[idx]
    print len(fdesc.data_lst), idx
    
    return desc_table

def open_enr_desc_file_write(flnm):
    """
    open a new file 
    
    Inputs:
    -----------------------------------------------
    1. flnm:<str>: file name 
    
    Returns:
    -----------------------------------------------
    1. fl:<file>: file handle 
    
    """
    fl=open(flnm, 'w')
    return fl

def close_new_enr_desc_file(fl):
    
    """
    close file 
    Inputs;
    ---------------------------------
    1. fl:<file>: file handle 
    
    
    """
    fl.close()


def write_enr_desc_colum_heads(fl, colnm_lst=enr_desc_column_lst, max_nmlen=20):
    """ write out the column names
    
    Inputs:
    -----------------------------------------------------------
    1. fl:<file>: handle to file 
    2. colnm_lst:<list, t:str>: list of column names
    3. max_nmlen:<int>: maximum length of column name

    Returns:
    ----------------------------------
    1. fl:<file>: file object
    
    """
    
    ncol=len(colnm_lst)
    nm_fmt='%-'+str(max_nmlen)+'s'
    sformat=""
    if (ncol>1):
        sformat=(ncol-1)*(nm_fmt+', ')
    sformat=sformat+nm_fmt+'\n'
    fl.write(sformat, colnm_lst)
    return fl



def write_enr_desc_colum_val(fl, col_vals, \
                                 coltype_lst=enr_desc_coltype_lst):
    
    """ write out the column values
    Inputs:
    --------------------------------------------------------
    1. fl:<file>: handle to run description file
    2. col_vals:<list/array>: value of each column
    
    Returns:
    -----------------------------------------------------------
    1. fl:<file>: file handle
    
    
    """
    
    ncol=len(colval_lst)
    sformat=""
    
    for icol in range(ncol):
        stype=coltype_lst[icol]
        if (stype[0]=='i'):
            cur_fmt='%8d'
        elif (stype[0]=='s'):
            cur_fmt='%-30s'
        elif (stype[0]=='f'):
            cur_fmt='%10.4f' 
        if (icol==ncol-1):
            sformat=sformat+cur_fmt+'\n'
        else:
            sformat=sformat+cur_fmt+','
        
    fl.write(sformat, colval_lst)
    return fl




    
        
#===< TESTS >====

if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    dd=2
    mm=3
    yyyy=2009
    flnm='enr_cfg_XYYYYX.dat'
    # S1: open file 
    
    fdesc=open_enr_desc_file(flnm, datapath, yyyy, mm, dd)
    
    #  S2: read in run description file 
    
    fdesc=read_enr_desc_file(fdesc)
    
    print 'nset', fdesc.nset
    
    gdata=fdesc.data_lst[0]
    
    print 'gdata--shape', npy.shape(gdata)
    
    
    # S3: list
    
    
    enr=get_enr_desc_table(fdesc, yyyy, mm, dd)
    
    
    
    
    
