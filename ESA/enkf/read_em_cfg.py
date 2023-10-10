"""
Functions for reading  configuration of emission ensemble   

Functions
=================================================

1. open_ens_cfg_file: set up class file_desc_cl for file access
2. read_ens_cfg_file: read ens_cfg   data into class file_desc_cl
3. close_ens_cfg_file: delete class file_desc_cl
4. get_data: get ens_cfg data at given locations

"""

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_descfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio


import numpy as npy 

ems_cfg_column_lst=['doy','year', 'npb', 'co2flnm']
ems_cfg_coltype_lst=['i4', 'i4',  'i4',  'S40']
def_delim=','
def_ems_header_start=0
def_ems_data_start='200'  
# if a string is given, it will be used to find data section beginning with this string 


ems_datapath='../surface/'  
ems_cfg_flnm_format='cfg_em_reg_144.XYYYYX'

ems_cfg_colnm_dict={'doy':'doy', 'year':'year', 'npb':'npb', 'co2flnm':'ems_flnm'}


def translate_column_names(colnm_lst, **keywords):
    """
    translate column names to the names used to access data 
    
    Inputs:
    -------------------------------
    1. colnm_lst:<list, t:str>: list of column names
    2. keywords:<dict>: extra inputs:
    ---reversed keywords:
    ---colnm_dict:<dict>: column name dictionary
    
    Returns:
    --------------------
    1. new_column_lst: <list, t:str>:New column names:
    """
    
    if ('colnm_dict' in keywords):
        colnm_dict=keywords[colnm_dict]
    else:
        colnm_dict=ems_colnm_dict
        
    

    if (oob.get_ot_type(colnm_lst)==oob.ot_list):
        ## if the names are in a list 
        new_colnm_lst=list(colnm_lst)
        nl=len(new_colnm_lst)
        for il in range(nl):
            keyname=new_colnm_lst[il]
            if (keyname in colnm_dict):
                new_colnm_lst[il]=colnm_dict[keyname]
    else:
        ## single values
        
        new_colnm_lst=colnm_lst
        keyname=new_colnm_lst
        
        if (keyname in colnm_dict):
            new_colnm_lst=colnm_dict[keyname]
            
    return new_colnm_lst


        
def get_ems_cfg_file_ref(yyyy, mm, dd, **keywords):
    """ check file reference
    
    for emsemble configuration. By default 
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


def open_ems_cfg_file(flnm, \
                          datapath, \
                          yyyy, \
                          mm, \
                          dd,\
                          **fio_keywords):
    
    """ construct fdesc to read ensemble  configurations data
    
    
    Inputs:
    ==============================================
    1. flnm:<str>: filename (reserved)
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

    
    # construct file_desc_cl
    
    ## which columns includes the file 

    if ('fl_colnm_lst' in fio_keywords):
        fl_colnm_lst=fio_keywords['fl_colnm_lst']
    else:
        ###c using the default one 
        fl_colnm_lst=ems_cfg_column_lst
    
        
    ## separator for the 

    if ('delim' in fio_keywords):
        delim=fio_keywords['delim']
    else:
        ###c using the default one 
        
        delim=def_delim
    
    
    ## data start 
    if ('data_start' in fio_keywords):
        data_start=fio_keywords['data_start']
    else:
        ###c using the default one 
        
        data_start=def_ems_data_start

    
        
    if ('header_start' in fio_keywords):
        header_start=fio_keywords['header_start']
    else:
        ###c using the default one 
        
        header_start=def_ems_header_start
        
    
    # ##c: name translation 
    
    colnm_dict=ems_cfg_colnm_dict
    
    if ('colnm_dict' in fio_keywords):
        colnm_dict=keywords['column_dict']

        
    # S3 construct fdesc
        
    ref=get_ems_cfg_file_ref(yyyy, mm, dd, **fio_keywords)

    
    
    fdesc=ofo.descfile_desc_cl(flnm, ref,\
                                  header_start, \
                                  data_start, \
                                  delim=delim, \
                                  colnm_lst=fl_colnm_lst,\
                                  colnm_dict=colnm_dict)
    
    
    # S4 set format 
    if (len(flnm)>0):
        flnm_format=flnm
    else:
        flnm_format=ems_cfg_flnm_format

    fdesc.set_filename_format(flnm_format)
    
    # S5 set path 
    
    fdesc.set_file_path(datapath)

    # S6 construct name 
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy)
    
    
    
    # S7 set file name 
    
    fdesc.set_filename(full_flnm, ref)
    
    
    return fdesc


def add_ems_cfg_data(fdesc, ref, **keywords):
    
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
    
    fdesc=read_ems_cfg_file(fdesc, **keywords)
    
    # S3 re check the index
    
    idx=fdesc.get_index(ref)
    
    return idx

def read_ems_cfg_file(fdesc, **fio_keywords):
    
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
        sel_colnm_lst=ems_cfg_column_lst 

    
    if ('sel_coltype_lst' in fio_keywords):
    
        sel_coltype_lst=keywords['sel_coltype_lst']
    else:
        # use default 
        sel_coltype_lst=ems_cfg_coltype_lst 
    
    
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



def close_ens_cfg_file(fdesc):
    
    """ Close fdesc file 
    
    Inputs:
    --------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access
    
    Returns:
    1. None
    
    """
    
    del fdesc
    
    return None

def get_ems_cfg_data(fdesc,  oyyyy, omm, odd,\
                         **keywords):
    """
    read enr_cfg data at given point 
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
    
    1. cur_ems_cfg:<recarray>: data for configuations 
    
    
    """
    
    ## get reference 
    
    iref=get_ems_cfg_file_ref(oyyyy,omm, odd, **keywords)
    
    ### find index for data set in fdesc 
    
    idx=fdesc.get_index(iref)
    
    if (idx==fdesc.mask_val):
        idx=add_ems_cfg_data(fdesc, iref)
        
    ### fetch the grid and data from fdesc with index=idx
        
    cur_ems_cfg=fdesc.data_lst[idx]
    
        
    
    return cur_ems_cfg


#===< TESTS >====

if (__name__=='__main__'):
    yyyy=2009
    mm=12
    dd=1
    
    datapath='../surface/'
    flnm=''
    fdesc=open_ems_cfg_file(flnm, datapath, yyyy, mm, dd)
    
    fdesc=read_ems_cfg_file(fdesc)
    
    print 'nset', fdesc.nset
    data=fdesc.data_lst[0]
    
    print 'data--shape', npy.shape(data)
    
    data=get_ems_cfg_data(fdesc, yyyy, mm, dd)
    
    print data['doy']
    print data['ems_flnm']
    
    

    
