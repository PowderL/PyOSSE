"""
   Functions for reading  GEOS-Chem time serie file from single tracer run or from ensemble runs 
   
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.26
   
   Functions
   =================================================
   # ts file
   
   1. open_gc_file:   setup fdesc for GEOS-Chem outputs
   2. setup_daily_field: set up  model fields from time-series files
   3. setup_daily_profile: set up  model profiles from time-series files
   
   # ensemble runs

   4. open_enr_file:  setup fdesc for GEOS-Chem ensemble runs
   5. setup_enr_daily_field: set up  model field ensemble from time-series outputs of ensemble runs
   6. setup_enr_daily_profile: set up  model profile ensemble from time-series outputs of ensemble runs
   
   # data access 

   7. get_mod_gp:  get model profile (field)for selected variables
   8. get_mod_sp:  get surface pressure
   9. get_data: find data according to criteria of traname, category and tracer ID
   10. get_mod_pres: computer or fetch model pressure
   
   """

import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_gcfile_io as ofo
import ESA.util.time_module as tm


#===< PARAMETERS >======

# #1: default file names for single tracer (forecast run)

def_ts_filename_format='ts_satelliteXEXTX.XYYYYXXMMXXDDX.bpch'
def_ts_ftracerinfo_format='tracerinfoXEXTX.dat'
def_ts_fdiaginfo_format='diaginfoXEXTX.dat'


# #2: default file names for ensemble runs 

def_enr_filename_format='ts_satellite.STXENRSTEPX.ENXEMSTX-ENXEMENDX.XYYYYXXMMXXDDX.bpch'
def_enr_ftracerinfo_format='tracerinfo.STXENRSTEPX.ENXEMSTX-ENXEMENDX.dat'
def_enr_fdiaginfo_format='diaginfo.STXENRSTEPX.ENXEMSTX-ENXEMENDX.dat'


def open_ts_file(datapath, flnm,  \
                     ftracerinfo,\
                     fdiaginfo,\
                     yyyy, mm,dd, \
                     categorys=None,\
                     tracers=None,\
                     taus=None,\
                     tranames=None,\
                     **fio_keywords):
    
    
    """ setup fdesc for GEOS-Chem outputs
    
    Inputs:
    --------------------------------------------------------------

    1. datapath:<str>: path 
    2. flnm:<str>: file name
    3, ftracerinfo:<str>: name for tracerinfo file
    4. fdiaginfo:<str>: name for diaginfo file

    5. categorys:<list, t:str>: list of categorys to be read
    6. tracers:<list,t:int>: list of tracers to be read 
    7. taus:<list, t:float>: list of taus to be read 
    8. tranames:<list, t:str>:list of tracer name
    9. fio_keywords: <dict>: extra inputs for file IO
    ---Reserved keywords:
    --
        
    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access

    """

    # S1 read in the list 
    
    
    fdesc=ofo.gcfile_desc_cl(flnm, ftracerinfo, \
                                 fdiaginfo,\
                                 yyyy, mm, dd,\
                                 is_enr_file=False,\
                                 categorys=categorys,\
                                 tracers=tracers,\
                                 taus=taus,\
                                 tranames=tranames,\
                                 mask_val=oob.fill_val,\
                                 fio_keywords=fio_keywords)
    
    # S2: set file format: 
    # ##T1: filename
    
    fdesc.set_file_path(datapath)
    
    if (len(flnm)>0):
        fdesc.set_filename_format(flnm)
        
    else:
        fdesc.set_filename_format(def_ts_filename_format)
    
    # ##T2: diaginfo 
      
    fdesc.set_diaginfo_file_path(datapath)
    fdiaginfo=fdiaginfo.strip()
    
    if (len(fdiaginfo)>0):
        fdesc.set_diaginfo_filename_format(fdiaginfo)
    else:
        fdesc.set_diaginfo_filename_format(def_ts_fdiaginfo_format)
    
    
    # ##T3: tracerinfo 
      
    fdesc.set_tracerinfo_file_path(datapath)
    ftracerinf=ftracerinfo.strip()
    
    if (len(ftracerinfo)>0):
        fdesc.set_tracerinfo_filename_format(ftracerinfo)
    else:
        fdesc.set_tracerinfo_filename_format(def_ts_ftracerinfo_format)
        
    
    
    # S4 set path 
    
    return fdesc
 

def setup_daily_profile(fdesc, yyyy, mm, dd, \
                          olon, olat,\
                           ext="",\
                            **keywords):
    
    
    """ set up  the model profiles from time-series files
    Inputs:
    ----------------------------------------------------------------------------------
    1. yyyy, mm,dd:<int>:  year, month, day
    2. olat, olon:<array, (nobs)>:  obsvervation location
    3. ext:<str>: extension for filname
    8. keywords:<dict>: extra inputs:
    
    
    """
    
    prof_lst=fdesc.read_file_to_profile_list(yyyy, mm, dd, \
                                                 olon, olat,\
                                                 ext=ext,\
                                                 **keywords)
        
    
    
    
    return fdesc


def setup_daily_field(fdesc, yyyy, mm, dd, \
                          olon, olat,\
                           ext="",\
                            **keywords):
    
    
    """ set up  the model field from time-series files
    Inputs:
    ----------------------------------------------------------------------------------
    1. yyyy, mm,dd:<int>:  year, month, day
    2. olat, olon:<array, (nobs)>:  obsvervation location
    3. ext:<str>: extension for filname
    8. keywords:<dict>: extra inputs:
    
    
    """
    
    field_lst=fdesc.read_file(yyyy, mm, dd, \
                                 olon, olat,\
                                 ext=ext,\
                                 **keywords)
        
    
    
    return fdesc





def open_enr_file(datapath, flnm,  \
                      ftracerinfo,\
                      fdiaginfo,\
                      yyyy, mm,dd, \
                      categorys=None,\
                      tracers=None,\
                      taus=None,\
                      tranames=None,\
                      is_enr_file=True,\
                      **fio_keywords):
    
    
    """ setup fdesc for GEOS-Chem ensemble runs
    
    
    Inputs:
    --------------------------------------------------------------

    1. datapath:<str>: path 
    2. flnm:<str>: file name
    3, ftracerinfo:<str>: name for tracerinfo file
    4. fdiaginfo:<str>: name for diaginfo file

    5. categorys:<list, t:str>: list of categorys to be read
    6. tracers:<list,t:int>: list of tracers to be read 
    7. taus:<list, t:float>: list of taus to be read 
    8. tranames:<list, t:str>:list of tracer name
    9. is_enr_file:<T/F>: True if it is for ensemble run file 
    10. fio_keywords: <dict>: extra inputs for file IO
    
    ---Reserved keywords:
    
    --->enr_yst_lst:<list, t:int>:list of years for ensemble run 
    --->enr_dst_lst:<list, t:int> list of doy for ensemble run 
    --->enr_step_lst:<list, int>: list of step for ensemble run 
    --->enr_emst_lst:<list, int>: list of start ensemble for ensemble run
    --->enr_emend_lst:<list, int>: list of start ensemble for ensemble run
    
    Returns:
    --------------------------------------------------------
    1. fdesc:<file_desc_cl>: class for file access

    """

    # S1 read in the list 
    
    fdesc=ofo.gcfile_desc_cl(flnm, ftracerinfo, \
                            fdiaginfo,\
                            yyyy, mm, dd,\
                            is_enr_file=is_enr_file,\
                            categorys=categorys,\
                            tracers=tracers,\
                            taus=taus,\
                            tranames=tranames,\
                            mask_val=oob.fill_val,\
                            fio_keywords=fio_keywords)
            
    # S2: set file format: 
    # ##T1: filename
    
    fdesc.set_file_path(datapath)
    if (len(flnm)>0):
        fdesc.set_filename_format(flnm)
        
    else:
        fdesc.set_filename_format(def_enr_filename_format)
    
    # ##T2: diaginfo 
      
    fdesc.set_diaginfo_file_path(datapath)
            
    if (len(fdiaginfo)>0):
        fdesc.set_diaginfo_filename_format(flnm)
    else:
        fdesc.set_diaginfo_filename_format(def_enr_fdiaginfo_format)
    
    
    # ##T3: tracerinfo 
      
    fdesc.set_tracerinfo_file_path(datapath)
    
    if (len(ftracerinfo)>0):
        fdesc.set_tracerinfo_filename_format(flnm)
    else:
        fdesc.set_tracerinfo_filename_format(def_enr_ftracerinfo_format)
    
    
    
    # S4 set path 
    
    return fdesc
 

def get_mod_pres(fdesc, **keywords):
    
    """
    
    Calculate or find model pressure 
    
    Inputs:
    ----------------------------------
    1. fdesc:<gcfdesc_cl>: file descript to GEOS-Chem 
    2. keywords:<dict>: for extra inputs
    
    ---reserved keys:
    ---traname:<str>: traname for pressure or sp 
    ---use_sp:<T/F>: True, surface pressure will be used to calculate pressures
    ---geos_ver:<int>: GEOS-Chem version
    ---use_reduced:<int>: 1 means reduced model vertical grid is used  
    ---Reserved keys:
    ---category:<str>: category of the surface pressure or pressure
    ---tracer:<int>: tracer ID of surface pressure or model pressure
    ---tau:<float>: time of surface pressure or model pressure
        
    returns:
    --------------------------------------------------
    1. mod_pres:<array>: model pres
        
    """
    mod_pres=fdesc.compute_mod_pres(**keywords)
    return mod_pres



    
    
def get_mod_sp(fdesc, sp_traname='PSURF', sp_category=None, sp_tracer=None):
    
    """
    get model surface pressure
           
    """
    
    traname_lst, category_lst, \
        tid_lst, sid_lst, found_data_lst=fdesc.get_data(traname=sp_traname,\
                                                            category=sp_category, tracer=sp_tracer)
    
    if (len(found_data_lst)>0):
        sp=found_data_lst[0]
        sp=npy.array(sp)
        return sp
    else:
        return None
    

def get_data(fdesc, traname, category, tracer):
    """
    find data according to criteria of traname, category and tracer ID

    
    Inputs:
    -------------------------------------
    1. fdesc:<gcfdesc_cl>: file access
    2. traname:<str>: tracer name
    3. category:<str>: category
    4. tracer:<int>: tracer ID 
    
    Returns:
    ---------------------------------------------
    found_data_lst:<list, t:array>: selected data set 
    

    """
    
    traname_lst, category_lst, \
        tid_lst, sid_lst, found_data_lst=fdesc.get_data(traname=traname,\
                                                            category=category, tracer=tracer)
    
    return found_data_lst

def get_mod_gp(fdesc, traname='CO2', category=None, tracer=None):
    
    """
    get model profile for selected variables
    
    
    Inputs:
    -------------------------------------
    1. fdesc:<gcfdesc_cl>: file access
    2. traname:<str>: tracer name
    3. category:<str>: category
    4. tracer:<int>: tracer ID 
    
    Returns:
    ---------------------------------------------
    gp_lst:<array>: selected data set 
    

    Notes:
    --------------------------------------------------
    1. if fdesc.data_lst contains list of profile,  gp_lst will be an array in the shape 
    of [nobs, nlvl, ntracer (nem)] or shape of [nobs, ntracer]
    
    
    
        
    """
    
    traname_lst, category_lst, \
        tid_lst, sid_lst, found_data_lst=fdesc.get_data(traname=traname,\
                                                            category=category, tracer=tracer)
    # print 'traname_lst', traname_lst
    
    # print 'sid', sid_lst
    # print 'tid_lst', tid_lst
    

    
    if (len(found_data_lst)>0):
        new_idx=range(len(sid_lst))
        
        
        # new_idx=npy.argsort(sid_lst)
        # print 'new_idx', new_idx
        # ttt=raw_input()
        
        gp_lst=list()
        for xid in new_idx:
            prof=found_data_lst[xid]
            dims=npy.shape(prof)
            
            gp_lst.append(prof)
        
        gp_lst=npy.array(gp_lst)
        dims=npy.shape(gp_lst)
        
        if (len(dims)==4):
            gp_lst=npy.transpose(gp_lst, [1,2,3,0])
        
        elif (len(dims)==3):
            gp_lst=npy.transpose(gp_lst, [1,2,0])
        else:
            gp_lst=npy.transpose(gp_lst, [1,0])
        
        return gp_lst
    
    else:
        
        return None
    
        
    
def setup_enr_daily_profile(fdesc, yyyy, mm, dd, \
                                olon, olat,\
                                enr_yst_lst,\
                                enr_dst_lst,\
                                enr_step_lst,\
                                enr_em_st_lst,\
                                enr_em_end_lst, **keywords):
    
    """ set up  model profile ensemble from time-series outputs of ensemble runs

    Inputs:
    ----------------------------------------------------------------------------------
    1. fdesc:<gcfdesc_cl>: file access

    2. yyyy, mm,dd:<int>:  year, month, day
    3. olat, olon:<array, (nobs)>:  obsvervation location
    4. enr_yyyy_lst:<int>: year for ensemble run start.
    5. enr_doy_lst:<int>:  doy for ensemble run start.
    6. enr_step_lst:<int>: step for ensemble run.
    7. enr_em_st_st:<int>: firt one in the ensemble. 
    8. enr_em_end_lst:<int>: last one in the ensemble. 
    9. keywords:<dict>: extra inputs:
    
    
    """
    nset=fdesc.read_enr_file_profile_list(yyyy, mm, dd, \
                                              olon, olat,\
                                              enr_yst_lst,\
                                              enr_dst_lst,\
                                              enr_step_lst,\
                                              enr_em_st_lst,\
                                              enr_em_end_lst,\
                                              **keywords)
        
    
    
    return fdesc


def setup_enr_daily_field(fdesc, yyyy, mm, dd, \
                              olon, olat,\
                              enr_yst_lst,\
                              enr_dst_lst,\
                              enr_step_lst,\
                              enr_em_st_lst,\
                              enr_em_end_lst, **keywords):
    
    """ set up  model field ensemble from time-series outputs of ensemble runs
    
    Inputs:
    ----------------------------------------------------------------------------------
    1. fdesc:<gcfdesc_cl>: file access

    2. yyyy, mm,dd:<int>:  year, month, day
    3. olat, olon:<array, (nobs)>:  obsvervation location
    4. enr_yyyy_lst:<int>: year for ensemble run start.
    5. enr_doy_lst:<int>:  doy for ensemble run start.
    6. enr_step_lst:<int>: step for ensemble run.
    7. enr_em_st_st:<int>: firt one in the ensemble. 
    8. enr_em_end_lst:<int>: last one in the ensemble. 
    9. keywords:<dict>: extra inputs:
    
    
    """
    nset=fdesc.read_enr_file(yyyy, mm, dd, \
                                 olon, olat,\
                                 enr_yst_lst,\
                                 enr_dst_lst,\
                                 enr_step_lst,\
                                 enr_em_st_lst,\
                                 enr_em_end_lst,\
                                 **keywords)
        
    
    
    return fdesc





if  (__name__=='__main__'):
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    yyyy, mm,dd=2009, 2, 20
    
    fdesc=open_enr_file(datapath, "",  \
                                    "",\
                                    "",\
                                    yyyy, mm,dd)
    
    enr_yst_lst=[2009, 2009]
    enr_dst_lst=[1,  1]
    enr_step_lst=[1, 1]
    enr_em_st_lst=[1, 81]
    enr_em_end_lst=[81, 146]
    olon=npy.arange(-40, 40, 20)
    olat=npy.arange(-40, 40, 20)
    olat=olat+10
    
    
    fdesc=setup_enr_daily_profile(fdesc, yyyy, mm, dd, \
                                      olon, olat,\
                                      enr_yst_lst,\
                                      enr_dst_lst,\
                                      enr_step_lst,\
                                      enr_em_st_lst,\
                                      enr_em_end_lst)
    
                                                
    sp=get_mod_sp(fdesc)
    
    # print npy.shape(sp)
    
    prof_co2=get_mod_gp(fdesc)
    
    # print npy.shape(prof_co2)
    
    


    datapath='/home/lfeng/local_disk/otool_data/enkf_rerun/'
    yyyy, mm,dd=2009, 2, 20
    
    fdesc=open_ts_file(datapath, 'ts_satelliteXEXTX.XYYYYXXMMXXDDX.bpch',  \
                           "",\
                           "",\
                           yyyy, mm,dd)
    
    olon=npy.arange(-40, 40, 20)
    olat=npy.arange(-40, 40, 20)
    olat=olat+10
    
    
    fdesc=setup_daily_profile(fdesc, yyyy, mm, dd, \
                                olon, olat, ext='.EN0001-EN0002')
    
                                                
    sp=get_mod_sp(fdesc)
    
    print npy.shape(sp)
    
    prof_co2=get_mod_gp(fdesc)
    
    print npy.shape(prof_co2)
    
    
