"""
 Asistent functions for generating ensemble fluxes 
 Authors: L. Feng, Edinburgh University
 History: v0.9, 2012.06.28
 History: v0.95, 2013.02.21
 
    
 Functions:
 ===============================================
 1. define_ens_step: time step for ensemble fluxes
    

"""

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.otool_var_io as ovio
import ESA.util.geo_constant as gc
import ESA.util.time_module as tm

def define_ens_step(yyyy_st, mm_st, dd_st, time_step, nstep, **keywords):
    
    """
    define  time step for ensemble fluxes 
    
    Inputs:
    --------------------------------------------------
    
    1. yyyy_st:<int>: start year 
    2. mm_st:<int>: start month 
    3. dd_st:<int>: start day 
    4. time_step:<int>: time interval for fluxes 
    5. nstep:<int>: numbe of step
    6. keywords:<dict>: dictionary for extra inputs
    ---Reserved entries:
    -->reset_at_new_year: <int>: 1=move doy to doy=1 for each new year 
    -->reset_at_new_month: <int>: 1=move doy to the fist day of each new month 
    
    
    Outputs:
    --------------------------------------------------
    1.ens_yyyy_lst:<list, t:int>: starting year
    2.ens_doy_lst:<list, t:int>: starting doy 
    
    """

    reset_at_new_year=1
    reset_at_new_month=1

    if ('reset_at_new_year' in keywords):
        reset_at_new_year=keywords['reset_at_new_year']
    
    
    if ('reset_at_new_month' in keywords):
        reset_at_new_month=keywords['reset_at_new_month']
    
    
    
    ens_doy_lst=list()
    ens_yyyy_lst=list()
    
    # S1: date for first step 

    
    cur_dd=dd_st
    cur_mm=mm_st
    cur_yyyy=yyyy_st
    
    cur_doy=tm.day_of_year(cur_yyyy, cur_mm, cur_dd)
    ens_doy_lst.append(cur_doy)
    ens_yyyy_lst.append(cur_yyyy)
    ndays_year=tm.days_in_year(cur_yyyy)
    
    # S2: loop over nstep , and end with nstep+1
    
    
    for istep  in range(nstep+1):
        
        new_doy=cur_doy+time_step
        # #c: check year begin
        
        if (new_doy>ndays_year):
            cur_yyyy=cur_yyyy+1
            if (reset_at_new_year):
                new_doy=1
            else:
                new_doy=new_doy-ndays_year
            
            ndays_year=tm.days_in_year(cur_yyyy)

        new_dd, new_mm, new_yyyy=tm.doy_to_time_array(new_doy, cur_yyyy)
        
        # #c: check month begin 
        
        if (new_mm<>cur_mm):
            cur_mm=new_mm
            if (reset_at_new_month):
                new_doy=tm.day_of_year(cur_yyyy, cur_mm, 1)
        
        ens_doy_lst.append(new_doy)
        ens_yyyy_lst.append(cur_yyyy)
    
        cur_doy=new_doy
        
        
    # loop istep end
    
    return ens_doy_lst, ens_yyyy_lst
    
                                   


