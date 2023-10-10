""" calculate jacobian from model outputs for dry-air columns measured satellite instruments. 


    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.22
    History: v0.95, 2013.03.02
    

    
    Classes:
    =============================================
    1. diag_info_cl: class for diaginfo file
    2. tracer_info_cl: class for tracer info file 
    3. gcfile_desc_cl: class for accessing GEOS-Chem bpch2 files
    

    Functions:
    ================================================
    
    1.  setup_vert_intpl: Setup vertical interpolate class 
    2.  interpl_mod_prof: Interpolate model profile (or ensemble) to the observation vertical grid 
    3. calculate_xgp: calculate dry-air column using averaging kernels for model profiles represented at observation 
    
    
"""
import ESA.util.vertical_interp_m as vintpl_m
import ESA.util.vertical_profile  as vpf
import numpy as npy
import ESA.util.otool_obj as oob

def setup_vert_intpl(obs_pres, mod_pres, is_in_log, vintpl=None):
    

    """
    Setup vertical interpolate class
    
    Inputs:
    ------------------------------------------------------------
    1. obs_pres:<array, (nobs, obs_nlvl)>: observation (retrieval) pressure levels 
    2. mod_pres:<array, (nobs, mod_nlvl)>: model pressure levels 
    3. is_in_log:<T/F>: True if the model and observation pressure grid are in log.  

    Returns:
    --------------------------------------------------------------
    1. vintpl: <vintpl_cl>: class for vertical interpolation 
    
    """
    
    # S1: check if the pressure (and profile) is given in the order of top-surface   
    if (vintpl==None):
        do_reverse=False
        
        print 'shape mod_pres:', npy.shape(mod_pres)
        
        if (mod_pres[0,2]>mod_pres[0,3]):
            do_reverse=True
            
        # #T: setup the vertical interpolation  class
    
        vintpl=vintpl_m.vintpl_cl(mod_pres, is_in_log, do_reverse, mask_val=oob.fill_val)
    
   # S2: configure the class for interpolations between obs and model vertical grid 
    
    do_obs_reverse=False
    
    
    if (obs_pres[0,2]>obs_pres[0,3]):
        
        do_obs_reverse=True
        
    vintpl.init_interp(obs_pres, is_in_log, do_ob_reverse=do_obs_reverse)
    
    
    return vintpl


def interpl_mod_prof(vintpl, mod_prof):
    
    """ 
    Interpolate model profile (or ensemble) to the observation vertical grid 
    
    Inputs:
    ------------------------------------------------------------
    1. vintpl: <vintpl_cl>: class for vertical interpolation 
    2. mod_prof:<array, (nob, mod_nlvl [,nem])>: model profiles 
    
    Returns:
    --------------------------------------------------------------
    1. prof_at_ob::<array, (nob, mod_nlvl [,nem])>: model profiles at observation grid 
    
    
    
    """        
    prof_at_ob=vintpl.interpolate_mod_prof(mod_prof)
    
    return prof_at_ob

def calculate_xgp(vintpl, obs_apr, obs_ak, prof_at_ob):
    """
    calculate dry-air column using averaging kernels for model profiles represented at observation 
    vertical grid. 
   
    Inputs:
     -------------------------------------------------------------
    1. vintpl:<vintpl_cl>: class for vertical interpolation 
    2. obs_apr:<array, (nobs, obs_nlvl)>: observation prior profile 
    3. obs_ak:<array, (nobs, obs_nlvl)>: observation averaging kernels 
    4. prof_at_ob:<array, (nobs, obs_nlvl)>: model (or ensemble)  profile at observation vertical grid           
    
    Returns:
    ------------------------------------------------------------------------
    1. xgp0:<array, (nobs)>: column values for (1-obs_ak)*obs_apr
    2. xgp:<array, (nobs, [nem])>: column values for model profile (see Note 1)
    
    Notes:
    1. xgp is the column value of  mod_prof=(1-obs_ak)*obs_apr+obs_ak*prof_at_ob
    
    """
    # S1: construct  obs_ak_res=1.0-obs_apr
     
    
    usd_idx=npy.where(obs_ak<>vintpl.mask_val)
    obs_ak_res=npy.array(obs_ak)
        
    obs_ak_res[usd_idx]=1.0-obs_ak[usd_idx]

   # S2: column values for (1-obs_ak)*obs_apr
    
    xgp0=vintpl.get_ak_ob_column(obs_apr, obs_ak_res)
    # S3: column values for obs_ak*prof_at_ob

    
    xgp=vintpl.get_ak_ob_column(prof_at_ob, obs_ak)
    
    # S4: xgp=xgp0+xgp 
    
    gp0_dims=npy.shape(xgp0)
    gp0_ndim=npy.size(gp0_dims)
    
    gp_dims=npy.shape(xgp)
    gp_ndim=npy.size(gp_dims)
    
    if (gp_ndim>gp0_ndim):
        # if it is for ensemble 
        # print 'prof_at_obs:', prof_at_ob[0, :, 20]
        # ttt=raw_input()
        
        xgp=xgp+xgp0[:,npy.newaxis]
        
    else:
        
        # single tracer 
        xgp=xgp+xgp0
    
    return xgp0, xgp

def reset_prior_prof(obs_pres,  obs_apr, sp, filling_value=oob.fill_val):
    
    """
    
    reset prior profile with pressure level be cut to surface pressure
    
    
    """
    
    new_pres, new_prof=vpf.prof_replace_1d(obs_pres, obs_apr, sp, mask_val=filling_value)
    return new_pres, new_prof

def get_xgp(obs_pres, obs_apr, obs_ak,  mod_pres, mod_prof, is_in_log=False):
    
    """
    calculate xgp for model profiles
    Inputs: 
    ---------------------------------------------------------------------------------
    1. obs_pres:<array, (nobs, obs_nlvl)>: observation (retrieval) pressure levels 
    2. obs_apr:<array, (nobs, obs_nlvl)>: observation prior profile 
    3. obs_ak:<array, (nobs, obs_nlvl)>: observation averaging kernels 
    4. mod_pres:<array, (nobs, mod_nlvl)>: model pressure levels 
    5. mod_prof:<array, (nobs, mod_nlvl, [nem])>: model profiles 
    6. is_in_log:<T/F>: True if the model and observation pressure grid are in log.  

    Returns:
    -------------------------------------------------------------------------------      
    1. xgp0:<array, (nobs)>: column values for (1-obs_ak)*obs_apr
    2. xgp:<array, (nobs, [nem])>: column values for model profile  (see Note 1)
    
    
    Notes:
    1. xgp is the column value of  mod_prof=(1-obs_ak)*obs_apr+obs_ak*prof_at_ob
              
    """
    # S1: setup the interpolation 
    
    vintpl=setup_vert_intpl(obs_pres, mod_pres, is_in_log)

    # S2: project profiles to observation space
    dims=npy.shape(mod_prof)
    
    # if (len(dims)==2):
    #    print 'mod_prof', mod_prof[0, :]
    # else:
    #    print 'mod_prof', mod_prof[0, :, 20]
    
    prof_at_ob=interpl_mod_prof(vintpl, mod_prof)

    # S3: calculate xgp
    
    xgp0, xgp=calculate_xgp(vintpl,  obs_apr, obs_ak, prof_at_ob)
    
    return xgp0, xgp



if (__name__=='__main__'):
    import os
    import satellite_operator as sop
    em_st_list=[1]
    em_end_list=[2]
    obs_datapath='./gosat_v28_obs/'
    viewtype='gosat_v28'
    viewmod='nadir'
    data_path='./enkf_rerun/'
    yyyy=2009
    doy_st=192
    doy_end=doy_st+8
    op=sop.satellite_xgp(viewtype, obs_datapath, 0, True, True)
    for doy in range(doy_st, doy_end):
        op.read_obs(yyyy, doy)
        yyyy, mm, dd=tmdl.doy_to_time_array(doy, yyyy)
        xgp=op.read_sel_xgp(yyyy, mm, dd,\
                            1, 2,\
                            1,\
                            data_path)
        
        
         
        
               
    
           

