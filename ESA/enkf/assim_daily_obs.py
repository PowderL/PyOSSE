"""
   this code is used to construct state vector 
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.15
   
   Functions:
   
   =================================================
   1. assim_satellite_obs:    
   
   
   


"""
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_descfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio
import ESA.enkf.run_desc_file_m as rdescio
import ESA.enkf.run_desc_m as rdesc
import ESA.enkf.state_vector as stv_m
import ESA.enkf.x2flux_m as x2flux
import ESA.enkf.construct_state_vector as gen_stv
import ESA.util.otool_menu_m as menu_m
import ESA.enkf.assim_def_m as defm


def assim_read_config(menu_flnm):
    """
    read assimilation configuration 
    
    Input:
    ---------------------------------
    1. menu_flnm: <str>: file name for menu 
    

    Returns:
    --------------------------------
    1. root_menu:<menu_cl>: menu for configurations 
    
    """
    root_menu=menu_m.txt_to_menu(menu_flnm)
    root_menu.print_menu()
    return root_menu 

def assim_create_objects(yyyy, mm, dd, root_menu):
    """
    Create class for assimilation 
    
    Inputs:
    -------------------------------------------
    1. yyyy, mm, dd:<int>: date
    2. root_menu:<menu_cl>: menu for configurations 

    Returns:
    -----------------------------------------
    1. cl_def:< assim_def_cl>: container for object 
    
    """
    cl_def=defm.assim_def_cl(yyyy, mm, dd, root_menu)

    return cl_def


    

def assim_digest_satellite_obs(cfg, yyyy, mm, dd, \
                                cl_stv, \
                                cl_sat_opr,\
                                cl_run_desc, \
                                hm_do_update=True,\
                                hm_save=True,\
                                mean_do_update=True,\
                                mean_save=True,\
                                **keywords):
    """
    
    Adding apriori for one step into the state vector
    
    Inputs:
    -------------------------------------------------

    1. yyyy, mm, dd:<int>: date
    2. cl_run_desc:<run_desc_cl>: class for ensemble run description 
    3. cl_stv: <state_vector>: class for state vector 
    4. cl_x2flux:<x2flux_cl>:  class for project x to fluxes 
    5. tag_fmt:<str>:  format for tag used to identify format 
    6. ftag   :<func>: function for generating tag 
    7. fscale_m:<func>: for generating correlation 
    8. keywords:<dict>: extra inputs
    

    Returns:
    --------------------------------------------------
    1. iadd_x:<int>: number of state vector added 
    
    """

    # S1: adjust state vector if necessary 
        
    
    
    
    # S2: read in observations
    
    nobs=cl_sat_opr.read_obs(yyyy, mm, dd, viwemode=cfg['sat_obs.viewmode'], \
                                 viewtype=cfg['sat_obs.viewtype'],\
                                 **keywords)
    
    if (nobs==0):
        return 0, None, None, None, None
    

    # S3: read mean upate 
    
    mod_xgp=cl_sat_opr.read_sel_xgp(yyyy, mm, dd, \
                                        gpname_lst=cfg['sat_opr.gpname_lst'],\
                                        sel_idx=cfg['sat_opr.sel_idx'],\
                                        sel_traname=cfg['sat_opr.traname'],\
                                        do_update=mean_do_update, \
                                        do_save=mean_save,\
                                        xgp_flnm=cfg['sat_opr.xgp_flnm'])
    
    # S4: read in hm 
    
    enr_yst_lst=cl_run_desc.st_year_lst
    enr_dst_lst=cl_run_desc.st_doy_lst
    enr_step_lst=cl_run_desc.step_lst
    enr_em_st_lst=cl_run_desc.st_no_lst
    enr_em_end_lst=cl_run_desc.end_no_lst

        
    hm=cl_sat_opr.read_jacobian(yyyy, mm, dd, \
                                    enr_yst_lst,\
                                    enr_dst_lst,\
                                    enr_step_lst,\
                                    enr_em_st_lst,\
                                    enr_em_end_lst,\
                                    gpname=cfg['sat_opr.traname'],\
                                    gpname_lst=cfg['sat_opr.gpname_lst'],\
                                    append_data=False,\
                                    do_update=hm_do_update, \
                                    do_save=True,\
                                    hm_flnm=cfg['sat_opr.hm_flnm'], **keywords)
    
    
    
    # S5: shift y and dy as a result of previous assimilations. 
    
    
    
    mean_y, dy=cl_stv.convolve_y_dy(mod_xgp, hm,  use_wnd_xtm=cfg['stv.use_wnd_xtm'], \
                                   use_wnd_inc_m=cfg['stv.use_wnd_inc_m'])
    
    

    # assimilation 

    yobs=cl_sat_opr.obs
    yerr=cl_sat_opr.oerr
    
    yerr=yerr*yerr
    
    xinc, xtm_inc_m, xtm, inc_m=cl_stv.do_assim(dy, mean_y, yobs, yerr)       
    
    return nobs, xinc, xtm_inc_m, xtm, inc_m





if (__name__=='__main__'):
    menu_flnm='./assim_def.cfg'
    root_menu=assim_read_config(menu_flnm)
    yyyy, mm, dd=2009, 1,1
    cl_def=assim_create_objects(yyyy, mm, dd, root_menu)
    
    cl_run_desc=cl_def.cl_run_desc
    cl_stv=cl_def.cl_stv
    cl_x2flux=cl_def.cl_x2flux
    fprior=root_menu['stv.fprior']
    cl_sat_opr=cl_def.cl_sat_opr
    
    
    
    fprior(yyyy, mm, dd, \
               cl_run_desc, cl_stv, \
               cl_x2flux=cl_x2flux,\
               tag_fmt=root_menu['stv.tag_fmt'],\
               ftag=root_menu['stv.ftag'], \
               fscale_m=root_menu['stv.fscale_m'])
    
    
    
    assim_digest_satellite_obs(root_menu, yyyy, mm, dd, \
                                   cl_stv, \
                                   cl_sat_opr,\
                                   cl_run_desc, \
                                   hm_do_update=True,\
                                   hm_save=True,\
                                   mean_do_update=True,\
                                   mean_save=True)
    
