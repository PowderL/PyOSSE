""" module for generating model observations at given locations of satellite (or dummy) observations 
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2013.01.21
   History: v0.95, 2013.04.23
    
   Notes:
   ------------------------------------------------------------------   
   1. Observation file are expected to  contain at least
   ---> lon, lat, time # location and 
    ---->obs_ak, obs_pres, obs_apr: averaging kernel, observation retrieval pressure level, and retrieval apriori.


    2. procedures:
    ---> decode menu (configuration file)
    ---> create configuration
    ---> create satellite observation operator
    ---> read in satellite dummy observation 
    ---> sample model profile and calculate model xgp
    ---> save model xgp to file
    

"""

import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_menu_m as menu_m
import ESA.observation.sobs_def_m as sobs_def_m

import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import numpy.random as rnd


# S1: read in configuration file 

menu_flnm='sob_def.cfg'

if ('.xml' in menu_flnm):
    # #: if the menu is given in xml  file 
    root_menu=menu_m.xml_to_menu(xml_menu_flnm)
else:
    # if the menu is given in format text file 
    root_menu=menu_m.txt_to_menu(menu_flnm)
    
# #c: show definition 
    
root_menu.print_menu()
        

# S2: date 

# #c: start time 
yyyy=root_menu['cfg.yst']
mm=root_menu['cfg.mst']
dd=root_menu['cfg.dst']

doy=tm.day_of_year(yyyy, mm, dd)

# #c: end time 
yend=root_menu['cfg.yend']
mend=root_menu['cfg.mend']
dend=root_menu['cfg.dend']

doy_end=tm.day_of_year(yend, mend, dend)

ndoy=0

if (yyyy==yend):
    # #c: in the same year
    
    ndoy=doy_end-doy+1
else:
    # not in the same year 
    if (npy.mod(yyyy, 4)==0):
        ndoy=366-doy+1
    else:
        ndoy=365-doy+1

    # days from year start  to year start 
    for cur_yyyy in range(yyyy+1, yend):
        if (npy.mod(cur_yyyy,4)==0):
            ndoy=ndoy+366
        else:
            ndoy=ndoy+365
    # final doy
    ndoy=ndoy+doy_end

# S3: create observation definition  class


odef=sobs_def_m.sobs_def_cl(yyyy, mm, dd, \
                                root_menu)

# S4: create observation operator 

cl_sat_opr=root_menu['sat_opr.fclass'](odef.cl_obs,\
                                           odef.cl_enr_prof,\
                                           odef.cl_fc_prof,\
                                           fxgp=root_menu['sat_opr.fxgp'],\
                                           fhm=root_menu['sat_opr.fhm'])

# S5 generate and save model observation day by day


days_in_year=365
if (npy.mod(yyyy, 4)==0):
    days_in_year=366

viewtype=root_menu['cfg.viewtype']
viewmode=root_menu['cfg.viewmode']
outpath=root_menu['cfg.outpath']
outname=root_menu['cfg.outname']

for  idd in range(ndoy):
    
    yyyy, mm, dd=tm.doy_to_time_array(doy, yyyy)
    
    # S3: read observation  
    
    cl_sat_opr.read_obs(yyyy, mm, dd, viewtype=viewtype, viewmode=viewmode)
    
    
    # S4: sample model output 
    cl_sat_opr.read_sel_xgp(yyyy, mm, dd, \
                                gpname_lst=['PSURF', 'CO2'],\
                                sel_idx=0,\
                                sel_traname='CO2',\
                                append_data=False,\
                                do_update=True, \
                                do_save=False,\
                                xgp_flnm=outpath+'/mod_xgp')
    

    
    # S5: generate random error
    
    nobs, nlvl=npy.shape(cl_sat_opr.obs_apr)
    
    # S6: save data

    # #c: filename
    
    daily_obsflnm=ovar.construct_filename(outpath, outname, \
                                              XYYYYX=yyyy,\
                                              XMMX=mm,\
                                              XDDX=dd,\
                                              XDOYX=idd)
    
    
    # #c: set dimensions
    xlvl=npy.arange(nlvl)
    xno=npy.arange(nobs)
    
    
    xno_info=ovar.io_var_cl('nobs', 'i', ['nobs'], xno, \
                                varattr=None)
    
    xlvl_info=ovar.io_var_cl('nlvl', 'i', ['nlvl'], xlvl, \
                                 varattr=None)
    
    # #c: save the data
    # #c: save model values  as both model and obs 
    
    mod_info=ovar.io_var_cl('mod', 'f', ['nobs'], cl_sat_opr.mod_xgp)
    
    obs_info=ovar.io_var_cl('obs', 'f', ['nobs'], cl_sat_opr.mod_xgp)
    
    # #c: observation uncertainty and random error
    ocnt=cl_sat_opr.ocnt
    oerr=1.0e-6*cl_sat_opr.oerr
    oerr=oerr/npy.sqrt(ocnt)
    rnd_err=rnd.normal(scale=oerr, size=nobs)
    
    oerr_info=ovar.io_var_cl('oerr', 'f', ['nobs'], oerr)
    
    rnd_err_info=ovar.io_var_cl('rnd_err', 'f', ['nobs'], rnd_err)
    
    
    olon_info=ovar.io_var_cl('lon', 'f', ['nobs'], cl_sat_opr.olon)
    olat_info=ovar.io_var_cl('lat', 'f', ['nobs'], cl_sat_opr.olat)
    sp_info=ovar.io_var_cl('sp', 'f', ['nobs'], cl_sat_opr.sp)
    cnt_info=ovar.io_var_cl('cnt', 'f', ['nobs'], ocnt)
    otime_info=ovar.io_var_cl('time', 'f', ['nobs'], cl_sat_opr.otime)
    
    ak_info=ovar.io_var_cl('obs_ak', 'f', ['nobs', 'nlvl'], cl_sat_opr.obs_ak)
    apr_info=ovar.io_var_cl('obs_apr', 'f', ['nobs', 'nlvl'], cl_sat_opr.obs_apr)
    pres_info=ovar.io_var_cl('obs_pres', 'f', ['nobs', 'nlvl'], cl_sat_opr.obs_pres)
    
    
    ncfio.ncf_save_var(daily_obsflnm, [olon_info,\
                                           olat_info,\
                                           sp_info,\
                                           cnt_info, \
                                           otime_info,\
                                           obs_info,\
                                           oerr_info,\
                                           rnd_err_info,\
                                           mod_info, ak_info, \
                                           apr_info, pres_info], \
                           [xno_info, xlvl_info], create_new=True)
    

    
    
    doy=doy+1
    if (doy>days_in_year):
        # reach
        yyyy=yyyy+1
        doy=1
        days_in_year=365
        if (npy.mod(yyyy, 4)==0):
            days_in_year=366
    

        
