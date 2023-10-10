""" calculate jacobian from model outputs for dry-air columns measured satellite instruments. 

    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.22
    History: v0.95, 2013.03.02
    
    Classes:
    ---------------------------------------------------
    1. satellite_xgp_cl: observation operator
    
"""


import os

import satellite_obs as sat_obs
import ESA.util.time_module as tm
import xgp_jacobian_m as xgp_m
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ts_slice_m as ts_io_m
import ESA.atmosphere.enr_slice_m as enr_io_m
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import numpy as npy


import ESA.util.gc_ts_file_m as tsm

class satellite_xgp_cl:
    """ observation operator for 
    1) read in satellite observations 
    2) sample model outputs along location of observation 
   
    Members:
    ========================================================
    1. attr_dict:<dict>: dictionary 
    2. enr_prof_cl:<enr_prof_cl>: class for reading ensmeble profile
    3. fc_prof_cl:<fc_prof_cl>: class for reading single profile
    
    4. obs_cl:<obs_cl>: class for reading observationd
    
    # storage for observation
        
    5. olon:<array, (nobs)>: observation longitude
    6. olat:<array, (nobs)>: observation latitude 
    7. obs:<array, (nobs)>:  observation 
    8. oerr:<array, (nobs)>:  observation error
    9. obs_ak:<array, (nobs, nlvl)>:  observation averaging kernel
    10.obs_apr:<array, (nobs, nlvl)>: retrieval prior 
    11.obs_pres:<array, (nobs, nlvl)>: retrieval pressure
    12.self.oflag:<array, (nobs)>: retrieval flag
    13. ocnt:<array, (nobs)>: actual observation number before gridding etc
    14. sp:<array, (nobs)>: surface pressure
   
    
        
    
    15. mod_sp:<array, (nobs)>: model surface pressure
    16. mod_pres:<array, (nobs, nmod_lvl)>: model pressure
    17. mod_prof:<array, ((nobs, nmod_lvl)>: model pressure
    
    
    15. mod_xgp:<array, (nobs)>: model observation
   
    16. hm:<array, (nobs, ne)>: response function to perturbation
    17. nem:<int>: size of ensemble 
    18. nobs:<int>: number of observation
    19. fxgp:<func>: function for convolving profiles using averaging kernel
    20. fhm:<func>: function for convolving ensemble profiles using averaging kernel 
    21.fsubtract:<func>: function for calculating response functions by subtracting the reference values
    

    Functions
    ======================================================
    1. __init__:initialization
    2. subset_obs:subset observations
    3. read_obs: Read observations into var list
    4. read_sel_xgp: read outputs from single (forecast) run, and calculate xgp
    5. read_jacobian:read outputs from ensemble run, and calculate response functions
                
        
    """
    
    def __init__(self, \
                     obs_cl,\
                     enr_prof_cl,\
                     fc_prof_cl,\
                     fxgp=xgp_m.get_xgp,\
                     fhm=xgp_m.get_xgp,\
                     fsubtract=None,\
                     **keywords):
        
        

        """
        Inputs:
        =========================================
        1. obs_cl:<satellite_obs_cl>: class to read satellite obs 
        2. enr_prof_cl:<enr_slice_cl>: class to read model profiles for ensemble run
        3. fc_prof_cl:<ts_slice_cl>: class to read forecast (single tracer) runs.
        4. fxgp:<func>: function for calculating xgp from model profile
        5. fhm:<func>: function for calculating hm from ensemble model profile 
        
        """
        # S1: Set attribute 
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_obs_op})
        
        
        # S2: assign initial values
       
        # #T1: class for data/output access
        self.enr_prof_cl=enr_prof_cl
        self.fc_prof_cl=fc_prof_cl
        self.obs_cl=obs_cl
        # #T2: storage for observation
        
        self.olon, self.olat, self.obs, self.oerr, \
            self.obs_ak, self.obs_apr, self.obs_pres, self.oflag=\
            None, None, None, None,\
            None, None, None, None
        # #c: true obs count
        
        
        self.ocnt=None
        self.sp=None
        
        
        # #T3: storage for model values

        self.mod_xgp=None
        self.mod_pres=None
        self.mod_prof=None
        
        self.mod_sp=None
        
        self.hm=None
        
        self.nem=0
        self.nobs=0
        self.fxgp=fxgp
        self.fhm=fhm
        self.fsubtract=fsubtract
        
        
        
        
    def subset_obs(self, usd_index_lst):
        
        """
        subset the observations
        Inputs:
        -----------------------------------------------------
        1. usd_index_lst:<list, t:int>: list of observation index to be used
        

        Returns:
        ------------------------------------------------------
        1.self.nobs:<int>: observation number left. 
        
        """
        self.olon, self.olat, self.obs, self.oerr, \
            sel.obs_ak, self.obs_apr, self.obs_pres, self.oflag=\
            self.olon[usd_index_lst], self.olat[usd_index_lst], self.obs[usd_index_lst], self.oerr[usd_index_lst], \
            sel.obs_ak[usd_index_lst], self.obs_apr[usd_index_lst], self.obs_pres[usd_index_lst], self.oflag[usd_index_lst]
        self.ocnt=self.ocnt[usd_index_lst]
        self.sp=self.sp[usd_index_lst]
        
        self.nobs=npy.size(self.olon)
        
    
    def read_obs(self, yyyy, mm, dd, **keywords):
        
        """
        Read observations into var list
        
        
        
        Inputs:
        ============================================================
        1. yyyy, mm, dd:<int>: date 
        2. keywords:<dict>: dictionary for extra input 
        

~       Returns:
        ============================================================
        1. self.nobs:<int>: number of observations found 
    

        Notes:
        ===========================================================
        """
        
        self.nobs=self.obs_cl.read_obs(yyyy, mm, dd, **keywords)
        self.olon=self.obs_cl['lon']
        self.olat=self.obs_cl['lat']
        self.obs=self.obs_cl['obs']
        self.oerr=self.obs_cl['err']
        self.obs_ak=self.obs_cl['obs_ak']
        self.obs_apr=self.obs_cl['obs_apr']
        self.obs_pres=self.obs_cl['obs_pres']
        self.oflag=self.obs_cl['flag']
        self.otime=self.obs_cl['time']
        
        self.ocnt=self.obs_cl['cnt']
        
        
        if (self.ocnt==None):
            self.ocnt=npy.ones(npy.size(self.olon), float)
        
            
        self.sp=self.obs_cl['sp']
        
        if (self.sp==None):
            self.sp=npy.zeros(npy.size(self.olon), float)
            self.sp[:]=oob.fill_val # not in use
            
            

        print 'read obs: shape(oerr)', npy.shape(self.oerr)
        
        
        
        return self.nobs
    
    def read_sel_xgp(self, yyyy, mm, dd, \
                         gpname_lst=['PSURF', 'CO2'],\
                         sel_idx=0,\
                         sel_traname='CO2',\
                         append_data=False,\
                         do_update=True, \
                         do_save=True,\
                         xgp_flnm='./mod_xgp',\
                         **keywords):
        
        
        """ read outputs from single (forecast) run, and calculate xgp
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. gpname_lst:<lst:str>: name of tracers to be read in 
        3. sel_idx:<int>: index of the selected tracer
        4. sel_traname:<str>: name of the selected tracer
        5. append_data:<T/F>: if True, read data will be appended to list end
        6. do_update:<T/F>: if true, mod_xgp will be recalculated. If false, it will be read from disk file        
        7. do_save:<T/F>: if true, mod_xgp will be saved to disk file
        8. xgp_flnm:<str>: file used to store model xgp data
        9. keywords:<dict>: extra inputs:
        
        ---reserved keywords:
        --->1. categorys: <list, t:str>: category for variables to be read 
        --->2. tracers:<list, t:int>: list for tracer ID 
        --->3. taus: <list, t:float>: list for starting time 
        
        
        
"""
        
        # S1: file name for xgp data
        
        sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        daily_obsflnm=xgp_flnm+'.'+sdate+'.nc'
        
        # S2: get model xgp
        
        
        if (do_update):
            # ##T1: read outputs forecast (single tracer runs)
            
            if (gpname_lst<>None):
                if (not sel_traname in gpname_lst):
                    gpname_lst=gpname_lst+[sel_traname]
            
            self.fc_prof_cl.load_mod_output(yyyy, mm, dd, \
                                                olon=self.olon, \
                                                olat=self.olat,\
                                                gpname_lst=gpname_lst,\
                                                append_data=append_data,\
                                                **keywords)
            
            
            # ##2: get model profile
            
            self.mod_pres=self.fc_prof_cl.set_mod_pres(**keywords)
            
            self.mod_sp=self.fc_prof_cl.set_mod_pres(use_sp=0, **keywords)
            self.mod_sp=npy.squeeze(self.mod_sp)
            
            
            print 'get_xgp: shape of mod_prfes, mod_sp', npy.shape(self.mod_pres), npy.shape(self.mod_sp)
            
            # replace sp with model sp for bad values
            
            idx=npy.where(self.sp<=0.0)
            
            if (npy.size(idx)>0):
                # this case is only for gen observations 
                
                self.sp=npy.where(self.sp<=0, self.mod_sp, self.sp)
                
                # adjust obs_pres, obs_ak, and obs_apr 
                
                new_pres, self.obs_apr=xgp_m.reset_prior_prof(self.obs_pres, self.obs_apr,  self.sp)
                new_pres, self.obs_ak=xgp_m.reset_prior_prof(self.obs_pres, self.obs_ak,  self.sp)
                self.obs_pres=new_pres
                
            
            
            self.mod_prof=self.fc_prof_cl.get_sel_profile(sel_traname, sel_idx, self.otime, **keywords)
            
            
            # ##3: calcualte xgp
            
            xgp0, mod_xgp=self.fxgp(self.obs_pres, self.obs_apr, self.obs_ak, self.mod_pres, self.mod_prof, self.fc_prof_cl.is_in_log)
            self.mod_xgp=mod_xgp
            
            # ##4: save data if required
            
            if (do_save):
                # ##T1: data to be save 
                nobs, nlvl=npy.shape(self.obs_apr)
                xlvl=npy.arange(nlvl)
                xno=npy.arange(nobs)
                xno_info=ovar.io_var_cl('nobs', 'i', ['nobs'], xno, \
                                            varattr=None)
                
                xlvl_info=ovar.io_var_cl('nlvl', 'i', ['nlvl'], xlvl, \
                                             varattr=None)
                
                
                
                mod_info=ovar.io_var_cl('mod', 'f', ['nobs'], self.mod_xgp)
                obs_info=ovar.io_var_cl('obs', 'f', ['nobs'], self.obs)
                oerr_info=ovar.io_var_cl('oerr', 'f', ['nobs'], self.oerr)
                olon_info=ovar.io_var_cl('lon', 'f', ['nobs'], self.olon)
                olat_info=ovar.io_var_cl('lat', 'f', ['nobs'], self.olat)
                otime_info=ovar.io_var_cl('time', 'f', ['nobs'], self.otime)
                sp_info=ovar.io_var_cl('sp', 'f', ['nobs'], self.sp)
                mod_sp_info=ovar.io_var_cl('mod_sp', 'f', ['nobs'], self.mod_sp)
                
                ak_info=ovar.io_var_cl('obs_ak', 'f', ['nobs', 'nlvl'], self.obs_ak)
                apr_info=ovar.io_var_cl('obs_apr', 'f', ['nobs', 'nlvl'], self.obs_apr)
                pres_info=ovar.io_var_cl('obs_pres', 'f', ['nobs', 'nlvl'], self.obs_pres)
                print 'sat_opr_cl.sel_xgp--flnm :', daily_obsflnm
                
                ncfio.ncf_save_var(daily_obsflnm, [olon_info,\
                                                       olat_info,\
                                                       otime_info,\
                                                       obs_info,\
                                                       oerr_info,\
                                                       sp_info,\
                                                       mod_sp_info,\
                                                       mod_info, ak_info, apr_info, pres_info], [xno_info, xlvl_info], create_new=True)
                
                
        else:
            
            varnames=['mod']
            mod_xgp=ncfio.ncf_read(daily_obsflnm, varnames)
            mod_xgp=npy.squeeze(mod_xgp)
        
        self.mod_xgp=mod_xgp
    
        return mod_xgp
    
    def read_jacobian(self, yyyy, mm, dd, \
                          enr_yst_lst,\
                          enr_dst_lst,\
                          enr_step_lst,\
                          enr_em_st_lst,\
                          enr_em_end_lst,\
                          gpname='CO2',\
                          gpname_lst=['PSURF', 'CO2'],\
                          append_data=False,\
                          do_update=True, \
                          do_save=True,\
                          hm_flnm='obs_hm',\
                          **keywords):
        
        """                                        
        
        read outputs from ensemble run, and calculate response functions
        
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. enr_yyyy_lst:<list, t:int>: year for ensemble run start.
        3. enr_doy_lst:<list, t:int>:  doy for ensemble run start.
        4. enr_step_lst:<list, t:int>: step for ensemble run.
        5. enr_em_st_st:<list, t:int>: firt one in the ensemble. 
        6. enr_em_end_lst:<list, t:int>: last one in the ensemble.
        7. gpname_lst:<list, t:str>: name of tracers to be read in 
        8. append_data:<T/F>: if True, read data will be appended to list end
        9. do_update:<T/F>: if true, hm will be recalculated. If false, it will be read from disk file
        
        10. do_save:<T/F>: if true, hm will saved to disk file
        
        11. hm_flnm:<str>: file used to store hm data
        
        12. keywords:<dict>: extra inputs:
        ---reserved keywords:
        --->1. categorys: <list, t:str>: category for variables to be read 
        --->2. tracers:<list, t:int>: list for tracer ID 
        --->3. taus: <list, t:float>: list for starting time 
        
        
        
        Returns:
        ==========================================================
        1, self.hm:<array, (nobs, nem)>: xgp for ensemble runs
        

        Notes:
        ----------------------------------------
        1. nem may be different from those used in inversion
        """
        



        # S1: file name for hm data
        
        out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        daily_obsflnm=hm_flnm+"."+out_sdate+".nc"
        # S2: get hm 
        
        if (do_update):
            # ##1: if need, read model profiles to class self.enr_gp_cl 
            
            self.enr_prof_cl.load_mod_output(yyyy, mm, dd, \
                                                 enr_yst_lst,\
                                                 enr_dst_lst,\
                                                 enr_step_lst,\
                                                 enr_em_st_lst,\
                                                 enr_em_end_lst,\
                                                 olon=self.olon, \
                                                 olat=self.olat,\
                                                 gpname_lst=gpname_lst,\
                                                 append_data=append_data,\
                                                 **keywords)
            
            # ##2: get the ensemble model profile 
            mod_prof=self.enr_prof_cl.get_sel_profile(gpname, None, self.otime, **keywords)
            # ##3: get model pressure
            
            mod_pres=self.enr_prof_cl.set_mod_pres()
            
            if (mod_pres==None):
                mod_pres=self.mod_pres
            else:
                self.mod_pres=mod_pres
            
            # ##4: caculcate hm
            
            xgp0, hm=self.fhm(self.obs_pres, self.obs_apr, self.obs_ak, self.mod_pres, mod_prof, self.enr_prof_cl.is_in_log)
            

            if (self.fsubtract<>None):
                hm=self.fsubtract(enr_yst_lst,\
                                      enr_dst_lst,\
                                      enr_step_lst,\
                                      enr_em_st_lst,\
                                      enr_em_end_lst, hm)
                
            
            # ##5: if required save data
            
            if (do_save):
                # ##T1: dimension
                sav_nobs, sav_ne=npy.shape(hm)
                xno, xne=npy.arange(sav_nobs), npy.arange(sav_ne)
                
                xno_info=ovar.io_var_cl('nobs', 'i', ['nobs'], xno, \
                                            varattr=None)
                xne_info=ovar.io_var_cl('ne', 'i', ['ne'], xne, \
                                            varattr=None)
                
                
                # ##T2: variables to be saved
                

                
                hm_info=ovar.io_var_cl('hm', 'f', ['nobs', 'ne'], hm)
                olon_info=ovar.io_var_cl('lon', 'f', ['nobs'], self.olon)
                olat_info=ovar.io_var_cl('lat', 'f', ['nobs'], self.olat)
                otime_info=ovar.io_var_cl('time', 'f', ['nobs'], self.otime)
                print 'sat_opr_cl.read_jacobian--flnm:', daily_obsflnm
                
                ncfio.ncf_save_var(daily_obsflnm, [olon_info,\
                                                       olat_info,\
                                                       otime_info,\
                                                       hm_info], [xno_info, xne_info], create_new=True)
                
        else:
            varnames=['hm']
            hm=ncfio.ncf_read(daily_obsflnm, varnames)
            hm=npy.squeeze(hm)
        
        self.hm=hm
        nobs, self.nem=npy.shape(self.hm)
        
        return self.hm

        
         
            
if (__name__=='__main__'):

    # observation 
    
    
    obs_path='/home/lfeng/local_disk/otool_data/obs/gosat_v210_obs/'
    viewtype='gosat_v29'
    
    
    flnm='XVIEWTYPEX.XYYYYXDXDOYX.nc'
    
    yyyy, mm, dd=2009, 1, 20
    
    sat_id=3
    obs_cl=sat_obs.satellite_obs_cl(obs_path, \
                                       flnm, \
                                       viewtype, \
                                       sat_id, \
                                       varname_lst=sat_obs.sat_obs_varname_lst,\
                                       varname_dict=sat_obs.sat_obs_varname_dict, \
                                       mask_val=oob.fill_val)
    
    
    


            
                
    
    datapath='/home/lfeng/local_disk/otool_data/enkf_rerun/'
    
    
    filename, fdiaginfo, ftracerinfo="", "", ""
    name='test'
    
    fc_prof_cl=ts_io_m.gc_slice_cl(datapath, filename, \
                                       fdiaginfo,\
                                       ftracerinfo,\
                                       name, \
                                       fopen=tsm.open_ts_file,\
                                       fread=tsm.setup_daily_profile,\
                                       fget=tsm.get_mod_gp,\
                                       fpres=tsm.get_mod_pres,\
                                       ext='.EN0001-EN0002'\
                                       )
    
    
    
    
    

    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    
    filename, fdiaginfo, ftracerinfo="", "", ""
    name='test'
    
    enr_prof_cl=enr_io_m.enr_slice_cl(datapath, filename, \
                                          fdiaginfo,\
                                          ftracerinfo,\
                                          name, \
                                          fopen=tsm.open_enr_file,\
                                          fread=tsm.setup_enr_daily_profile,\
                                          fget=tsm.get_mod_gp,\
                                          fpres=tsm.get_mod_pres\
                                          )
    
    # create the operator 


    yyyy, mm,dd=2009, 2, 20
    
    sat_operator=satellite_xgp_cl(obs_cl,\
                                      enr_prof_cl,\
                                      fc_prof_cl,\
                                      fxgp=xgp_m.get_xgp,\
                                      fhm=xgp_m.get_xgp)

    
    print 'yyyy, mm, dd', yyyy, mm, dd
    
    viewmode=''
    sat_operator.read_obs(viewmode, yyyy, mm, dd)

    print 'read fc XCO2'

    sat_operator.read_sel_xgp(yyyy, mm, dd, \
                                  gpname_lst=['PSURF', 'CO2'],\
                                  sel_idx=0,\
                                  sel_traname='CO2',\
                                  append_data=False,\
                                  do_update=True, \
                                  do_save=True,\
                                  xgp_flnm='./mod_xgp')

        
    
    enr_yst_lst=[2009, 2009]
    enr_dst_lst=[1,  1]
    enr_step_lst=[1, 1]
    enr_em_st_lst=[1, 81]
    enr_em_end_lst=[81, 146]

    print 'read ensemble  XCO2'
    
    sat_operator.read_jacobian(yyyy, mm, dd, \
                                   enr_yst_lst,\
                                   enr_dst_lst,\
                                   enr_step_lst,\
                                   enr_em_st_lst,\
                                   enr_em_end_lst,\
                                   gpname_lst=['PSURF', 'CO2'],\
                                   append_data=False,\
                                   do_update=True, \
                                   do_save=True,\
                                   hm_flnm='obs_hm')
    
    
