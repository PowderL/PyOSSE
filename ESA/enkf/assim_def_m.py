""" calculate jacobian from model outputs for dry-air columns measured satellite instruments. 

    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.22
    History: v0.95, 2013.03.02
    
    Classes:
    ---------------------------------------------------
    1. satellite_xgp_cl: observation operator
    
"""


import os
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ts_slice_m as ts_io_m
import ESA.atmosphere.enr_slice_m as enr_io_m
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import numpy as npy
import ESA.util.gc_ts_file_m as tsm

import ESA.enkf.run_desc_file_m as rdescio
import ESA.enkf.run_desc_m as rdesc
import ESA.enkf.state_vector as stv_m
import ESA.enkf.x2flux_m as x2flux



class assim_def_cl:
    """
    
    class for configuring assimilation process 
    
    Members:
    =============================================================
    1. cl_obs: class for reading observation data
    2. cl_fc_prof: class for sampling single tracer model run
    3. cl_enr_prof: class for sampling 
    4. cl_sat_opr: class for observation operator (sampling models)
    5. cl_x2flux: class for projecting coefficients to fluxes
    6. cl_run_desc: class for description outputs from ensemble runs 
    

    Functions:
    1. __init__: initialization 
    
    
    """
    def __init__(self, \
                     yyyy, mm, dd, \
                     root_menu, **keywords):
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_def})
        
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        # S1: load menu 
        self.root_menu=root_menu
        
        
        # S2: create observation class
                
        fc=self.root_menu['sat_obs.fclass']
        self.cl_obs=fc(self.root_menu['sat_obs.path'], \
                           self.root_menu['sat_obs.flnm'], \
                           self.root_menu['sat_obs.viewtype'], \
                           self.root_menu['sat_obs.sat_id'], \
                           varname_dict=self.root_menu['sat_obs.varname_dict'],\
                           varname_lst=self.root_menu['sat_obs.varname_lst'],\
                           mask_val=oob.fill_val)
        
        # S3: class for single profile  
        
        self.cl_fc_prof=self.root_menu['fc_prof.fclass'](self.root_menu['fc_prof.path'],\
                                                             self.root_menu['fc_prof.flnm'], \
                                                             self.root_menu['fc_prof.ftracerinfo'],\
                                                             self.root_menu['fc_prof.fdiaginfo'],\
                                                             self.root_menu['fc_prof.vname'], \
                                                             fopen=self.root_menu['fc_prof.fopen'],\
                                                             fread=self.root_menu['fc_prof.fread'],\
                                                             fget=self.root_menu['fc_prof.fget'],\
                                                             fpres=self.root_menu['fc_prof.fpres'],\
                                                             ext=self.root_menu['fc_prof.ext']\
                                                             )
                                                     
        
    
        # S4: class for ensemble profiles 
        
        self.cl_enr_prof=self.root_menu['enr_prof.fclass'](self.root_menu['enr_prof.path'],\
                                                               self.root_menu['enr_prof.flnm'], \
                                                               self.root_menu['enr_prof.ftracerinfo'],\
                                                               self.root_menu['enr_prof.fdiaginfo'],\
                                                               self.root_menu['enr_prof.vname'], \
                                                               fopen=self.root_menu['enr_prof.fopen'],\
                                                               fread=self.root_menu['enr_prof.fread'],\
                                                               fget=self.root_menu['enr_prof.fget'],\
                                                               fpres=self.root_menu['enr_prof.fpres']\
                                                               )
        
        
        # S5: observation operator
        
        
        self.cl_sat_opr=root_menu['sat_opr.fclass'](self.cl_obs,\
                                                        self.cl_enr_prof,\
                                                        self.cl_fc_prof,\
                                                        fxgp=root_menu['sat_opr.fxgp'],\
                                                        fhm=root_menu['sat_opr.fhm'], \
                                                        fsubtract=self.root_menu['enr_prof.fsubtract']\
                                                        )
        
        
        
        # S5: create class for x to flux transform 
        
        
        
        # create class
        
        datapath=root_menu['x2flux.path']
        flnm=root_menu['x2flux.flnm']
        
        self.cl_x2flux=root_menu['x2flux.fclass'](datapath, flnm, \
                                                     yyyy, mm, dd, \
                                                     bf_get=root_menu['x2flux.bf_get'], \
                                                     bf_cl=root_menu['x2flux.bf_cl'], \
                                                     lon_nm=root_menu['x2flux.lon_nm'], \
                                                     lat_nm=root_menu['x2flux.lat_nm'], \
                                                     flux_nm=root_menu['x2flux.flux_nm'])
        
    
        # S6: create class for run_desc
        
        datapath=root_menu['run_desc.path']
        flnm=root_menu['run_desc.flnm']
        
        self.cl_run_desc=root_menu['run_desc.fclass'](flnm,\
                                                         datapath, \
                                                         yyyy, mm, dd,\
                                                         fopen=root_menu['run_desc.fopen'],\
                                                         fread=root_menu['run_desc.fread'],\
                                                         fclose=root_menu['run_desc.fclose'],\
                                                         fget=root_menu['run_desc.fget'])
                                                     
        
        # S7: state vector 
        fc=root_menu['stv.fclass']
        self.cl_stv=fc(yyyy,  mm, dd,\
                           root_menu['stv.max_step'],\
                           root_menu['stv.lag_window'],\
                           mod_bias=root_menu['stv.mod_bias'],\
                           mod_bias_err=root_menu['stv.mod_bias_err'])
        
        
                                                     
        

                                                  
 
