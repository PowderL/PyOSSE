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


class sobs_def_cl:
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

        
        self.viewtype=self.root_menu['cfg.viewtype']
        self.viewmode=self.root_menu['cfg.viewmode']
        # S2: create observation class
                
        
        self.cl_obs=self.root_menu['sat_obs.fclass'](self.root_menu['sat_obs.path'], \
                                                         self.root_menu['sat_obs.flnm'], \
                                                         self.viewtype, \
                                                         self.root_menu['sat_obs.sat_id'], \
                                                         varname_lst=self.root_menu['sat_obs.varname_lst'],\
                                                         varname_dict=self.root_menu['sat_obs.varname_dict'], \
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
        
        
    
        # S3: class for ensemble profiles 
        
        self.cl_enr_prof=self.root_menu['enr_prof.fclass'](self.root_menu['enr_prof.path'],\
                                                               self.root_menu['enr_prof.flnm'], \
                                                               self.root_menu['enr_prof.ftracerinfo'],\
                                                               self.root_menu['enr_prof.fdiaginfo'],\
                                                               self.root_menu['enr_prof.vname'], \
                                                               fopen=self.root_menu['enr_prof.fopen'],\
                                                               fread=self.root_menu['enr_prof.fread'],\
                                                               fget=self.root_menu['enr_prof.fget'],\
                                                               fpres=self.root_menu['enr_prof.fpres'])
        
        
        
