""" Class for  providing parameters used in generating virtual (dummy) observations

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. vobs_def_cl: Class for generating dummy observations 
    
"""
import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_menu_m as menu_m
import ESA.instrument.orbit_m as orb_m
import ESA.instrument.cloud_m as cld_m
import ESA.instrument.aod_m as aod_m
import ESA.instrument.ak_m as ak_m
import ESA.instrument.err_m as err_m


class vobs_def_cl:

    def __init__(self, \
                     yyyy, mm, dd, \
                     root_menu, **keywords):

        """
        Initialization 
        Inputs:
        ------------------------------------------
        1. yyyy, mm, dd:<int>: year, month, day
        2. root_menu:<menu_cl>:  menu for configuration 
        3. keywords:<dict>: extra inputs 
        
        """
        # S0 set values
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_def})
        
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        
        # S1: load menu 
        
        
        self.root_menu=root_menu

        # S2: set up orbit access 
        
        self.viewtype=self.root_menu['cfg.viewtype']
        self.viewmode=self.root_menu['cfg.viewmode']
        

        path=self.root_menu['orbit.path']


        
        self.cl_orb=self.root_menu['orbit.fclass'](self.viewtype, \
                                                     self.viewmode, \
                                                     self.root_menu['orbit.path'],\
                                                     self.root_menu['orbit.flnm'],\
                                                     yyyy, \
                                                     mm, \
                                                     dd,\
                                                     varname_dict=self.root_menu['orbit.dict'],\
                                                     fopen=self.root_menu['orbit.fopen'],\
                                                     fread=self.root_menu['orbit.fread'],\
                                                     fclose=self.root_menu['orbit.fclose'],\
                                                     fget=self.root_menu['orbit.fget'],\
                                                     mask_val=oob.fill_val, \
                                                     fio_keywords=self.root_menu['orbit.keywords'])
        
        # S3: Setup cloud climatology access 
        
        self.cl_cld=self.root_menu['cloud.fclass'](self.root_menu['cloud.path'],\
                                                       self.root_menu['cloud.flnm'], \
                                                       yyyy, mm, dd, \
                                                       varname_dict=self.root_menu['cloud.dict'],\
                                                       fopen=self.root_menu['cloud.fopen'], \
                                                       fread=self.root_menu['cloud.fread'],    \
                                                       fclose=self.root_menu['cloud.fclose'], \
                                                       fget=self.root_menu['cloud.fget'], \
                                                       fio_keywords=self.root_menu['cloud.keywords'])
        
        # S4: Setup aod accesss
        
        self.cl_aod=self.root_menu['aod.fclass'](self.root_menu['aod.path'], \
                                                     self.root_menu['aod.flnm'], \
                                                     yyyy, mm, dd, \
                                                     varname_dict=self.root_menu['aod.dict'],\
                                                     fopen=self.root_menu['aod.fopen'], \
                                                     fread=self.root_menu['aod.fread'],    \
                                                     fclose=self.root_menu['aod.fclose'], \
                                                     fget=self.root_menu['aod.fget'],\
                                                     fio_keywords=self.root_menu['aod.keywords'])

        # S5: Setup observation error access

        
        self.cl_err=self.root_menu['err.fclass'](self.viewtype, \
                                                     self.root_menu['err.path'], \
                                                     self.root_menu['err.flnm'], \
                                                     yyyy, mm, dd, \
                                                     viewmode_dict=self.root_menu['err.viewmode_dict'],\
                                                     surf_type_dict=self.root_menu['err.stype_dict'],\
                                                     fopen=self.root_menu['err.fopen'], \
                                                     fread=self.root_menu['err.fread'],    \
                                                     fclose=self.root_menu['err.fclose'], \
                                                     fget=self.root_menu['err.fget'],\
                                                     fio_keywords=self.root_menu['err.keywords'])
       
        # S6: setup averaging kernel access
        
        self.cl_ak=self.root_menu['ak.fclass'](self.viewtype, \
                                                   self.root_menu['ak.path'],\
                                                   self.root_menu['ak.flnm'], \
                                                   yyyy, mm, dd, \
                                                   viewmode_dict=self.root_menu['ak.viewmode_dict'],\
                                                   surf_type_dict=self.root_menu['ak.stype_dict'],\
                                                   fopen=self.root_menu['ak.fopen'], \
                                                   fread=self.root_menu['ak.fread'],    \
                                                   fclose=self.root_menu['ak.fclose'], \
                                                   fget=self.root_menu['ak.fget'],\
                                                   fio_keywords=self.root_menu['ak.keywords'])
        


        # S7: setup landcover  access
        
        self.cl_lc=self.root_menu['lc.fclass'](self.root_menu['lc.path'],\
                                                   self.root_menu['lc.flnm'], \
                                                   yyyy, mm, dd, \
                                                   varname_dict=self.root_menu['lc.dict'],\
                                                   fopen=self.root_menu['lc.fopen'], \
                                                   fread=self.root_menu['lc.fread'],    \
                                                   fclose=self.root_menu['lc.fclose'], \
                                                   fget=self.root_menu['lc.fget'],\
                                                   fio_keywords=self.root_menu['lc.keywords'])

        

    
        # S8: setup prior  access
        
        self.cl_apr=self.root_menu['apr.fclass'](self.root_menu['apr.path'],\
                                                     self.root_menu['apr.flnm'], \
                                                     yyyy, mm, dd, \
                                                     self.root_menu['apr.gpname'],\
                                                     varname_dict=self.root_menu['apr.dict'],\
                                                     fopen=self.root_menu['apr.fopen'], \
                                                     fread=self.root_menu['apr.fread'],    \
                                                     fclose=self.root_menu['apr.fclose'], \
                                                     fget=self.root_menu['apr.fget'],\
                                                     fio_keywords=self.root_menu['apr.keywords'])
        
        
        # S8: update attributes
        
                                                   
        # S9: sample functions 

        self.fcld_sample=self.root_menu['sample.fcld_sample']
        self.fcld_penalty=self.root_menu['sample.fcld_penalty']
        self.fcld_keywords=self.root_menu['sample.fcld_keywords']
        
        # #c:AOD

        self.faod_sample=self.root_menu['sample.faod_sample']
        self.faod_uplimit=self.root_menu['sample.faod_uplimit']
        self.faod_keywords=self.root_menu['sample.faod_keywords']

        self.flc_sample=self.root_menu['sample.flc_sample']
        self.flc_stype_lst=self.root_menu['sample.flc_stype_lst']
        self.flc_keywords=self.root_menu['sample.flc_keywords']
        

        
        # S10: grid 
        
        self.lon=self.root_menu['grid.lon']
        self.lat=self.root_menu['grid.lat']
        
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
    
    def set_attr(self, attr_name, attr_val):
        """
        set attribute
        Inputs:
        ---------------------------------------------
        1. attr_name:<str>: name
        2. attr_value:<ob>: attribute
        
        """
        self.attr_dict.update({keyname:keywords})
    def get_attr(self, attr_name):
        """
        Get attribute
        Inputs:
        ---------------------------------------------
        1. attr_name:<str>: name
        
        Returns:
        1. attr_value:<ob>: attribute
        
        """
        if (attr_name in self.attr_dict):
            return self.attr_dict[attr_name]
        else:
            return None

        
    
if (__name__=='__main__'):
    viewtype='aqua'
    viewmode='nadir'
    yyyy, mm, dd=2009, 3, 1
    menu_flnm='vob_def.cfg'
    odef=vobs_def_cl(viewtype, viewmode, \
                         yyyy, mm, dd, \
                         menu_flnm)
    
