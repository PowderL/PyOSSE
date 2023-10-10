""" Class for GEOS-Chem profiles 

    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.10.28
    
    Classes:
    --------------------------------------------
    1. enr_slice_cl:  class derived from ctm_slice_cl to read and store GEOS-Chem ensemble run outputs
""" 
import numpy as npy
# from scipy import *

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.vertical_interp_m as vintpl_m
import ESA.util.otool_gcfile_io as ogco
import ESA.util.pres_m   as pm

import ctm_grid_3d as cgrd
import gc_grid_3d as gcgrd
import ctm_field_m as cfm
import ctm_slice_m as slm
import ctm_profile_m as cpm

import ESA.util.gc_ts_file_m as tsm

class enr_slice_cl(slm.ctm_slice_cl):
    
    """ class derived from ctm_slice_cl to read and store ensemble run outputs
    
    Members (overridden or new) 
    --------------------------------------------------------
    

     1.  fdesc: <gcfdesc_cl>: file access
     2.  fopen: <func>:  open file to read
     3.  fread:<func>: function to read data to fdesc and self.data_lst
     4.  fget:<func>: function to get selected profile 
     5.  fpres:<func>: function to calculate or fetch model pressure
     6.  flnm:<str>: file name (template) for GEOS-Chem outputs
     7.  datapath:<str>: path for GEOS-Chem outputs
     8.  fdiaginfo:<str>: file name (template) for diaginfo file
     9.  ftracerinfo:<str>:file name (template) for tracerinfo file
     10. yyyy:<int>: year
     11  mm:<int>: monthq
     12. dd:<int>: day
     
     
        
    
    Functions (new or overridden)
    ---------------------------------------
    
    1. __init__: initialization
    2. load_mod_output:read GEOS-Chem outputs from ensemble run into self.data_lst
    3. get_sel_profile: select data according to name, id, and time
    4. compute_mod_pres: calculate model pressures from surface pressure 
       or model pressure read in from disk file 
    5. set_mod_pres: set pressure grid from given pressure or data read in from disk file 
    
    """
    
    def __init__(self,datapath, filename, \
                     ftracerinfo,\
                     fdiaginfo,\
                     name, \
                     fopen=tsm.open_enr_file,\
                     fread=tsm.setup_enr_daily_profile,\
                     fget=tsm.get_mod_gp,\
                     fpres=tsm.get_mod_pres,\
                     **keywords):
        
        
        """Initialize 
        
        Inputs:
        --------------------------------
        
        1. datapath:<str>: path 
        2. flnm:<str>: file name
        3, ftracerinfo:<str>: name for tracerinfo file
        4. fdiaginfo:<str>: name for diaginfo file
        5. geos_ver:<int>:geos-chem version (4 or 5)
        6. use_reduced:<int>: 1 if reduced vertical grid is used
        
        7. fopen: <func>:  open file to read
        8. fread:<func>: function to read data to fdesc and self.data_lst
        9. fget:<funct>: function to get selected profile 
        10. keywords: <dict>: dictionary for attributes
        
        """
        
        # #S1: initialize parent class
        
        slm.ctm_slice_cl.__init__(self,name, [])
        
        
        
        # #S2: set values
        
        self.fdesc=None
        self.fopen=fopen
        self.fread=fread
        self.fget=fget
        self.flnm=filename
        self.datapath=datapath
        self.fdiaginfo=fdiaginfo
        self.ftracerinfo=ftracerinfo
        self.fpres=fpres

        self.yyyy=2009
        self.mm=1
        self.dd=1
        self.attr_dict.update(keywords)
        
        
    def load_mod_output(self, yyyy, mm, dd, \
                            enr_yst_lst,\
                            enr_dst_lst,\
                            enr_step_lst,\
                            enr_em_st_lst,\
                            enr_em_end_lst,\
                            olon=None, \
                            olat=None,\
                            gpname_lst=['PSURF', 'CO2'],\
                            append_data=False,\
                            **keywords):
        
        """                                        
        read outputs from ensemble run into self.data_lst
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. enr_yyyy_lst:<int>: year for ensemble run start.
        3. enr_doy_lst:<int>:  doy for ensemble run start.
        4. enr_step_lst:<int>: step for ensemble run.
        5. enr_em_st_st:<int>: firt one in the ensemble. 
        6. enr_em_end_lst:<int>: last one in the ensemble.
        7. olon:<array>: longitude
        8. olat:<array>:latitude 
        9. spname:<str>: tracer name for surface pressure outputs
        10.spcategory:<str>: tracer category for surface pressure outputs
        11. gpname_lst:<list, t:str>: name of tracers to be read in 
        12. append_data:<T/F>: if True, read data will be appended to list end
        
        13. keywords:<dict>: extra inputs:
        ---reserved keywords:
        --->1. categorys: <list, t:str>: category for variables to be read 
        --->2. tracers:<list, t:int>: list for tracer ID 
        --->3. taus: <list, t:float>: list for starting time 
        
        

        Returns:
        ==========================================================
        

        Notes:
        ----------------------------------------
        1. if olon and olat are not given, self.lon, and self.lat will be used 
        --- if olon and olat are given, self.lon, and self.lat will be replaced 
        
        """
        
        # #S1: delete previous file descript
            
        if (self.fdesc<>None):
            del self.fdesc
            self.fdesc=None
        
        # ##T: if not append, re-set the list 
        
        if (not append_data):
            if (self.data_lst<>None):
                del self.data_lst
                self.data_lst=[]
        else:
            if (self.data_lst==None):
                self.data_lst=[]
        
            
         
        # S2 set up new fdesc
        
        if (gpname_lst==None):
            tranames=None
        else:
            tranames=gpname_lst
        
        if ('categorys' in keywords):
            categorys=keywords['categorys']
            del keywords['categorys']
        else:
            categorys=None
            
            
        if ('tracers' in keywords):
            tracers=keywords['tracers']
            del keywords['tracers']
        else:
            tracers=None
        
        if ('taus' in keywords):
            taus=keywords['taus']
            del keywords['taus']
        else:
            taus=None
            
        # S2: open file
        
        self.fdesc=self.fopen(self.datapath, self.flnm,\
                                  self.ftracerinfo,  \
                                  self.fdiaginfo,\
                                  yyyy, mm,dd,\
                                  categorys=categorys,\
                                  tranames=tranames,\
                                  tracers=tracers,\
                                  taus=taus,\
                                  is_enr_file=True,\
                                  **keywords)
            
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd                                      
        
        if (olon<>None):
            self.lon=olon
        
        if (olat<>None):
            self.lat=olat
        
        
            
        # S3: read in data 
        
        self.fdesc=self.fread(self.fdesc, \
                                  self.yyyy, \
                                  self.mm, self.dd, \
                                  self.lon, self.lat,\
                                  enr_yst_lst,\
                                  enr_dst_lst,\
                                  enr_step_lst,\
                                  enr_em_st_lst,\
                                  enr_em_end_lst)
        
            
        # S4: collect data, and append to the list 
        
        cprof_lst=[]            
        
        if (tranames <>None):
            #  '##put things under their name 
            for gpname in tranames:
                found_data_lst=self.fget(self.fdesc, gpname)
                
                if (found_data_lst<>None):
                    found_data_lst=npy.array(found_data_lst)
                    cprof=cpm.ctm_profile_cl(gpname, found_data_lst, \
                                                 lon=self.lon, lat=self.lat, \
                                                 traname=gpname)
                    cprof_lst.append(cprof)
        else:
            # ##F:  put every one into the list      
        
            idx=0
            
            for iref in self.fdesc.ref_lst:  # loop over data set 
                
                bp_lst=self.fdesc.data_lst[idx]
                
                for bpdata in bp_lst:
                    gpname=bpdata.get_attr('traname')
                    prof_attr=dict(bpdata.attr_dict)
                    
                    # remove ot_type to avoid possible conffict. 
                    
                    if ('ot_type' in prof_attr):
                        del prof_attr['ot_type']

                    if ('name' in prof_attr):
                        del prof_attr['name']

                    prof_attr.update({'ref':iref})
                    
                    # ##T: construct profile 
                    
                    cprof=cpm.ctm_profile_cl(gpname, bpdata.data, \
                                                 lon=self.lon, lat=self.lat, **prof_attr)
                            
                    cprof_lst.append(cprof)
                idx=idx+1
        
        self.data_lst=self.data_lst+cprof_lst
    
    def get_sel_profile(self, name, sel_id, time, **keywords):
        
        """
        
        select data according to name, id, and time
        
        Inputs:
        
        ======================================================
        1.name:<str>: name of the profile
        2.sel_id:<int>: index for selected profile (see Notes). 
        3.time:<array>: time (not used for GEOS-Chem outputs)
        4.keywords:<dict>: extra inputs

        Returns:
        ======================================================
        1. out_lst:<list/array>: select data set (Note 1)
        
        Notes:
        ===================================================
        1. When sel_id is set, out_lst could be:
        ---if cprof.data is in shape of [nobs, nlvl, ntracer (nem)], out_lst=cprof.data[:,:, sel_id].
        ---if cprof.data is given as [nobs, nlvl],  out_lst is cprof.data when  
        
        """
        
        out_lst=list()
        
        idx=0
        # search over the data lst 
        
        for cprof in self.data_lst:
            
            if (cprof.name==name):
                dims=npy.shape(cprof.data)
                
                if (len(dims)==3):
                    # #C: cprof.data given as [nobs, nlvl, ntracer (nem)]
                    
                    if (sel_id<>None):
                        return cprof.data[:,:, sel_id]
                    else:
                        # return the full set 
                        
                        return cprof.data
                    
                else:
                    # #C: data given as [nobs, nlvl]
                    
                    if (idx==sel_id):
                        return cprof.data
                    elif (sel_id==None):
                        out_lst=out_lst+[cprof.data]
                idx=idx+1
                    
        
        if (sel_idx==None):
            return out_lst
        else:
            return None
    

    def compute_mod_pres(self, **keywords):
        """
        calculate or fetch model pressure
        
        Inputs:
        --------------------------------------------
         
        """
        
        mod_pres=self.fpres(self.fdesc, **keywords)
        return mod_pres
    
    def set_mod_pres(self, mod_pres=None, is_in_log=False, do_reverse=True, **keywords):
        """
        set up pressure grid 
        Inputs:
        ---------------------------------------------
        1. mod_pres:<array>: model pressure
        2. is_in_log:<T/F>: True if model pressure if given as log 
        3. do_reverse:<T/F>: True if model pressure is in descending order
        
        """
        
        if (mod_pres==None):
            mod_pres=self.fpres(self.fdesc, **keywords)
        
            if (mod_pres<>None):
                self.set_pressure(mod_pres, is_in_log, do_reverse=do_reverse, mask_val=self.mask_val)
            
        return mod_pres
    


                                  
if (__name__=="__main__"):
    
    print 'For test, see ctm_world_m.py'
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    
    yyyy, mm,dd=2009, 2, 20
    filename, fdiaginfo, ftracerinfo="", "", ""
    name='test'
    
    gsf=enr_slice_cl(datapath, filename, \
                         fdiaginfo,\
                         ftracerinfo,\
                         name, \
                         fopen=tsm.open_enr_file,\
                         fread=tsm.setup_enr_daily_profile,\
                         fget=tsm.get_mod_gp,\
                         fpres=tsm.get_mod_pres\
                         )
    
    
    enr_yst_lst=[2009, 2009]
    enr_dst_lst=[1,  1]
    enr_step_lst=[1, 1]
    enr_em_st_lst=[1, 81]
    enr_em_end_lst=[81, 146]
    olon=npy.arange(-40, 40, 10)
    olat=npy.arange(-40, 40, 10)
    olat=olat+10

    # read in the data
    
    gsf.load_mod_output(yyyy, mm, dd, \
                            enr_yst_lst,\
                            enr_dst_lst,\
                            enr_step_lst,\
                            enr_em_st_lst,\
                            enr_em_end_lst, \
                            olon=olon, olat=olat,\
                            gpname_lst=['PSURF', 'CO2'])
    
    # mod_pres=gsf.pressure
    
    mod_co2_prof_lst=gsf.find_profile_cl(traname='CO2')
    
    print len(mod_co2_prof_lst)
    mod_co2_prof=mod_co2_prof_lst[0]
    print npy.shape(mod_co2_prof.data)
    
    pres=gsf.set_mod_pres()
    print npy.shape(gsf.pressure)
    print gsf.pressure[0,:]
    
    
    
