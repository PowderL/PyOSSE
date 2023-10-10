""" Class for reading GEOS-Chem outputs to profiles 

    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.10.28
    
    Classes:
    --------------------------------------------
    1. ts_slice_cl:  class derived from ctm_slice_cl to read and store time-serie outputs

""" 

import numpy as npy
# from scipy import *

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.vertical_interp_m as vintpl_m
import ESA.util.otool_gcfile_io as ogco

import ctm_grid_3d as cgrd
import gc_grid_3d as gcgrd
import ctm_field_m as cfm
import ctm_slice_m as slm
import ctm_profile_m as cpm

import ESA.util.gc_ts_file_m as tsm

class gc_slice_cl(slm.ctm_slice_cl):
    
    """ class derived from ctm_slice_cl to read and store GEOS-Chem outputs
    
    Members (overridden or new) 
    --------------------------------------------------------
    

     1.  fdesc: <gcfdesc_cl>: file access
     2.  fopen: <func>:  open file to read
     3.  fread:<func>: function to read data to fdesc and self.data_lst
     4.  fget:<func>: function to get selected profile 
     5.  fpres:<func>: function to calculate or find model pressure
     6.  flnm:<str>: file name (template) for GEOS-Chem outputs
     7.  datapath:<str>: path for GEOS-Chem outputs
     8.  fdiaginfo:<str>: file name (template) for diaginfo file
     9.  ftracerinfo:<str>:file name (template) for tracerinfo file
     10. yyyy:<int>: year
     11  mm:<int>: month
     12. dd:<int>: day
     13. sp_name:<str>: surface pressure name 
     
          
     
        
    
    Functions (new or overridden)
    ---------------------------------------
    
    1. __init__: initialization
    2. load_model_output: read model outputs into slice
    3. get_sel_profile: find the selected data
    4. compute_mod_pres: calculate model pressures from surface pressure 
       or model pressure read in from disk file 
    5. set_mod_pres: set pressure grid from given pressure or data read in from disk file 
    

    """
    
    def __init__(self,datapath, filename, \
                     ftracerinfo,\
                     fdiaginfo,\
                     name, \
                     fopen=tsm.open_ts_file,\
                     fread=tsm.setup_daily_profile,\
                     fget=tsm.get_mod_gp,\
                     fpres=tsm.get_mod_pres,\
                     ext="",\
                     sp_name='PSURF',\
                     **keywords):
        
        
        """Initialize 
        
        Inputs:
        --------------------------------
        
        1. datapath:<str>: path 
        2. flnm:<str>: file name
        3, ftracerinfo:<str>: name for tracerinfo file
        4. fdiaginfo:<str>: name for diaginfo file
        
        6. fopen: <func>:  open file to read
        7. fread:<func>: function to read data to fdesc and self.data_lst
        8. fget:<func>: function to get selected profile 
        9. fpres:<func>: function to calculate pressure 
        
        10. ext:<str>: extension for file name 
        
        11. keywords: <dict>: dictionary for attributes
        
        """
        
        # #S1: initialize parent class
        
        slm.ctm_slice_cl.__init__(self,name, [])
                
        
        # #S2: set values
        
        self.fdesc=None
        self.fopen=fopen
        self.fread=fread
        self.fget=fget
        self.fpres=fpres
        
        self.flnm=filename
        self.datapath=datapath
        self.fdiaginfo=fdiaginfo
        self.fext=ext
        self.sp_name=sp_name
        
        self.ftracerinfo=ftracerinfo
        self.yyyy=2009
        self.mm=1
        self.dd=1
        
        # #S3: get coefficients for pressure levels 
    
        self.attr_dict.update(keywords)
    
    def get_attr(self, varname):
        if (varname in self.attr_dict):
            return self.attr_dict[varname]
        else:
            return None
    
    def load_mod_output(self, yyyy, mm, dd, \
                            olon=None, \
                            olat=None,\
                            gpname_lst=['SPURF', 'CO2'],\
                            append_data=False,\
                            **keywords):
        
        """                                        
        read GEOS-Chem outputs from ensemble run into self.data_lst
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. olon:<array>: longitude
        3. olat:<array>:latitude 
        4. spname:<str>: tracer name for surface pressure outputs
        5.spcategory:<str>: tracer category for surface pressure outputs
        6. gpname_lst:<list, t:str>: name of tracers to be read in 
        7. append_data:<T/F>: if True, read data will be appended to list end
        8. keywords:<dict>: extra inputs:
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
                                  ext=self.fext)
        
        # S4: collect data, and append to the self.data_lst 
        
        cprof_lst=[]            
        
        if (tranames <>None):
            # 
            for gpname in tranames:
                found_data_lst=self.fget(self.fdesc, gpname)
                if (found_data_lst<>None):
                    cprof=cpm.ctm_profile_cl(gpname, found_data_lst, \
                                                 lon=self.lon, lat=self.lat, \
                                                 traname=gpname)
                    cprof_lst.append(cprof)
        else:
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
         

        outputs:
        ---------------------------------------------
        1. mod_pres:<array, ([nx], [ny], [nlvl])>: model pressure
        
        """
        
        new_keywords=dict(keywords)
        
        if ('traname' in new_keywords):
            traname=new_keywords['traname']
            del new_keyords['traname']
            
        else:
            traname=self.sp_name
        
        if ('use_sp' in new_keywords):
        
            use_sp=new_keywords['use_sp']
            del new_keyords['use_sp']
            if (use_sp==0):
                use_sp=False
            
        else:
            use_sp=True
        
        
            
        mod_pres=self.fpres(self.fdesc, traname=traname, \
                                use_sp=use_sp, \
                                **new_keywords)
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
    datapath='/home/lfeng/local_disk/otool_data/enkf_rerun/'
    
    yyyy, mm,dd=2009, 2, 20
    filename, fdiaginfo, ftracerinfo="", "", ""
    name='test'
    
    gsf=gc_slice_cl(datapath, filename, \
                        fdiaginfo,\
                        ftracerinfo,\
                        name, \
                        fopen=tsm.open_ts_file,\
                        fread=tsm.setup_daily_profile,\
                        fget=tsm.get_mod_gp,\
                        fpres=tsm.get_mod_pres,\
                        ext='.ST001.EN0001-EN0002'\
                        )
    
    olon=npy.arange(-40, 40, 10)
    olat=npy.arange(-40, 40, 10)
    olat=olat+10

    # read in the data
    
    gsf.load_mod_output(yyyy, mm, dd, \
                            olon=olon, olat=olat,\
                            gpname_lst=['PSURF', 'CO2'])
    
    # mod_pres=gsf.pressure
    
    mod_co2_prof_lst=gsf.find_profile_cl(traname='CO2')
    
    print len(mod_co2_prof_lst)
    mod_co2_prof=mod_co2_prof_lst[0]
    print npy.shape(mod_co2_prof.data)
    mod_co2_prof2=gsf.get_sel_profile('CO2', 0, 0)
    
    pres=gsf.set_mod_pres()
    
    
    
    print npy.shape(mod_co2_prof2)
    print npy.shape(gsf.pressure)
    print gsf.pressure[0,:]
    
    
    pres=gsf.set_mod_pres(use_sp=0)
    print npy.shape(pres)
    print pres
    
    
