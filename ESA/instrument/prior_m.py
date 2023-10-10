""" Class for pre-calculated orbit 
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. prior_cl: Class for reading and filtering pre-calculated priors 
    
"""
    
import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.vertical_profile as vpf  # vertical regrid
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import prior_file_m as apr_io

import ESA.util.vertical_interp_m as vintpl_m



class prior_cl:

    """
    Class for reading priors 
    

    Members:
    ====================================================================
    1. flnm    :<str>: Name of prior file (reserved for future use)
    2. datapath:<str>: file path 
    3. yyyy    :<int>: year 
    4. mm      :<int>: month 
    5. dd      :<int>: day 
    
    6. fopen   :<func>: construct class file_desc_cl for accessing prior file  (see Note 5)
    7. fread   :<func>: read prior data into fdesc (see Note 6)
    8. fclose :<func>: prior file close function  (see Noye 7)
    9. fget :<func>:  get prior profile at given locations

    10. var_name_dict:<dict>: translation between column name in the prior file 
    --- and variable name in class 
    
    11. fio_keywords:<dict>: extra parameter for opening and reading prior files 
   
    12. fdesc:<file_desc_cl>: file access class
    
    

    Functions
    ===============================================================
    1. __init__: initialization
    2. load_prior_from_file: read prior from file 
    3. interp_prior: interpolate prior profile to observation location
    4. get_prior:  get prior data at observation time and location
    5. get_var: get variable stored in self.fdesc

    Notes:
    ================================================================
    
    1) After prior file is read, self.orb and self.var_name_dict must  
    contain the following entries:
    --- 1) 'pres' for pressure
    --- 2) 'vmr' for mixing ratio

    
         
  

    
    """
    
    def __init__(self,   datapath,\
                     flnm,\
                     yyyy, \
                     mm, \
                     dd,\
                     gpname,\
                     fopen=apr_io.open_prior_file,\
                     fread=apr_io.read_prior_file,\
                     fclose=apr_io.close_prior_file,\
                     fget=apr_io.get_prior_data,\
                     varname_dict=apr_io.apr_varname_dict,\
                     mask_val=oob.fill_val,\
                     fio_keywords={} ):
        
        """
        Initialization
        
        Inputs:
        ---------------------------------------------
        1. datapath:<str>: file path 
        2. flnm    :<str>: Name of prior file (reserved for future use)
        3. yyyy    :<int>: year 
        4. mm      :<int>: month 
        5. dd      :<int>: day 
        6. gpname  :<str>: name of the tracer
        7. fopen   :<func>: prior file open function (Note 5). 
        8. fread   :<func>: prior file read function (Note 6).
        9. fclose :<func>: prior file close function (Note 7).
        10. fget:  <func>: get prior data
        11. varname_dict:<dict>: translation from variable name to column name in prior file 
        
        12. fio_keywords:<dict>: extra parameter for opening and reading prior files 
        

        Notes:
        ----------------------------------------------------------------------
        1. variable names used to access data
        
        --->pres: pressure
        --->vmr:  profile for mixing ratio 
        
        
        """
        
        # S1 initialize members 

        
        ## time 
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        ## prior file 
        self.flnm=flnm
        self.datapath=datapath
        self.gpname=gpname
        
        # file function 
        
        self.fopen=fopen
        self.fread=fread
        self.fclose=fclose
        self.fget=fget
        self.fio_keywords={}
        self.fdesc=None # file access class
        
        
        self.varname_dict=varname_dict
        
        self.mask_val=mask_val
        

                
        # S2 read prior file 
        if (yyyy<>None and yyyy>0):
            self.load_prior_from_file(yyyy, mm, dd, \
                                          **self.fio_keywords)
        
            
    
    def load_prior_from_file(self, yyyy, mm, dd, \
                                 **keywords):
        
        """ read prior from file 
        
        Inputs:
        ---------------------------------------------
        1. yyyy    :<int>: year 
        2. mm      :<int>: month 
        3. dd      :<int>: day 
        4. keywords :<dict>: extra parameters for file access 
                
        Returns
        -----------------------------------------------
        1. len(fdesc.dtable):<int>: number of record read in 
        
        
        """
        
        # S1 reset the class members 
        
        ## time 
        
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        # S1 open the file 
        
        doy=tm.day_of_year(yyyy, mm, dd)
        self.fclose(self.fdesc)
        self.fdesc=self.fopen(self.flnm, self.datapath, \
                                  self.gpname,\
                                  yyyy, mm,dd,\
                                  **keywords)
        
        # S2 read into prior data to fdesc.dtable 
        
        self.fdesc=self.fread(self.fdesc, **keywords)
        return self.fdesc
    
    def get_var(self, varname):
        
        """
        get the variable stored in self.fdesc

        Inputs:
        ------------------------------------
        1. varname:<str>: variable name 
        

        Returns:
        -------------------------------
        1. var:<array> variable 
        
        
        """
        if (varname in self.varname_dict):
            varname_fl=self.varname_dict[varname]
            var=self.fdesc.dtable[varname_fl]
        else:
            var=self.fdesc.dtable[varname]
            
        return var

        
    def interp_prior(self, olon, olat):
        """
        
        interpolate prior profile to observation location
        
        Inputs:
        ---------------------------------------------------------------
        
        1. olon    :<array, (nobs)>: observation (sampling) longitude
        2. olat    :<array, (nobs)>: observation (sampling) latitude 
        
        Returns:
        ---------------------------------------------------------------
        1. mod_pres:<array, (nobs, nz)>: press at observation location
        2. mod_prof:<array, (nobs, nz)>: (original) prior profile at observation location

        """
        # S1: get pressure and profile
        
        # #c:get pressure 

        pres=self.get_var('pres')
        
        # #c: get  profile 
    
        prof=self.get_var('vmr')
        
        # S2: interpolate to locations q
        

        nobs=npy.size(olon)
        nz=npy.size(pres)
        
        # #T1: consruct mod_pres and mod_prof from prior values
        
        if (oob.get_ot_type(olon)==oob.ot_array):
            # #c: if olon is array, mod_pres and mod_prof are in the shape of [nobs, nz]
            
            mod_pres=npy.zeros([nobs, nz], float)
            mod_prof=npy.zeros([nobs, nz], float)
            
            mod_pres[:,:]=pres[npy.newaxis, :]
            mod_prof[:,:]=prof[npy.newaxis, :]
        else:
            
            mod_pres=pres
            mod_prof=prof
        
        return mod_pres, mod_prof
    
    
    def get_prior(self, yyyy, mm, dd,\
                      olon, olat, opres,\
                      is_opres_in_log=False, \
                      do_ob_reverse=False,\
                      **fio_keywords):
        
        """    get prior data at observation time and locations
    
        
        Inputs:
        ---------------------------------------------------------------
        
        1. yyyy, mm, dd:<int>: year, month, day
        2. olon    :<array, (nobs)>: observation (sampling) longitude
        3. olat    :<array, (nobs)>: observation (sampling) latitude 
        3. opres   :<array, (nobs, obs_nz)>: observation (retrieval) pressure grid 
        4. is_opres_in_log:<T/F>: True if the pressure is given in log10. 
        5. do_ob_reverse: <T/F>: True, if opres is given in decending order. 
        6. fio_keywords :<dict>: extra parameters for file access 
        
        Returns:
        -------------------------------------
        1. opres:<array (nobs, obs_nz)>: observation (retrieval) pressure grid 
        2. oprof::<array (nobs, obs_nz)>: observation prior profile 
        3. col_obs::<array (nobs)>: observation column values
        

        Notes:
        -------------------------------------
        1, oprof are calculated through horizontal and vertical interpolations
        
        
        """
    
        # S1: reload the data
    
        
        self.fdesc=self.load_prior_from_file(yyyy, mm, dd, \
                                                 **fio_keywords)
        
        
        
        # S2: get prior pressure and profile at observation location
        
        mod_pres, mod_prof=self.interp_prior(olon, olat)
        
        # S3: create vertical interpolation class
    
        do_mod_reverse=False
    
        if (mod_pres[0,2]<mod_pres[0,1]):
            do_mod_reverse=True
        
        is_mod_in_log=True
        if (npy.max(mod_pres[0,:])>10):
            is_mod_in_log=False
            
        cl_vitpl=vintpl_m.vintpl_cl(mod_pres, is_mod_in_log, do_mod_reverse)
    
        # #T3: initialize interpolation to observation grid 
    
        
        cl_vitpl.init_interp(opres, is_in_log=is_opres_in_log, do_ob_reverse=do_ob_reverse)
        
        # S4: interpolate profiles to observation (retrieval) pressure grid 
        
        oprof=cl_vitpl.interpolate_mod_prof(mod_prof)
    
        # S5: get column 
        
        col_obs=cl_vitpl.get_ob_column(oprof)
        
        return opres, oprof, col_obs


    

    
        
if (__name__=='__main__'):
    datapath="/home/lfeng/local_disk/otool_data/clim_dat/"
    
    flnm=""
    yyyy, mm, dd=2009, 2,1
    cl_apr=prior_cl(datapath,\
                        flnm,\
                        yyyy, \
                        mm, \
                        dd,\
                        'CO2',\
                        varname_dict=apr_io.apr_varname_dict,\
                        fopen=apr_io.open_prior_file,\
                        fread=apr_io.read_prior_file,\
                        fclose=apr_io.close_prior_file,\
                        fget=apr_io.get_prior_data,\
                        mask_val=oob.fill_val,\
                        fio_keywords={} )
    
    
    
    
    
            
                  
