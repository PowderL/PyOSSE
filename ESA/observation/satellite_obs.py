"""
Functions for reading observations from netcdf files.  

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.21
    
    
    Classes:
    ===============================================
    1. satellite_obs_cl: class for reading satellite observations 
    


"""

import os
import pylab as plb
import numpy as npy 
import ESA.util.time_module as tm
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.otool_var_io as ovio
import ESA.util.otool_obj as oob 
import ESA.util.message_m as msm




## <PARAMTERS>

# names for variable in observation file 
# they are the names to be read


sat_obs_varname_lst=['lon', 'lat', \
                         'time', 'obs', \
                          'err', 'corrected_obs', \
                          'aod', 'osza', \
                          'psurf_ap', \
                          'psurf_ret', \
                          'master_flag', \
                          'cloud_flag', \
                          'land_fraction', \
                          'obs_ak', \
                          'obs_apr', \
                          'obs_pres', \
                          'ostype']



## translate varname to names in observation file 


sat_obs_varname_dict={'lon':'lon', \
                          'lat':'lat', \
                          'time':'time', \
                          'obs':'obs', \
                          'err':'err', \
                          'corrected_obs':'corrected_obs', \
                          'aod':'aod', 'osza':'osza', \
                          'psurf_ap':'psurf_ap', \
                          'psurf_ret':'psurf_ret', \
                          'master_flag':'master_flag', \
                          'cloud_flag':'cloud_flag', \
                          'land_fraction':'land_fraction', \
                          'obs_ak':'obs_ak', \
                          'obs_apr':'obs_apr', \
                          'obs_pres':'obs_pres', \
                          'ostype':'ostype'}






class satellite_obs_cl:
    """ class for reading satellite observations 
    
    Members:
    ===================================================
    1. fdesc: <ncfile_desc_cl>: access to netcdf file 
    2. nvar:<int>: number of variables
    3. nobs:<int>: number of observations
    4. var_lst:<list, t:io_var_cl>: list of variables
    
    5. attr_dict:<dict> : dictionary for attributes. 
    6. viewtype:<str>: View type 
    7. viewmode:<str>: viewmode
    
    8. flnm:<str>: observation file name (format)
    9. datapath:<str>: obs_datapath
    
    10. varname_dict:<dict>: (varname_dict) # dictionary
    11. mask_val:<float>: filling value for missing or bad data 
    12. block_obs_flag:<list>: list of observation to be block
    13. obs_platform:<str>: observation platform
    
    

    Functions:
    ========================================================
    1. __init___: initialization
    2. read_obs: read variable into self.var_lst 
    3. make_var_list: set up list for var to be read
    4. get_obs_list:list of observation to be used
    5. filter_obs: filtering observations. 
    6. __getitem__: overriden index function


    """
    
    def __init__(self, \
                     obs_datapath, \
                     flnm, \
                     viewtype, \
                     sat_id, \
                     viewmode='',\
                     varname_lst=sat_obs_varname_lst,\
                     varname_dict=sat_obs_varname_dict, \
                     mask_val=oob.fill_val,\
                     **keywords):
   
        """
        1. obs_datapath:<str>: obs_datapath
        
        2. flnm:<str>:  file name  (format) 
        3. viewtype:<str>: view type 
        4. sat_id:<str>:ID of satellite
        5. viewmode:<str>: viewmode
        6. varname_lst:<list>: nam of variables to be read 
        7. varname_dict:<dict>: dictionary for translating variable used in the code to variable stored in the NetCDF file 
        
        8. mask_val:<float>: filling value for missing or bad data. 
        
        

        """

        # S1: set the attribute 
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_obs_op})
        self.attr_dict.update(keywords)
        self.viewtype=viewtype
        self.viewmode=viewmode
        
        
        # S2: set up fdesc to access the file 
        
        self.flnm=flnm  # format 
        
        self.varname_dict=dict(varname_dict) # dictionary
        self.nvar=None  
        self.nobs=0
        self.var_lst=None
        self.mask_val=mask_val
        self.make_var_list(varname_lst)
        
        self.datapath=obs_datapath
        
        
        self.fdesc=ncfio.ncfile_desc_cl(var_lst=self.var_lst, \
                                            varname_dict=self.varname_dict,\
                                            mask_val=self.mask_val,\
                                            read_to_array=True)
        
        self.fdesc.set_file_path(self.datapath)
        self.fdesc.set_filename_format(self.flnm)

        
        
        self.block_obs_flag=None
        self.sat_id=sat_id
        self.obs_platform="satellite"
        
    
    def read_obs(self, yyyy, mm,dd, **keywords):
        """
        Read into var list
        
        Inputs:
        ============================================================
        1. viewmode:<str>: view mode 
        2. yyyy, mm, dd:<int>: date 
        3. keywords:<dict>: dictionary for extra input 
    

        Returns:
        ============================================================
        1. self.nobs:<int>: number of observations found 
    

        Notes:
        ===========================================================
        
        1. observation variable assumed to be in shape of [nobs, ....]
        
        """
        
        if ('viewmode' in keywords):
            self.viewmode=keywords['viewmode']
        
        if ('viewtype' in keywords):
            self.viewtype=keywords['viewtype']
        
        
        doy=tm.day_of_year(yyyy, mm, dd)
        # S1: construct file name 
        
        full_flnm=self.fdesc.construct_filename(XVIEWTYPEX=self.viewtype,\
                                                    XVIEWMODEX=self.viewmode, \
                                                    XYYYYX=yyyy, XMMX=mm, XDDX=dd, XDOYX=doy)
        
        
        # read in to variables if file existing 
        
        print 'read obs from:', full_flnm
        
        if (os.path.isfile(full_flnm)):
            # print 'I am here'
            self.fdesc.read_var_from_file(full_flnm)
            # S2: var_lst will be shared bewteen fdesc and its 
            self.var_lst=self.fdesc.var_lst
            
            # print len(self.var_lst)
            
            for gvar in self.var_lst:

                dims=gvar.shape()
                print dims
                
                if (len(dims)>=1):
                    nobs=dims[0]
                    break
            
            self.nobs=nobs
            
        else:
            msg='file not found: '+full_flnm
            msm.show_err_msg(msg)
            self.nobs=0
        
        return self.nobs
   

    def make_var_list(self, varname_lst, **keywords):
        """
        Make a list for var to be read 

        Inputs:
        ======================================================
        1. varname_lst: <list, t:str>: list of variables to read in 
        2. keywords: <dict>: extra inputs
        ---if varname:io_var_cl in keywords, it will be used to replace the stand setting. 


        Returns:
        ================================================
        1, self.nvar:<int>: number of var 
        Notes:
        """
        
        self.varname_lst=list(varname_lst)
        if (self.var_lst<>None):
            del self.var_lst
        
        self.nvar=len(self.varname_lst)
        # loop over vars

        self.var_lst=[]
        
        for varname in self.varname_lst:
            
            if (varname in keywords):
                new_var=keywords[varname]
            else:
                if (varname in self.varname_dict):
                    # translate to name stored in netcdf file
                    varname=self.varname_dict[varname]
                
                new_var=ovio.io_var_cl(varname, 'f', [], None)
            
            self.var_lst.append(new_var)
        print 'nvar fatre create', len(self.var_lst)
        print 'varname_lst', self.varname_lst
        
        return self.nvar
    
    def get_obs_list(self, **keywords):
        
        """
        index list for filtering the observations 
        (to be overriden by user)
        
        Inputs:
        ============================================================
        1. keywords:<dict>: extra inputs 
        
        """
        usd_idx_lst=range(self.nobs)
        
        
        return usd_idx_lst

    def filter_obs(self, use_idx_lst, **keywords):
        """
        Remove bad obs 
        
        Inputs:
        =========================================================
        1. use_idx_lst:<list, t:int>: index for obs to be kept 


        Outputs:
        ===========================================================
        1. self.nobs:<int>: number of observations left 
        
        """
        
        self.nobs=len(use_idx_lst)
        for gvar in self.var_lst:
            gvar.var=gvar.var[use_idx_lst]
            
        return self.nobs
    

    def __getitem__(self, index):
        
        """
        Over-ridding the index function 
        
        Inputs:
        =============================================================
        1. index:<int/str>: name of the list in the  
        
        Returns
        1. data:<array>: data 
        =============================================================
        
        """
        data=None
        gpvar=self.fdesc[index]
        if (gpvar<>None):
            data=gpvar.var
            
        return data 
    


#<<<TEST>>>

if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/obs/gosat_v210_obs/'
    viewtype='gosat_v29'
    flnm='XVIEWTYPEX.XYYYYXDXDOYX.nc'
    yyyy, mm, dd=2010, 03, 21
    sat_id=3
    obs_op= satellite_obs_cl(datapath, \
                                 flnm, \
                                 viewtype, \
                                 sat_id, \
                                 varname_lst=sat_obs_varname_lst,\
                                 varname_dict=sat_obs_varname_dict, \
                                 mask_val=oob.fill_val)
    
    nobs=obs_op.read_obs(yyyy, mm, dd)
    
    lon=obs_op['lon']
    lat=obs_op['lat']
    err=obs_op['err']
    obs=obs_op['obs']
    
    obs_ak=obs_op['obs_ak']
    obs_apr=obs_op['obs_apr']
    obs_ps=obs_op['obs_pres']
    
    
    
    print npy.shape(lon)
    print npy.shape(obs_ak)
    print npy.shape(obs)
    
    
    use_obs_idx=range(30)
    
    obs_op.filter_obs(use_obs_idx)
    print 'after filter'
    
    lon=obs_op['lon']
    lat=obs_op['lat']
    obs_ak=obs_op['obs_ak']
    obs_pres=obs_op['obs_pres']

    obs=obs_op['obs']

    print 'lon:', npy.shape(lon)
    print 'obs_ak:', npy.shape(obs_ak)
    print 'obs:', npy.shape(obs)
    
    print obs[0:5]
    print lat[0:5]
    print err[0:5]
    print obs_pres[0,:]
    print obs_pres[1,:]
    

    # read in another day 
    
    yyyy, mm, dd=2010, 8, 21
    
    nobs=obs_op.read_obs(yyyy, mm, dd)
    print nobs
    
