""" Class for pre-calculated orbit 
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. orbit_cl: Class for reading and filtering pre-calculated orbits 
    
"""
    
import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.vertical_profile as vpf  # vertical regrid
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import orbit_file_m as orb_io




class orbit_cl:

    """
    Class for reading and filtering pre-calculated orbits 
    

    Members:
    ====================================================================
    1. viewtype:<str>: instrument type
    2. viewmode:<str>: view mode such as 'nadir', 'glint' etc
    3. flnm    :<str>: Name of orbit file (reserved for future use)
    4. datapath:<str>: file path 
    5. yyyy    :<int>: year 
    6. mm      :<int>: month 
    7. dd      :<int>: day 
    
    8. fopen   :<func>: construct class file_desc_cl for accessing orbit file  (see Note 5)
    9. fread   :<func>: read orbit data into fdesc (see Note 6)
    10. fclose :<func>: orbit file close function  (see Noye 7)
    11. fget   :<func>:  get origin orbit table
    
    12. var_name_dict:<dict>: translation between column name in the orbit file 
    --- and variable name in class 
    
    13. fio_keywords:<dict>: extra parameter for opening and reading orbit files 
   
    14. fdesc:<file_desc_cl>: file access class
    
    15. orb:<dict>:  orbit information 
    

    Functions
    ===============================================================
    1. __init__: initialization
    2. load_orbit_from_file: read orbit from file 
    3. grid_orbit: grid obit data into grid boxes defined by rlon, rlat
    4. get_var: get variable stored in self.fdesc
    Notes:
    ================================================================
    
    1) After orbit file is read, self.orb and self.var_name_dict must  
    contain the following entries:
    --- 1) 'lon' for longitude
    --- 2) 'lat' for latitude
    --- 3) 'sza' for solar zenith angle 
    --- 4) 'time' for time in seconds

    2) self.fesc.dtable contains all data read from orbit files.  
    3) self.orb contains only those with their names in var_name_dict. 
    4) self.orb could have less records than fdesc after being regridded or being filtered 
    5) fopen is expected to have an interface like:
    ---fdesc=fopen(flnm, datapath, \
       viewtype, viewmode,\
       yyyy, doy,\
       **keywords)
   
    6) fread is expected to have an interface like:
    ---fdesc=fread(fdesc, **keywords)
    
    
    7) fclose is expected to have an interface like:
    ---fclose(fdesc)
    
    
    
     
  

    
    """
    
    def __init__(self, viewtype, \
                     viewmode, \
                     datapath,\
                     flnm,\
                     yyyy, \
                     mm, \
                     dd,\
                     fopen=orb_io.open_orbit_file,\
                     fread=orb_io.read_orbit_file,\
                     fclose=orb_io.close_orbit_file,\
                     fget=orb_io.get_orbit_data,\
                     varname_dict=orb_io.orb_varname_dict,\
                     mask_val=oob.fill_val,\
                     fio_keywords={} ):
        
        """
        Initialization
        
        Inputs:
        ---------------------------------------------
        1. viewtype:<str>: instrument type
        2. viewmode:<str>: view mode such as 'nadir', 'glint' etc
        3. datapath:<str>: file path 
        4. flnm    :<str>: Name of orbit file (reserved for future use)
        5. yyyy    :<int>: year 
        6. mm      :<int>: month 
        7. dd      :<int>: day 

        8. fopen   :<func>: orbit file open function (Note 5). 
        9. fread   :<func>: orbit file read function (Note 6).
        10. fclose :<func>: orbit file close function (Note 7).
        11. varname_dict:<dict>: translation between column name in orbit file 
        --- and variable name in class (Note 1). 

        12. fio_keywords:<dict>: extra parameter for opening and reading orbit files 
        

        
        """
        
        # S1 initialize members 

        ## instrument 
        self.viewtype=viewtype
        self.viewmode=viewmode
        
        ## time 
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        ## orbit file 
        self.flnm=flnm
        self.datapath=datapath

        self.fopen=fopen
        self.fread=fread
        self.fclose=fclose
        self.fget=fget
        self.fio_keywords={}
        self.fdesc=None # file access class
        
        
        self.varname_dict=varname_dict
        
        ## data storage  
        self.orb=None  # variables needed to be kept
        self.mask_val=mask_val
        

                
        # S2 read orbit file 
        if (yyyy<>None and yyyy>0):
            self.load_orbit_from_file(yyyy, mm, dd, \
                                          **self.fio_keywords)
        
    
    
    def load_orbit_from_file(self, yyyy, mm, dd, \
                                 **keywords):
        
        """ read orbit from file 

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
        
        
        self.fdesc=self.fopen(self.flnm, self.datapath, \
                                  self.viewtype, \
                                  self.viewmode,\
                                  yyyy, mm, dd,\
                                  **keywords)
        
        # S2 read into orbit data to fdesc.dtable 
        
        self.fdsec=self.fread(self.fdesc, **keywords)
        
        # S3 fill orb will required variables 
        
        self.orb={}
        for varname in self.varname_dict:
            varname_fl=self.varname_dict[varname]
            var=self.fdesc.dtable[varname_fl]
            
            
            if ('lon' in varname):
                var=npy.where(var>180, var-360.0, var)
            elif('Lon' in varname):
                var=npy.where(var>180, var-360.0, var)
                
            self.orb.update({varname:var})
            
        
        return len(self.fdesc.dtable)

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


    def grid_orbit(self, rlon, rlat, update_orbit=True):
        
        """
        grid orbit data into grid boxes defined by rlon, rlat 
        
        Inputs:
        --------------------------------------------------
        1. rlon: <array, (nlon)>: longitude grid in ascending order. 
        2. rlat: <array, (nlat)>: latitude grid  in descending order.

        3. reset_orbit:<T/F>: if true, self.orb will be set to new gridded data
        
        
        Returns:
        -------------------------------------
        1. all_fields:<dict>: orbit data summaried at grid boxes defined by (rlonxrlat)
        
        
        """
        
        
        # S1 retrieve orbit locations 
        olon=self.orb['lon']
        olat=self.orb['lat']
        
        
        # fill the regrid obs
        
        ## setup a grid for 
        all_fields={}
        ocnt=self.orb['cnt']
        for varname in self.orb:
            odata=self.orb[varname]
            if (varname<>'cnt'):
                odata=ocnt*odata
            
            all_fields.update({varname:odata})
            
            grd_counts, grd_sum, grd_sum_2=\
                vpf.grid_profile_2d(olon,olat,odata,rlon,\
                                        rlat,0, self.mask_val)
            
            grd_sum=npy.squeeze(grd_sum)
            all_fields.update({varname:grd_sum})
        
        
        # S4 change back to single array 
    
        if (update_orbit):
          ## get the count of observation 
            grd_cnt=all_fields['cnt']
            usd_idx=npy.where(grd_cnt>0)
            # usd_idx=npy.squeeze(usd_idx)
            grd_cnt=grd_cnt[usd_idx]
            print npy.shape(grd_cnt)
            
            for varname in all_fields:
                grd_sum=all_fields[varname]
                grd_sum=grd_sum[usd_idx]
                if (varname=='cnt'):
                    self.orb.update({varname:grd_sum})
                else:
                    grd_sum=grd_sum/grd_cnt
                    self.orb.update({varname:grd_sum})
            
                print 'orbit_grid: (name), size', varname, npy.size(grd_sum)
            
        
    
        return all_fields


    
    
        
if (__name__=='__main__'):

    datapath='/home/lfeng/local_disk/otool_data/oco/aqua_0.25x0.25/'
    viewtype='aqua'
    viewmode='nadir'
    mm=1
    yyyy=2006
    dd=5
    flnm=''
    
    cl_orb=orbit_cl(viewtype, \
                        viewmode, \
                        datapath,\
                        flnm,\
                        yyyy, \
                        mm, \
                        dd,\
                        fopen=orb_io.open_orbit_file,\
                        fread=orb_io.read_orbit_file,\
                        fclose=orb_io.close_orbit_file,\
                        varname_dict=orb_io.orb_varname_dict,\
                        mask_val=oob.fill_val,\
                        fio_keywords={} )
    
    print cl_orb.orb.keys()
    olon=cl_orb.orb['lon']
    olat=cl_orb.orb['lat']
    time=cl_orb.orb['time']
    print npy.size(olon)
    print olon[0:5]
    print olat[0:5]
    print npy.max(olon), npy.min(olon)
    
    
    
    
    rlon=npy.arange(-180, 180, 60.0)
    rlat=npy.arange(-90, 90, 30.0)

    fields=cl_orb.grid_orbit(rlon, rlat)
    print fields.keys()
    
    
    olon=cl_orb.orb['lon']
    olat=cl_orb.orb['lat']
    time=cl_orb.orb['time']
    print npy.size(olon)
    print olon[0:5]
    print olat[0:5]
    

    
    

            
                  
