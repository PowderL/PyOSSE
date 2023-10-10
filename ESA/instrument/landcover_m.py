"""
  class for LC PDF 
  
  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.06.28
  History: v0.95, 2013.02.24
  
  Classes:
  ==========================================
  1. lc_cl: class for LC PDF data
  
  
  Functions:
  =========================================
  1. sample_lc: sample lc PDF, and return PDF
  
"""

import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.atmosphere.compute_gc_grid as cgrd
import ESA.util.gp_grid_m as grid
import ESA.util.gp_axis_m as axis_m
import landcover_file_m as lc_io

import numpy as npy

#===<CLASSES>===

class lc_cl:
    
    """
    class for lc data 

    
    Members
    ======================================
    
    1. datapath:<str>: file path for lc data 
    2. flnm: <flnm>: file name (for future use)
    
    3. yyyy: <int>:  year 
    4. mm:<int>:     month 
    5. dd:<int>:     day 
    6. fio_keywords:<dict>: key words for reading lc data
    7. varname_dict:<dict>: translations between variablr names  
       and names in the file.
    8. fdesc:<grdfile_desc_cl>: file access to lc data
    9. fopen:<func>: open data file 
    10. fread:<func>: read data in 
    11. fclose:<func>: close file 
    12. fget: <func>:  get value at given coordinates
    13. attr_dict:<dict>: attributes
    
    


    Functions
    =======================================================
    1. __init___: intialization
    2. load_data: load data from file
    3. get_data : get lc at given locations
    
    
    """
    
    
    
    def __init__(self, datapath, flnm, \
                     yyyy, mm, dd, \
                     varname_dict=lc_io.lc_varname_dict,\
                     fopen=lc_io.open_lc_file, \
                     fread=lc_io.read_lc_file,    \
                     fclose=lc_io.close_lc_file, \
                     fget=lc_io.get_lc_data,\
                     fio_keywords={},\
                     **keywords):
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for lc data 
        2. flnm: <flnm>: file name for lc file
        3. yyyy: <int>:  year 
        4. mm:<int>:     month 
        5. dd:<int>:     day 
        6. varname_dict:<dict>: tranlation between var names 
        
        7. fopen:<func>: open data file 
        8. fread:<func>: read data into 
        9. fclose:<func>: close file 
        10. fget: <func>:  get lc at given locations
        11. keywords:<dict>: attribute 
        
        """
        
        # file 
        
        self.flnm=flnm
        self.datapath=datapath
        
        # current time
       
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        self.is_ready=False
        
        self.varname_dict=varname_dict
        
        self.fio_keywords=fio_keywords
        
        
        # clound file access function
        
        self.fopen=fopen
        self.fread=fread
        self.fclose=fclose
        self.fget=fget
        
        # load data 

        self.fdesc=None
        
        self.is_ready=self.load_data(self.datapath, \
                                         self.flnm, self.yyyy, \
                                         self.mm, self.dd \
                                         )
        
        # setup attribute  dictionary 
        
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_lc})
        
        
        self.attr_dict.update(keywords)
        
        
    def get_data(self, olon, olat, oyyyy, omm, odd,osec, **keywords):
        """
        get lc coverage at given location 
        
        Inputs:
        ---------------------------------------------
        1.olon:<array, (nob,)>: longitude
        2.olat:<array, (nob,)>: latitude  
        3.oyyyy:<int/array, (nob,)>:  year  (Note 2)
        4.omm:<int/array, (nob,)>:  month   (Note 2)
        5.odd:<int/array, (nob,)>:  day      (Note 2)
        6.osec:<int/array, (nob,)>:  seconds in the day  (Note 2).
        
        7. keywords:<dict>: extra inputs
        ---Reserved Entries
        -->common_ref:<str/numeric>: the reference shared by observations.
        
        
        
        Returns:
        ------------------------------------------
        1. lc: <array, (nob, nlc)>: land cover  
        
        
        Notes:
        --------------------------------------------
        1. keywords passed to self.fget will be self.fio_keywords+keywords
        2. time (yyyy, mm, dd, sec) can be single values while olon and olat 
        two arrays of same sizes. 
        
        """

        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        
        
        lc=self.fget(self.fdesc,  olon, olat, \
                         oyyyy, omm, odd, \
                         **all_keywords)
        
        return lc
        
        
    
    def load_data(self, datapath, flnm, yyyy, mm, dd,  **keywords):
    
        """
        read data from file to array/table 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for lc data 
        2. flnm: <flnm>: file name for lc file
        3. yyyy: <int>:  year 
        4. mm:<int>:     month 
        5. dd:<int>:     day 
        6. keywords:<dict>: extra inputs for reading file 
        
        
        

        Returns:
        -------------------------------------
        1. is_ready:<T/F>: True if data are read successifully. 


        Notes:
        -----------------------------------------------
        
        """
        if (self.fdesc<>None):
            self.fclose(self.fdesc)
        
        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        
        self.fdesc=self.fopen(datapath, flnm, \
                                  yyyy, mm, dd,\
                                  **all_keywords)
        
        self.fdesc=self.fread(self.fdesc, **all_keywords)
        
        return True

#===<FUNCTIONS>===

    
def sample_lc(lc, lc_st_lst):
    
    """
    Get surface type from  land coverage percentage
    
    
    Inputs:
    ------------------------------------------
    1. lc:<array, (nobs, nlc)>:   land cover
    2. lc_type_lst:<list, t:str>: landcover to surface type used in averaging kernel or observation errors
 

    Returns:
    ----------------------------------------
    1. stype_lst:<list, t:str>: list of surface type

    """
    
    dims=npy.shape(lc)
    
    if (len(dims)==1):
    
        lc_pb=lc[:]
        idx=npy.argmax(lc_pb)
        stype_lst=lc_st_lst[idx]
        return stype_lst    
    
    else:
        
        nobs, nlc=npy.shape(lc)
        stype_lst=list()
        for iobs in range(nobs):
            lc_pb=lc[iobs, :]
            idx=npy.argmax(lc_pb)
            stype_lst.append(lc_st_lst[idx])
        return stype_lst    
        
    

if (__name__=='__main__'):
     
    datapath='/home/lfeng/local_disk/otool_data/clim_dat/'
    flnm='glc2000.nc'
    dd=2
    mm=3
    yyyy=2006
    
    
    lc_m=lc_cl(datapath, flnm, \
                       yyyy, mm, dd)
    
    

    
    olon=npy.arange(-180, 180, 60.0)
    olat=npy.arange(-90, 90, 30.0)

    oyyyy=2006
    omm=3
    odd=20
    
    osec=3600.0
    
    
    lc=lc_m.get_data(olon, olat, \
                           oyyyy, omm, odd, osec)
    
    print 'lc for 200603:', lc[:,19]
    
    
    
    
        



