"""
   Class for cloud PDF data
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.26
   
   Classes:
   ============================================
   1. cloud_cl: class for cloud PDF data
   

"""
import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.atmosphere.compute_gc_grid as cgrd
import ESA.util.gp_grid_m as grid
import ESA.util.gp_axis_m as axis_m


import cloud_file_m as cld_io
import numpy.random as rnd



class cloud_cl:
    
    """
    class for cloud data 

    
    Members
    ======================================
    
    1. datapath:<str>: file path for cloud data 
    2. flnm: <flnm>: file name (for future use)
    
    3. yyyy: <int>:  year 
    4. mm:<int>:     month 
    5. dd:<int>:     day 
    6. fio_keywords:<dict>: key words for reading cloud data
    7. varname_dict:<dict>: translations between names in the class 
    8. and column names in the file. 
    9. fdesc:<grdfile_desc_cl>: file access to cloud data
    10. fopen:<func>: open data file 
    11. fread:<func>: read data in 
    12. fclose:<func>: close file 
    13. fget: <func>:  get value at given coordinates
    14. attr_dict:<dict>: attributes

    


    Functions
    =======================================================
    1. __init___: intialization
    2. load_data: load data from file
    3. get_data : get cloud at given locations
    
    
    """
    
    
    
    def __init__(self, datapath, flnm, \
                     yyyy, mm, dd, \
                     varname_dict=cld_io.cld_varname_dict,\
                     fopen=cld_io.open_cloud_file, \
                     fread=cld_io.read_cloud_file,    \
                     fclose=cld_io.close_cloud_file, \
                     fget=cld_io.get_cloud_data,\
                     fio_keywords={},\
                     **keywords):
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for cloud data 
        2. flnm: <flnm>: file name for cloud file
        3. yyyy: <int>:  year 
        4. mm:<int>:     month 
        5. dd:<int>:     day 
        6. read_to_array:<T/F>: if True, data will be read from file as arrays. 
        7. use_current_time:<T/F>: 
        
        ---if True, input time (yyyy, mm, dd) will be used to read data.
        ---if False, default time will be used to read data 
        
       
        8. fopen:<func>: open data file 
        9. fread:<func>: read data into 
        10. fclose:<func>: close file 
        11. fget: <func>:  get cloud at given locations
        12. keywords:<dict>: attribute 
        
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
        self.attr_dict.update({'ot_type':oob.ot_cloud})
        
        
        self.attr_dict.update(keywords)
        
        
    def get_data(self, olon, olat, oyyyy, omm, odd, osec,**keywords):
        """
        get cloud coverage at given location 
        
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
        1. pcld: <array, (nob)>: cloud coverage percentage
        
        Notes:
        --------------------------------------------
        1. keywords passed to self.fget will be self.fio_keywords+keywords
        2. time (yyyy, mm, dd, sec) can be single values while olon and olat 
        two arrays of same sizes. 
        
        """

        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        
        
        pcld=self.fget(self.fdesc,  olon, olat, \
                           oyyyy, omm, odd, osec, \
                           **all_keywords)
        
        
        return pcld
    
    
    
        
    
    def load_data(self, datapath, flnm, yyyy, mm, dd,  **keywords):
    
        """
        read data from file to array/table 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for cloud data 
        2. flnm: <flnm>: file name for cloud file
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
    
    

def cloud_sample(olon, olat, otime, ocnt, cloud_pdf, fpenalty=None, **keywords):
    """
    
    Set cloud flag according cloud PDF 
    
    Inputs:
    ------------------------------------------
    1. olon:<array, (nobs,)>: observation longitude 
    2. olat:<array, (nobs,)>: observation latitude
    3. otime:<array, (nobs,)>: observation time
    4. ocnt:<array, (nobs,)>: observation number 
   
    5. cloud_pdf:<array>: cloud percentage
    6. fpenalty:<funct>:cloud  PDF penality function.   
    7. keywords:<dict>: extra contral variables
    
    Returns:
    ----------------------------------------------------------
    1. clear_cnt:<array>: the count number of cloud free scene
    
    """
     
    
    if (oob.get_ot_type(ocnt)==oob.ot_array):
        # if given a list 
        nobs=npy.size(ocnt)
        clear_cnt=npy.zeros(nobs)
        for iobs in range(nobs):
            clear_cnt[iobs]=ocnt[iobs]
            cld_pb=cloud_pdf[iobs]
            if (ocnt[iobs]>0):
                # ## add size penality etc 
                if (fpenalty<>None):
                    cld_pb=fpenalty(olon, olat, otime, cld_pb, **keywords)
                 
                clear_pb=1.0-cld_pb
                cloud_free_s=rnd.random(ocnt[iobs])
                usd_idx=npy.where(cloud_free_s<=clear_pb)
                usd_idx=npy.squeeze(usd_idx)
                clear_cnt[iobs]=npy.size(usd_idx)

        #  loop iobs
        
        return clear_cnt
    else:
        cld_pb=cloud_pdf

        if (fpenalty<>None):
            cld_pb=fpenalty(olon, olat, otime, cld_pb, **keywords)
        clear_pb=1.0-cld_pb
        cloud_free_s=rnd.random()
        if (cloud_free_s<clear_pb):
            clear_cnt=1
        else:
            clear_cnt=0
            
        return clear_cnt

         
     
             
                               
    



if (__name__=='__main__'):
     
    datapath='/home/lfeng/local_disk/otool_data/clim_dat/'
    viewtype='aqua'
    viewmode='nadir'
    dd=1
    mm=1
    yyyy=2009
    flnm=''
    
    cld_m=cloud_cl(datapath, flnm, \
                       yyyy, mm, dd)
    
    
    
    
    olon=npy.arange(-180, 180, 60.0)
    olat=npy.arange(-90, 90, 30.0)

    oyyyy=2006
    omm=3
    odd=20
    
    osec=3600.0
    
    
    cld=cld_m.get_data(olon, olat, \
                           oyyyy, omm, odd, osec)
    
    print 'cld for 200603:', cld
    
        
        



