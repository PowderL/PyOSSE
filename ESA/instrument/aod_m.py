"""
  class for AOD PDF 
  
  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.06.28
  History: v0.95, 2013.02.24
  
  Classes:
  ==========================================
  1. aod_cl: class for AOD PDF data
  
  
"""

import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.atmosphere.compute_gc_grid as cgrd
import ESA.util.gp_grid_m as grid
import ESA.util.gp_axis_m as axis_m
import numpy.random as rnd


import aod_file_m as aod_io

import numpy as npy


class aod_cl:
    
    """
    class for aod data 

    
    Members
    ======================================
    
    1. datapath:<str>: file path for aod data 
    2. flnm: <flnm>: file name (for future use)
    
    3. yyyy: <int>:  year 
    4. mm:<int>:     month 
    5. dd:<int>:     day 
    6. fio_keywords:<dict>: key words for reading aod data
    7. varname_dict:<dict>: translations between names in the class 
    8. and column names in the file. 
    9. fdesc:<grdfile_desc_cl>: file access to aod data
    10. fopen:<func>: open data file 
    11. fread:<func>: read data in 
    12. fclose:<func>: close file 
    13. fget: <func>:  get value at given coordinates
    14. attr_dict:<dict>: attributes
    
    


    Functions
    =======================================================
    1. __init___: intialization
    2. load_data: load data from file
    3. get_data : get aod at given locations
    
    
    """
    
    
    
    def __init__(self, datapath, flnm, \
                     yyyy, mm, dd, \
                     varname_dict=aod_io.aod_varname_dict,\
                     fopen=aod_io.open_aod_file, \
                     fread=aod_io.read_aod_file,    \
                     fclose=aod_io.close_aod_file, \
                     fget=aod_io.get_aod_data,\
                     fio_keywords={},\
                     **keywords):
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for aod data 
        2. flnm: <flnm>: file name for aod file
        3. yyyy: <int>:  year 
        4. mm:<int>:     month 
        5. dd:<int>:     day 
        6. varname_dict:<dict>: tranlation between var names 
        7. fopen:<func>: open data file 
        8. fread:<func>: read data into 
        9. fclose:<func>: close file 
        10. fget: <func>:  get aod at given locations
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
        self.attr_dict.update({'ot_type':oob.ot_aod})
        
        
        self.attr_dict.update(keywords)
        
        
    def get_data(self, olon, olat, oyyyy, omm, odd, osec,**keywords):
        """
        get aod coverage at given location 
        
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
        1. aod: <array, (nob, naod)>: AOD values 
        2. paod: <array, (nob, naod)>:  Probability for AOD values
        
        
        
        Notes:
        --------------------------------------------
        1. keywords passed to self.fget will be = self.fio_keywords+keywords
        2. time (yyyy, mm, dd, sec) can be single values while olon and olat 
        two arrays of same sizes. 
        
        """

        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        
        
        aod, paod=self.fget(self.fdesc,  olon, olat, \
                           oyyyy, omm, odd, osec, \
                           **all_keywords)
        
        
        return aod, paod
    
    
        
        
    
    def load_data(self, datapath, flnm, yyyy, mm, dd,  **keywords):
    
        """
        read data from file to array/table 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for aod data 
        2. flnm: <flnm>: file name for aod file
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
    
def aod_sample(olon, olat, otime, ocnt, aod, paod, aod_limit, mask_val=oob.fill_val):
    """
    sampling AOD 

    Inputs:
    ----------------------------------------
    1. olon:<array, (nobs,)>: observation longitude 
    2. olat:<array, (nobs,)>: observation latitude
    3. otime:<array, (nobs,)>: observation time
    4. ocnt:<array, (nobs,)>: observation number 
    5. aod:<array, (nobs,naod)>: AOD values
    6. paod:<array, (nobs,naod)>: Probability for AOD values.
    7. aod_limit:<float>: upper limit 
    8. mask_val:<float>: masking value 
    
    
    Returns:
    --------------------------------------------
    1. mean_aod:<array, (nobs,)>: mean AOD  values for measurements of each observation (see Note 1).
    2. mean_cnt:<array, (nobs,)>: counts for good measurements of each observation.
    
    Notes:
    -----------------------------------------------
    1. mean_aod will be calculated only for those measurement below aod_limit, 
    if there are some good (small) AOD values sampled 
    
    
    """
    
    # S1: check observation number, and set mean_aod

    
    
    
    nobs, nval=npy.shape(aod)
    
    mean_aod=npy.zeros(nobs, float)
    mean_aod[:]=mask_val
    mean_cnt=npy.zeros(nobs)
    
    # S2: loop over observations

    sum_pb=npy.sum(paod, axis=1)
    pb_acc=npy.add.accumulate(paod, axis=1)
        
            
    for iobs in range(nobs):
        if ((sum_pb[iobs]<0.9)):
            pass
        
        elif (ocnt[iobs]>0):
            
            # #T1: sample AOD probability 
            
            pb_prof=paod[iobs,:] # probability 
            aod_prof=aod[iobs, :]  # value
            
            # ##c: normalizing  
            
            
            
            # ##c: accumulating probability 
            # ##c: as a result, pb_acc[i] stands probability for  
            # ##c:  AOD being within {aod_prof[0],  aod_prof[i]}
            # ##c: in other words, pb_acc[i]-pb_acc[i-1] is probability of AOD within aod_prof[i] and aod_prof[i-1]
            
            
            
            # ##c: generating random (0 to 1) number for measurement counts
            
            rnd_num=rnd.random(ocnt[iobs])
            
            # ##c: sampling pb_acc intervals: pb_acc[i]-pb_acc[i-1]
            
            
            
            idx=npy.searchsorted(pb_acc[iobs,:], rnd_num) # probability  
            idx=npy.where(idx>=nval, nval-1, idx)
            
            rnd_aod=aod_prof[idx] # aod value for each individual measurment
            rnd_pb=pb_prof[idx]  #  probabilty 
            
            # S3: find mean for good (aod<limit) measurements if available
            
            usd_idx=npy.where(rnd_aod<=aod_limit)
            good_pb=npy.sum(rnd_pb[usd_idx])
            
            if (good_pb>0):
                # #T1: good measurement
                
                obs_aod=npy.mean(rnd_aod[usd_idx])
                obs_cnt=npy.size(usd_idx)
                
            else:
                # #T2: no good measurement
                
                obs_aod=npy.mean(rnd_aod)
                obs_cnt=0
                
            
                
            mean_aod[iobs]=obs_aod
            mean_cnt[iobs]=obs_cnt


    return mean_aod, mean_cnt

            
          
     
    
    
    


if (__name__=='__main__'):
     
    datapath='/home/lfeng/local_disk/otool_data/clim_dat/'
    dd=2
    mm=1
    yyyy=2009
    flnm=''
    
    aod_m=aod_cl(datapath, flnm, \
                       yyyy, mm, dd)
    
    

    
    olon=npy.arange(-180, 180, 60.0)
    olat=npy.arange(-90, 90, 30.0)

    oyyyy=2009
    omm=1

    odd=20
    
    osec=3600.0
    
    
    aod, paod=aod_m.get_data(olon, olat, \
                                 oyyyy, omm, odd, osec)
    
    print 'aod axis for 200901:', aod[2,:]
    
    
    print 'paod for 200901:', paod[2,:]
    
    
        



