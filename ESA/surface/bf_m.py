""" class for  basis functions as bridge to access flux file 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.01.21
    
    
    Classes:
    ===============================================
    1. bf_cl: class for basis functions 
    
    """
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_var_io as ovio
import ESA.util.otool_ncfile_io as nfio

import ESA.util.time_module as tm
import ESA.util.process_nf_array as pnf

import ESA.util.gp_grid_m as grid
import ESA.util.gp_axis_m as axis_m
import bf_file_m as bf_ncio

class bf_cl:
    """
    Classes for basis functions

    Members:
    ----------------------------------------
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


    
    Functions:
    --------------------------------------------
    1. __init__: initialization 
    2. get_bf: get basis function 
    3. load_bf: read data from file to array
    
    
    """
    
    def __init__(self, datapath, flnm, yyyy, mm, dd,  \
                     varname_lst=bf_ncio.bf_varname_lst,\
                     varname_dict={},\
                     fopen=bf_ncio.open_bf_file, \
                     fread=bf_ncio.read_bf_file, \
                     fget=bf_ncio.get_bf, \
                     fclose=bf_ncio.close_bf_file, \
                     fio_keywords={}, **keywords):
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for aod data 
        2. flnm: <flnm>: file name for aod file
        3. yyyy: <int>:  year 
        4. mm:<]int>:     month 
        5. dd:<int>:     day 
        6. varname_list:<list, t:str>: list of variables to be read
        7. varname_dict:<dict>: tranlation between var names 

        8. fopen:<func>: open data file 
        9. fread:<func>: read data into 
        10. fclose:<func>: close file 
        11. fio_keywords: <dict>:  extra inputs for file accessing
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
        self.varname_lst=varname_lst
        self.varname_dict=varname_dict
        self.fio_keywords=fio_keywords
        
        
        # clound file access function
        
        self.fopen=fopen
        self.fread=fread
        self.fclose=fclose
        self.fget=fget
        
        # load data 

        self.fdesc=None
        
        self.is_ready=self.load_bf(self.datapath, self.flnm,\
                                       self.yyyy, \
                                       self.mm, self.dd \
                                       )
        
        # setup attribute  dictionary 
        
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_bf})
        
        
        self.attr_dict.update(keywords)
        
        
    def get_bf(self, ilayer,**keywords):
        """
        get basis function at given layer
        
        Inputs:
        ---------------------------------------------
        1.ilayer:<integer>: layer index

        
        2. keywords:<dict>: extra inputs
        ---reserved keys:
        ---varname_lst:<list, t:str>: name of variables.
        ---flux_nm:<str>:name of flux
        ---map_nm:<str>: name of map
        ---area_nm:<str>: name of area
        ---lon_nm:<str>: name of lon
        ---lat_nm:<str>: name of lat
        
        Returns:
        -------------------------------------------------------------------
        1. var_lst:<list>: list of variable read (see Note 1)
        
        Notes:
        ---------------------------------------
        1. if varname_lst is not set in keywords. A default list of variable will be returned 
        This list includes:
        ---1. lon:<array, (nlon)>: longitude
        ---2. lat:<array, (nlat)>: latitude
        ---3. sel_flux:<array, (nlon, nlat)>: basis functions 
        ---4. sel_map:<array, (nlon, nlat)>: map define the region  
        ---5. sel_area:<array, (nlon, nlat)>: areas for the region 
    
        """

        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        
        
        var_lst=self.fget(self.fdesc,  ilayer, \
                              **all_keywords)
        
        
        return var_lst
    
    
    
    
    def load_bf(self, datapath, flnm, yyyy, mm, dd,  **keywords):
    
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
        
        flnm=flnm.strip()
        
        if (len(flnm)==0):
            flnm=self.flnm
        else:
            self.flnm=flnm
        
        
        datapath=datapath.strip()
        if (len(datapath)==0):
            datapath=self.datapath
        else:
            self.datapath=datapath
            
        
        all_keywords={}
        all_keywords.update(self.fio_keywords)
        all_keywords.update(keywords)
        print 'in load bf', self.varname_lst
        print flnm
        
        self.fdesc=self.fopen(datapath, flnm, \
                                  yyyy, mm, dd,\
                                  varname_lst=self.varname_lst,\
                                  var_ncname_dict=self.varname_dict,\
                                  read_to_array=False,\
                                  fio_keywords=all_keywords)
        
        self.fdesc=self.fread(self.fdesc, **all_keywords)
        
        return True



    
#<<< TESTS >>>

if (__name__=='__main__'):
    print 'See gen_ensemble_flux.py'
    

    
