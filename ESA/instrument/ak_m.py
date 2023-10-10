"""
  class for scene-dependent  averaging kernel   
  
  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.06.28
  History: v0.95, 2013.02.24
  
  Classes:
  ==========================================
  1. ak_cl: class for ak data
  
  """

import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.otool_obj as  oob
import ESA.util.geo_constant as gc
import ESA.util.gp_grid_m as grid

import ak_file_m as ak_io



class ak_cl:
    
    """
    class for ak data 
    Members:
    -------------------------------------------
    1. flnm:<str>: file name (reserved for future use)
    2 .datapath:<str>: file path 
    
    # current time
    3. yyyy=yyyy
    4. mm=mm
    5.dd=dd
    
    # observation
    
    6. viewmode_dict:<dict>: dictionary for viewmode and its ID
    7. surf_type_dict:<dict>: dictionary for surface type and its ID
    
    # AK file access function
        
    8.  fopen:<func>: open file 
    9.  fread;<func>: read file 
    10. fclose:<func>: close file 
    11. fget:<func>: read data in 
         
    
    """
    
    
    
    def __init__(self, viewtype, \
                     datapath, flnm, \
                     yyyy, mm, dd, \
                     viewmode_dict=ak_io.def_ak_viewmode_dict,\
                     surf_type_dict=ak_io.def_ak_surf_type_dict,\
                     fopen=ak_io.open_ak_file, \
                     fread=ak_io.read_ak_file,    \
                     fclose=ak_io.close_ak_file, \
                     fget=ak_io.get_ak_data,\
                     fio_keywords={},\
                     **keywords):
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. viewtype:<str>: instrument type (or name)
        2. datapath:<str>: file path for AK data 
        3. flnm: <flnm>: file name for AK data
        4. yyyy: <int>:  year 
        5. mm:<int>:     month 
        6. dd:<int>:     day 
        7. viewmode_dict:<dict>: dictionary for view modes
        8. surf_type_dict:<dict>: dictionary for surface type

        9. fopen:<func>: open data file 
        10. fread:<func>: read data into 
        11. fclose:<func>: close file 
        12. fget: <func>: get ak values
        
        13. rw_keywords:<dict>: additional parameters for AK file read/write
        14. keywords:<dict>: attributes
        
        """
        
        # file 
        
        self.viewtype=viewtype 
        self.flnm=flnm
        self.datapath=datapath
        
        # current time
       
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        
        # dictionary on view mode and surface type
        
        
        self.viewmode_dict=viewmode_dict
        self.surf_type_dict=surf_type_dict
        self.fio_keywords=fio_keywords
        

        
        # clound file access function
        
        self.fopen=fopen
        self.fread=fread
        self.fclose=fclose
        self.fget=fget
        
        # load data 
        
        self.viewtype=viewtype
        viewmodes=self.viewmode_dict.keys()
        surf_types=self.surf_type_dict.keys()
        
        viewmode=viewmodes[0]
        surface=surf_types[0]
        
        # construct fdesc
        self.fdesc=None
        
        
        self.is_ready=self.load_data(flnm, datapath, \
                                         self.viewtype,\
                                         viewmode,\
                                         surface,\
                                         yyyy, mm,dd)
        
        

        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_ak})
                                      
        self.attr_dict.update(keywords)
        
        
        # setup the dictionary 
        
        
        
    def get_viewmode_id(self, viewmode, viewmode_dict=None):
        """ find ID for viewmode 
        Inputs:
        -----------------------------------------
        1.viewmode:<array/obj>: viewmodes to be checked against viewmode_dict
        2.viewmode_dict:<dict>: dictionary in the form of {viewmode:id}
        ---if viewmode_dict==None, self.viewmode will be used
        
        
        Returns:
        ------------------------------------------
        vid_lst:<array/obj>: id of the given viewmodes

        """
        
        
        if (viewmode_dict==None):
            viewmode_dict=self.viewmode_dict
        
        vid_lst=oob.translate_dict_ID(viewmode_dict, viewmode)
        
        return vid_lst
    

    
    def get_surf_type_id(self, surf_type, surf_type_dict=None):
        
        """ find ID for surface type
        
        
        Inputs:
        -----------------------------------------
        1.viewmode:<array/obj>: viewmodes to be checked against viewmode_dict
        2.viewmode_dict:<dict>: dictionary in the form of {viewmode:id}
        ---if surf_type_dict==None, self.surf_type_dict will be used
        
        Returns:
        ------------------------------------------
        vid_lst:<array/obj>: id of the given surface types
        
        """
        if (surf_type_dict==None):
            surf_type_dict=self.surf_type_dict
        
        vid_lst=oob.translate_dict_ID(surf_type_dict, surf_type)
        
        

        return vid_lst
    
            
            
    def get_data(self, osza, oaod, iviewmode, \
                     isurf_type, **keywords):
        
        """
        get ak for observation conditions (osza and oaod)
        
        Inputs:
        -------------------------------------------------------
        1.osza:<array, (nob)>: solar zenith angle 
        2.oaod:<array, (nob)>: aerosol optical depth 
        3.iviewmode:<array, (nob)/integer>:  view mode
        4.isurf_type:<array, (nob)/integer>: surface type 
        5. keywords:<dict>: 
        ---Reserved Entries
        ---common_ref:<str/numeric>: the reference shared by observations. 
        --->'ak_viewmode_dict':<dict>: dictionary for view mode
        --->'ak_surf_type_dict': <dict>: dictionary for surface type
        ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
        ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
        ---> zval: <array, (nz,)>: vertical axis 
        --->'replaces': <dict>: words to be replaced when decoding line 
    
        
                
        
        Returns:
        --------------------------------------------------
        1. lgp:<array, (nob, nz)>: log10(pres) values
        2. ak: <array,  (nob, nz)>: averaging kernels
        
        """
        
        all_keywords=dict(self.fio_keywords)
        all_keywords.update(keywords)
        
        lgp, ak=self.fget(self.fdesc, osza, oaod, \
                              iviewmode, isurf_type, \
                              **all_keywords)
        
        return lgp, ak
    
    
    
    def load_data(self, flnm, datapath, \
                      viewtype,\
                      viewmode,\
                      surface,\
                      yyyy, mm,dd,\
                      **keywords):
        """
        
        Construct fdesc and read AK file in 
        
        
        1. flnm:<str>:   name of file (reserved for future use )
        2. datapath:<str>: file path
        3. viewtype:<str>: instrument type 
        4. viewmode:<str>: nadir or glint view
        5. surface:<str>: surface type
        6. yyyy: <int>: year
        7. mm: <int>: month
        8. dd:<int>: day 
        9. keywords:<dict>: extra inputs for file reading 
        ----Reserved words
        --->'fl_colnm_lst':<list, t:str>: names of all columns in the
        --->'ak_viewmode_dict':<dict>: dictionary for view mode
        --->'ak_surf_type_dict': <dict>: dictionary for surface type
        --->'delim':<str>: splitor in averaging kernel (text) file
         ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
         ---> zval: <array, (nz,)>: vertical axis 
         --->'replaces': <dict>: words to be replaced when decoding line 
         
         Returns:
         -----------------------------------------
         1. is_ready:<T/F>: True when the file has been opened successfully 


         Notes:
         1. self.fdesc will be constructed by functions self.fopen. 

         2. Its data memebers will be updated by self.fread. 

        """
        all_keywords=dict(self.fio_keywords)
        all_keywords.update(keywords)
        
        if (self.fdesc<>None):
            self.fclose(self.fdesc)
            
        self.fdesc=self.fopen(flnm, datapath, \
                                  viewtype, \
                                  viewmode, surface, \
                                  yyyy, mm, dd,\
                                  **all_keywords)
        
        self.fdesc=self.fread(self.fdesc, **all_keywords)
        
        return True
    
    
    

# <<< TEST>>>

if (__name__=='__main__'):
    

    datapath='/home/lfeng/local_disk/otool_data/oco/oco_ak/'
    viewtype='aqua'
    viewmode='nadir'
    surface='snow'
    
    mm=3
    dd=2
    yyyy=2006
    flnm=''
    
    
    ak_m=ak_cl(viewtype, \
                   datapath, flnm, \
                   yyyy, mm, dd\
                   )
    
    
    osza=npy.array([10.0,15.0])
    oaod=npy.array([0.2, 0.10])
    iviewmode=[1,0]
    isurf_type=[1,0]
    
    lgp, ak=ak_m.get_data(osza, oaod, \
                              iviewmode, \
                              isurf_type)
    
    print 'shape(lgp):',npy.shape(lgp)
    print 'shape(ak):',npy.shape(ak)
    
    print 'lgp[0,:]:', lgp[0,:]
    print 'ak[0,:]',ak[0,:]
    
    








        
        



