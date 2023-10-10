""" class  for scene-dependent observation error 
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.19
    
    Classes:
    =================================================
    1.  err_cl: class for tablized observation errors
    
"""
import  numpy as npy

import ESA.util.time_module as  tm
import ESA.util.otool_obj as  oob
import ESA.util.geo_constant as gc
import ESA.util.gp_grid_m as grid

import err_file_m as err_io



class err_cl:
    
    """
    class for err data 
    
    Members:
    ==================================================================
    1. flnm:<str>: file name (reserved for future use)
    2 .datapath:<str>: file path 
    
    # current time
    3. yyyy=yyyy
    4. mm=mm
    5.dd=dd
    
    # observation
    6. viewmode_dict:<dict>: dictionary for viewmode and its ID
    7. surf_type_dict:<dict>: dictionary for surface type and its ID
    
    # ERR file access function
        
    8.  fopen:<func>: open file 
    9.  fread;<func>: read file 
    10. fclose:<func>: close file 
    11. fget:<func>: read data in 
         

    Functions
    =======================================================================
    1. __init__: initialization.
    2. get_data: get observation errors at given observation scenes.
    3. load_data: read data from observation error files.
    4. get_viewmode_id: find IDs for viewmodes. 
    5. get_surf_type_id: find IDs for surface types. 
    
    
    """
    
    
    
    def __init__(self, viewtype, \
                     datapath, flnm, \
                     yyyy, mm, dd, \
                     viewmode_dict=err_io.def_err_viewmode_dict,\
                     surf_type_dict=err_io.def_err_surf_type_dict,\
                     fopen=err_io.open_err_file, \
                     fread=err_io.read_err_file,    \
                     fclose=err_io.close_err_file, \
                     fget=err_io.get_err_data,\
                     fio_keywords={},\
                     **keywords):
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. viewtype:<str>: instrument type (or name)
        2. datapath:<str>: file path for ERR data 
        3. flnm: <flnm>: file name for ERR data
        4. yyyy: <int>:  year 
        5. mm:<int>:     month 
        6. dd:<int>:     day 
        7. viewmode_dict:<dict>: dictionary for view modes
        8. surf_type_dict:<dict>: dictionary for surface type

        9. fopen:<func>: open data file 
        10. fread:<func>: read data into 
        11. fclose:<func>: close file 
        12. fget: <func>: get err values
        
        13. rw_keywords:<dict>: additional parameters for ERR file read/write
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
        
        print viewmode
        print surface
        
        self.is_ready=self.load_data(flnm, datapath, \
                                         self.viewtype,\
                                         viewmode,\
                                         surface,\
                                         yyyy, mm,dd, \
                                         **self.fio_keywords)
        


        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_err})
                                      
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
        get err for observation conditions (osza and oaod)
        
        Inputs:
        -------------------------------------------------------
        1.osza:<array, (nob)>: solar zenith angle 
        2.oaod:<array, (nob)>: aerosol optical depth 
        3.iviewmode:<array, (nob)/integer>:  view mode
        4.isurf_type:<array, (nob)/integer>: surface type 
        5. keywords:<dict>: 
        ---Reserved Entries
        ---common_ref:<str/numeric>: the reference shared by observations. 
        --->'err_viewmode_dict':<dict>: dictionary for view mode
        --->'err_surf_type_dict': <dict>: dictionary for surface type
        ---> ax_colnm_lst: <list, t:str>: column name for horizontal axis
        ---> zname: <str>: axis name for vertical (e.g.,  log10(pressure)) axis
        ---> zval: <array, (nz,)>: vertical axis 
        --->'replaces': <dict>: words to be replaced when decoding line 
    
        
                
        
        Returns:
        --------------------------------------------------
        1. err: <array,  (nob, nz)>: observation errors

        
        """
        
        all_keywords=dict(self.fio_keywords)
        all_keywords.update(keywords)
        
        err=self.fget(self.fdesc, osza, oaod, \
                          iviewmode, isurf_type, \
                          **all_keywords)
        
        return err
    
    
    
    
    
    
    
    
    def load_data(self, flnm, datapath, \
                      viewtype,\
                      viewmode,\
                      surface,\
                      yyyy, mm,dd,\
                      **keywords):
        """
        
        Construct fdesc and read ERR file in 
        
        
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
        --->'err_viewmode_dict':<dict>: dictionary for view mode
        --->'err_surf_type_dict': <dict>: dictionary for surface type
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
        if (self.fdesc<>None):
            self.fclose(self.fdesc)
            
        self.fdesc=self.fopen(flnm, datapath, \
                                  viewtype, \
                                  viewmode, surface, \
                                  yyyy, mm, dd,\
                                  **keywords)
        
        self.fdesc=self.fread(self.fdesc, **keywords)
        
        return True
    
    
    


if (__name__=='__main__'):
    
    
    datapath='/home/lfeng/local_disk/otool_data/oco/oco_err/'
    viewtype='aqua'
    viewmode='nadir'
    surface='snow'
    
    mm=3
    dd=2
    yyyy=2006
    flnm=''
    
    
    err_m=err_cl(viewtype, \
                   datapath, flnm, \
                   yyyy, mm, dd\
                   )
    
    
    osza=npy.array([10.0,15.0])
    oaod=npy.array([0.2, 0.10])
    iviewmode=[1,0]
    isurf_type=[1,0]
    
    err=err_m.get_data(osza, oaod, \
                           iviewmode, \
                           isurf_type)
    
    print 'shape(err):',npy.shape(err)
    
    
    print 'err[:]',err[:]
    
    








        
        



