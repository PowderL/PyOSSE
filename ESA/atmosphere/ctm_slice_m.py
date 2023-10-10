""" Class for model profiles 

    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.10.28
    
    Classes:
    --------------------------------------------
    1. ctm_slice_cl:  general container for set of profiles for 
    tracers (or geo-physical variables) at set of locations
""" 
import numpy as npy
# from scipy import *

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.vertical_interp_m as vintpl_m


import ctm_grid_3d as cgrd
import gc_grid_3d as gcgrd
import ctm_field_m as cfm



         
class ctm_slice_cl:
    
    """ general container for set of profiles for 
    tracers (or geo-physical variables) at set of locations. 
    
    Members:
    ------------------------------------------------------------------------------
    1. name: <str>: name (container name etc)
    2. lon:<array, nob>: longitude
    3. lat:<array, nob>: latitude
    4. pressure:<array, (nob, nz)>: pressure level at ob locations 
    5. is_in_log:<T/F>:  if pressure given in log space
    6. do_reverse:<T/F>: T, if pressure is given in decending order
    7. mask_val:<float>: filling value for missing and bad data 
    8. data_lst:<list, t:ctm_prof>: list of profiles
    
    9. attr_dict:<dict>: dictionary for name and ndim etc
    
    10. ndata: <integer>: length of the data set
    11. time_tag:<numerical/string>: tag for time step
    
        
    
    Functions:
    ---------------------------------------
    
    1. __init__: initialization
    
    # member retrieve 
    
    2. find_profile_cl:   search data matching selection criteria of name, tracer_id, tracer_name, \
    category, group and time range etc 
    3. __getitem__: over-riding index 

    # member manage 
    4. add_data: add data to the list 
    5. del_data: delete data from the list
    6. update_data: update data
    
    # grid configuration
    7. set_pressure: set up pressure levels 
    8. set_lon: set longitude
    9. set_lat: set latitude
    
    # attribute
    10. add_attr: add attribute 
    11. get_attr: get attribute


    """
    
    
    def __init__(self,name, data_lst, \
                     lon=None,\
                     lat=None,\
                     pressure=None,\
                     is_in_log=False,\
                     do_reverse=False,\
                     time_tag="", \
                     mask_val=oob.fill_val,\
                     extra_axis=None,\
                     **keywords):
        
        
        """Initialize 
        
        Inputs:
        ---------------------------------------------

        1. name:<str>: Name of the object
        2. data_lst<list, t:array>: list of (ctm_profile_cl ) profiles 
        3. lon: <array, (nob,)>: longitude 
        4. lat: <array, (nob,)>: latitude
        5. pressure:<array, (nob, nlvl)>: pressure
        6. is_in_log:<T/F>: Ture if pressure is given as log10(pressure) 
        7. do_reverse:<T/F>: True if pressure is given in decending order
        8. time_tag:<str/float>: time 
        9. mask_val:<float>: filling values for missing and bad data
        10. extra_axis:<gp_axis>: extra axis for time or ensemble etc
        11. keywords:<dict>: attributes 
        
        """
        
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type': oob.ot_slice})
        
        self.name=name
        self.attr_dict.update({'name':self.name})
        
        self.mask_val=mask_val
        self.attr_dict.update({'mask_val':self.mask_val})
        
        
        self.ndata=0
        self.data_lst=[]
        
        for cdata in data_lst:
            self.add_data(cdata)
        
        # set grid type
        
        self.lon=None
        self.set_lon(lon)
        
        self.lat=None
        self.set_lat(lat)
        
        self.pressure=None
        
        self.is_in_log=is_in_log
        self.do_reverse=do_reverse
        
        self.vintpl=None
        # when pres is available, vintp will be set 
        self.set_pressure(pressure, is_in_log, do_reverse, self.mask_val)
        

        self.time_tag=time_tag
        self.attr_dict.update({'time_tag':self.time_tag})
        
        
        # S4 assign extral axis (such as time, or ensemble index etc)
        
        self.extra_axis=None
        
        if (extra_axis==None):
            self.extra_axis=None
        else:
            if (oob.check_type(extra_axis, oob.ot_axis)):
                self.extra_axis=extra_axis
            else:
                msm.show_err_msg('extra axis', msm.msm_wrong_type)
                
                
        
        for keyname in keywords:
        
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
        
    
        
    def set_lon(self, olon):
        self.lon=olon
        self.attr_dict.update({'lon':self.lon})
    
    
    def set_lat(self, olat):
        self.lat=olat
        self.attr_dict.update({'lat':self.lat})
    

    def set_pressure(self, pressure, is_in_log, do_reverse, mask_val):
        
        """ define the grid pressure level

        Inputs: 
        1.pressure:<array, (nob,nz)>: pressure at each level
        2. is_in_log:<T/F>: True when pressure is given in log10
        3. do_reverse:<T/F>: True when pressure is given at decending order
        
        """
        
        
        self.is_in_log=is_in_log
        self.pressure=pressure
        self.do_reverse=do_reverse
        self.mask_val=mask_val
        
        if (oob.check_type(pressure, oob.ot_array)):
            
            self.vintpl=vintpl_m.vintpl_cl(pressure, \
                                               self.is_in_log, \
                                               self.do_reverse, \
                                               mask_val=self.mask_val)
            
        else:
            # remove the  
            self.vinterp=None
        
        
        
    
        
        
    
    def add_data(self, data):
        
        """
        Add data (ctm_field_cl) into the class 
        Inputs:
        ------------------------------
        1. data:<ctm_field_cl>: class for ctm data
        """
        
        if (oob.check_type(data, oob.ot_profile)):
            self.data_lst.append(data)
            self.ndata=self.ndata+1
            print 'len', len(self.data_lst)
            
        else:
            msm.show_err_msg('data', msm.msm_wrong_type)
        
        
    def find_profile_cl(self, **keywords):
        """
        find the data meeting the requirements 
        
        Inputs:
        ----------------------------------------
        
        1. keywords:<dict>: list of criteria 
        
        Returns:
        ------------------------------------------
        1. found_data:<list, t:ctm_field_cl>: datasets meeting the criteria
        
        """
        
        found_data=[]
        ifound=0
        
        for data in self.data_lst:
            # check if the data meeting the criteria
            matched=data.is_matched(**keywords)
            if (matched):
                ifound=ifound+1
                found_data.append(data)
        if (ifound>0):
            return found_data
        else:
            return None
        
    def add_attr(self, name, val):
        
        """
        add or replace attributes
        Inputs
        ---------------------------------
        1. name:<str>: attribute name
        2. val: <obj>: value of the attribute
        
        
        """
        
        if (name=='mask_val'):
            self.mask_val=val
            self.attr_dict.update({name:val})
            
        
        elif (name=='time_tag'):
            self.time_tag=val
            self.attr_dict.update({name:val})
        
        elif (name=='name'):
            self.name=name
            self.attr_dict.update({name:val})
            
        elif ('lon'==name):
            # make things constant
            self.set_lon(val)
        
        elif ('lat'==name):
            # make things constant
            self.set_lat(val)
            
        else:
            self.attr_dict.update({name:val})
            
        
    
    def get_attr(self, name):
        
        """
        check attribute
        Inputs:
        ---------------------------------------------
        1. name:<str>: attribute name
        
        Returns:
        ----------------------------------------------
        1. val:<obj>: value of the attribute
        
        """
    
        if (name in self.attr_dict):
            return self.attr_dict[name]
        else:
            msm.show_err_msg(name, msm.msm_no_attr)
            return None
    
    
    def get_lon(self):
        
        """
        get lon
        Returns:
        1. lon:<array, (nob,)>: longitude
        
        """
        return self.lon
    

    def get_lat(self):
        
        """
        get lat
        Returns:
        ----------------------------------------------------------------------------------
        1. lat:<array, (nob,)>: latitude
        
        """
        return self.lat

    
    def get_pressure(self):
        
        """
        get pressure
        Returns:
        -------------------------------------------------------------------------------------------------
        1. pressure:<array, (nob, nz)>: longitude
        
        """
        return self.pressure
    
    def init_veritcal_interp(self, opres, is_in_log, do_ob_reverse):
        
        """
        initialize vertical interpolation
        
        Inputs
        ----------------------------------------------------------------------------------------------------
        1. opres:<array, (nlon, ob_nz)>: ob pressure grid to be projected. 
        2. is_in_log:<T/F>: optional. Ture if opres is given in log10 space
        3. do_ob_reverse:<T/F>: optional. True of opres is given at descending order 

        """
        
        if (oob.check_ot_type(self.vintpl)==oob.ot_vintpl):
           
            self.vintpl.init_interp(opres)
            
        else:
        
            msg='the vertical pressure is not set correctly'
            msm.show_err_msg(msg)
        
        return  self.vintpl
        
    
    
    
    def is_data_in_list(self, data, keys):
        """
        check whether a data is already in the self.data_lst 
        
        
        Inputs:
        -----------------------------------------------------------------------------
        1. data:<ctm_profile_cl>: profile fileds
        2. keys:<list, t:string>: a list of key should be checked. 
        
        Returns:
        ----------------------------------------------------------------------------------------
        1. is_found:<integer>: if a data set with the same attributes of data input  
        is found in the list, its index will be given; 
        Otherwise, oob.fill_val_int will be given 
        """
        
        is_found=oob.fill_val_int
        search_list={}
        
        if  (oob.get_ot_type(data)<>oob.ot_profile):
            return is_found
        
        for keyname in keys:
            if (keyname in data.attr_dict):
                keyval=data.get_attr(keyname)
                search_list.append({keyname:keyval})
                
            else:
                return is_found
        i=0
        for data_obj in self.data_lst:
            is_match=data_obj.is_matched(data.name, **search_list)
            if (is_match):
                is_found=i
                return is_found
            
            i=i+1
        
        
            

    def __delitem__(self, index):
        """
        del item given at index
        Inputs: 
        1. index: <integer/string>: index or name of the data to be deleted 

        Returns:
        
        1. ndel:<integer>: number of the data have been deleted 
        """
        
        otype=oob.get_ot_type(index)
        old_ndata=len(self.data_lst)
        
        if (otype==oob.ot_string):
            idx_st=0
            idx=0
            
            # T1 delete all items with the given name
            
            while (idx_st<self.ndata):
                idx=idx_st
                
                for data in self.data_lst[idx_st:]:
                    matched=data.is_matched(name=index)
                    print idx, data.name, matched
                    
                    if (matched):
                        
                        tracer_id=data.get_attr('tracer_id')
                        print 'to be deleted',  tracer_id
                        
                        del self.data_lst[idx]
                        
                        self.ndata=len(self.data_lst)
                        
                        print self.ndata
                        
                        
                        idx_st=idx
                        idx=0
                        break
                    else:                    
                        idx=idx+1
                
                if (idx>idx_st):
                    idx_st=idx
        
        elif (otype==oob.ot_int): # by name 
            del  self.data_lst[idx]
            
        elif (otype==oob.ot_slice): # by slice
            del self.data_lst[index]
            
        else:
            msm.show_err_msg(index, msm.msm_wrong_index)
            
        self.ndata=len(self.data_lst)
        
        return old_ndata-self.ndata
    
    

    def __getitem__(self, index):
        """ retrieve data given at index
        Inputs: 
        1. index: <integer/string>: index or name of the data to be retrieve 
                
        Returns:
        
        1.data_lst:<list, t:ctm_profile_cl>: list of profiles with index and name 
        
        
        """
        
        otype=oob.get_ot_type(index)
        if (otype==oob.ot_string): # by name 
            data_lst=self.find_profile_cl(name=index)
            if (data_lst<>None):
                return data_lst
            else:
                msm.show_err_msg(index, msm.msm_wrong_index)
                return None
            
        elif (otype==oob.ot_int): # integer
            
            if (index<self.ndata):
                data_lst=self.data_lst[index]
                return data_lst
            
            else:
                msm.show_err_msg(index, msm.msm_wrong_index)
                return None
        elif (otype==oob.ot_slice): # slice
            return self.data_lst[index]
        
        else:
            msm.show_err_msg(index, msm.msm_wrong_index)
            return None
        
    
        


                                  
if (__name__=="__main__"):
    print 'For test, see ctm_world_m.py'
    
