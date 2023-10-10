""" Class for two,  three or four dimenisonal (lon-lat-press-time etc) model data set

    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.12.16
    
    Classes:
    ----------------------------------------------------
    1. ctm_field_cl: class for model 2D/3D/4D field. 
    
""" 

import numpy as npy

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm


# the model grid can be different types. 

import ctm_grid_3d as cgrd
import gc_grid_3d as gcgrd
import ctm_profile_m as cprof


             
class ctm_field_cl:
    
    """ general container for geo-physical variables in a high-dimension grid 
    
    It can be  used to store multiple species storage species fields, 
    such as T, O3, H2O, etc.  It also contains member functions for 
    field regrid and interploations etc
    
    Members:
    --------------------------------------------------------------------   
    1. name: <str>: name (container name etc)
    2. ctm_grid:<ctm_grid_cl>: the grid for ctm
    3. is_xyz_data:<T/F>: True when data is in form of  [lon, lat, lz,...]. 
    ---Otherwise, data is inform of [lon, lat, ...].
    
    4. ndim: <integer>: number of dimensions.
    5. dims: <array>: sizes in each directions.
    6. extra_axis:<gp_axis>: extra axis.
    7. attr_dict:<dict>: dictionary for attributes.
    8. data: <array>: data
    9. unit:<str/float>: unit of the data
    

    
    Attributes (expected or reserved)
    ---------------------------------------------------
    1. tracer_name:<str>: tracer name 
    2. tracer_id:<integer>: tracer number
    3. category:<str>: tracer category
    4. group:<str/integer>: group id for the tracer
    
    5. tau_st:<float>: start time (hours since 1985.01.01 00:00)
    6. tau_end:<float>:end time (hours since 1985.01.01 00:00)
    
    
    
    Functions:
    ------------------------------------------------------

    1. __init__:Initialization
    2. get_grid: get CTM grid 
    
    # data retieval and maintainence. 
    3. update_data: update the data set and its attributes
    4. copy_data: get the copy of the data set.
    5. get_data: return the pointer of data member.
    5. __getitem__: (over-riding) index function 
    6. get_profile: get profiles at horizontal locations
    
    # attributes
    7. set_attr: set attribute
    8. get_attr: get attribute
    9. is_matched: compare attributes
    10. copy_attr_dict:  get attriute dictionary 
    
    11. copy: make a copy of the class
    
    
    """ 
    
    def __init__(self, name, data,  \
                     unit="",\
                     ctm_grd=None, \
                     is_xyz_data=True,\
                     extra_axis=None,\
                     **keywords):
        
        """ initialization
        
        Input:
        
        --------------------------------------
        1. name:<str>:name of the object 
        2. data:<array>:  (2D, 3D or 4D) data over ctm_grid
        ---Data is assumed to be in the form of
        
        (a). (lon, lat, lz) for 3D spatial distributions  or
        (b). (lon, lat, ntime) for 2D time serie, 
        (c). (lon, lat, ne) 2D data ensemble
        (d). (lon, lat, lz, ntime):3D time series
        (e). (lon, lat, lz, ne): 3D data ensemble 
        
        3. ctm_grd:<ctm_grid_cl>: grid 
        4. is_xyz_data:<T/F>: form of the field 
        
        5. keywords:<dict>:  attributes (expected):
        
        --->tracer_name:<str>: tracer name 
        --->tracer_id:<integer>: tracer number
        --->category:<str>: tracer category
        --->group:<str/integer>:tracer group ID. 
        --->tau0:<float>: start time (Hours since 1985.01.01 00:00)
        --->tau1:<float>:end time (Hours since 1985.01.01 00:00)
        
        """
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_field})
        
        # S1: 'Name and unit 
        name=name.strip()
        self.name=name
        self.attr_dict.update({'name':self.name})
        
        self.unit=unit
        
        self.attr_dict.update({'unit':self.unit})
        
        # S2: 'Data'

        self.data=None
        self.dims=None 
        self.ndim=None
        
        self.update_data(data)
        
        self.is_xyz_data=is_xyz_data
        
        # S3: ctm_grd
        self.ctm_grd=None
        
        self.set_grid(ctm_grd)
        
        # S4: assign extral axis (such as time, or ensemble index etc)
        self.extra_axis=None
        
        if (oob.check_type(extra_axis, oob.ot_axis)):
            self.extra_axis=extra_axis
        else:
            if (extra_axis<>None):
                msm.show_err_msg('extra axis', msm.msm_wrong_type)
        
        # S5: attribute     
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
        
            
        
    def set_grid(self, ctm_grd):
        """ define the grid 
        Inputs: 
        -------------------------------
        1.ctm_grid:<ctm_grid_cl>: grid
        """ 
        
        
        if (ctm_grd==None):
            #  re-set model grid to None
            
            if (self.ctm_grd<>None):
                del self.ctm_grd
            
            self.ctm_grd=None
        else:
            #  set model grid 
            
            if (oob.check_type(ctm_grd,  oob.ot_grid)):
                self.ctm_grd=ctm_grd.copy()
                
            else:
                msm.show_err_msg('ctm grd', msm_wrong_type)
        
        
        
    def get_grid(self):

        """ retrieve CTM grid
        
        Returns: 
        -----------------------------------
        1. ctm_grd:<ctm_grd_cl>: class
        """
        
        return self.ctm_grd
    
    
    def update_data(self, data, \
                        **keywords):
        """
        change  data 
        
        Inputs: 
        ----------------------------------------
        1.data:<array, ([nlon], [nlat], [nz], ,,,)>: 
        ---data array with shape assumed to be consistent with ctm_grid
        
       
        """
        if (data==None):
        
            self.data=None
            self.dims=None
            self.ndim=None
        
        else:
            if (oob.check_type(data, oob.ot_array) | oob.check_type(data, oob.ot_list)):
                self.data=npy.array(data) # 'copy' data
                self.dims=npy.shape(self.data)
                self.ndim=npy.size(self.dims)
            else:
                self.data=data
                self.dims=npy.shape(data)
                self.ndim=size(self.dims)
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
        
    def set_attr(self, name, val):
        """
        set attribute
        
        Inputs:
        -----------------------------
        1. name:<str>: attribute name
        2. val: <obj>: attribute value
        
        """
        # make dictionary consistent with member variables. 

        
        if (name=='name'):
            self.name=val
            self.attr_dict.update({'name':self.name})
        
        elif (name=='unit'):
            self.unit=val
            self.attr_dict.update({'unit':self.unit})
        else:
            self.attr_dict.update({name:val})
        
    def get_attr(self, name):
        
        """  get attribute
        
        Inputs:
        -----------------------------
        1. name:<str>: attribute name
        
        Returns:
        -------------------------------
        1. val: <obj>: attribute value
        """
        
        
        if (name in self.attr_dict):
            return self.attr_dict[name]
        else:
            msm.show_err_msg(name, msm.msm_no_attr)
            return None
        

        
    def copy_data(self):
        """
        get a copy of the data 
        
        Outputs:
       -------------------------------------------
       1. data:<array,>: one copy of the data stored here
       
       """
        if (oob.get_ot_type(self.data)==oob.ot_array):
            data=npy.array(self.data)

        elif (oob.get_ot_type(self.data)==oob.ot_list):
            data=list(self.data)
        else:
            data=self.data
        
        
        return data
    

    def get_data(self):
        
        """
        get the data 
        
        Outputs
        ----------------------------------
        1. self.data:<array,>: the data member 

        """
        return self.data
    
    



    def copy_attr_dict(self):
        """ get the attribute dictionary 
        """
        return dict(self.attr_dict)
    
    
    def __getitem__(self, index):
        """ over-ridding index function
        Inputs:
        ---------------------------------
        1. index:<array,(ndim, )>: element indexs of data array
        
        Returns:
        ---------------------------------
        1. data[index]:<array>: data values at the  elements specified by index. 
        
        """
        return self.data[index]
    
    
    def get_profile(self, olon, olat):
        
        """ get profiles at the horizontal locations  (olon, olat)
        
        Inputs:
        ---------------------------------
        1. olon:<array, (nob,)>: longitudes at nob locations
        2. olat:<array, (nob,)>: latitudes  at nob locations

        Returns:
        ---------------------------------
        1. prof:<array, (nob, nz)>: profiles at nob locations
        
        """
        
        if (self.ctm_grd==None):
            msm.show_err_msg('No grid has been set')
            return None
        else:
            if (oob.check_type(self.ctm_grd, oob.ot_grid)):
                hintpl=self.ctm_grd.get_hintpl(olon, olat)
                prof=hintpl.get_profile(self.data)
                if (self.ndim==2): # just 1-D data will be return 
                    
                    prof=prof[:,0] 
                return prof
            else:
                
                msm.show_err_msg('ctm_grid has not been set properly')
                return None
        
        
    
    def is_matched(self, **keywords):
        """
        Check whether attribute is matched what is required
        Inputs:
        -------------------------------------------------
        1. keywords:<dictionary>: Attributes to be compared
        
        """
        
        matched=oob.is_attr_matched(self, **keywords)
        return matched
    
    def copy(self):
        
        """ make a copy of the class
        Returns:
        --------------------------------------
        1. new_field: a copy of field
        
        """
        
        fld=self.data
        
        # S1: attributes
        attrs=self.copy_attr_dict()
        
        # #T: remove possible duplicated attributes
        
        
        if ('name' in attrs):
            del attrs['name']
        elif ('ot_type' in attrs):
            del attrs['ot_type']
            
        
        
            
        # S2: duplicate the class
        
        new_field=ctm_field_cl(self.name, fld,  \
                                   is_xz_data=self.is_xz_data,\
                                   extra_axis=self.extra_axis, **attrs)
        return new_field
    

if (__name__=="__main__"):
    print '1. form a grid 4x5 with 47 levels:'
    cm=gcgrd.gc_grid_cl(0, 0, 47, mod_res='4x5')
    nx=cm.nlon
    ny=cm.nlat
    ps=npy.zeros([nx, ny], float)
    ps[:,:]=1000.0
    mod_pres=cm.compute_mod_pres(ps)
    
    print npy.shape(mod_pres)
    cf=ctm_field_cl('mod_pres', mod_pres, tra_name='pressure', tracer_id=3, 
                    unit='hPa')
    cf.set_grid(cm)
    
    olat=npy.arange(-60, 60, 10.0)
    is_match=cf.is_matched(tra_name='pressure', tracer_id=3)
    print 'is_macth', is_match
    
    print cf[3,4,:]
    
    olon=npy.zeros(npy.shape(olat), float)
    prof=cf.get_profile(olon, olat)
    print npy.shape(prof), npy.shape(olon)
    print prof[0,:]
    

    ctm2=cf.get_grid()
    print ctm2.nlon, ctm2.nlat, ctm2.nz
    print ctm2.ix, ctm2.jx, ctm2.lx
    sp=ctm2.get_attr('surf_pres')
    print npy.shape(sp)
    print 'profile interpolation'

    hinterp=ctm2.get_hintpl(olon, olat)
    ops=hinterp.get_profile(ps)
    
    print npy.shape(ops)

    
    print 'regrid'
    
    regrid_mtx=ctm2.get_regrid_mtx(olon, olat)
    ops=regrid_mtx.get_profile(mod_pres)
    print npy.shape(mod_pres)
    print npy.shape(ops)
    
    print ops[0,0,:]
    
 

    
    
    
               
