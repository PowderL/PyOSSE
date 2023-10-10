""" container (class ) for two or three dimenisonal (lon-lat-time etc) data set
    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.10.28
    
    It is simplified version of ctm_field_m.py
    
    Classes:
    --------------------------------------------------
    1. ctm_field_2d: general container for geo-physical variables defined at a 2D-grid
   

""" 
import numpy as npy

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
# model grid can have different type

import ctm_grid_2d as cgrd
import gc_grid_2d as gcgrd
import ctm_profile_m as cprof



class ctm_field_2d:
    
    """ general container for geo-physical variables defined in a 2D-grid 

    It can be used to store multiple surface fields, and contains  transform and interploations etc
   
    Members:
    ------------------------------------------------------------------------------
    1. name: <str>: name (container name etc)
    2. ctm_grid:<ctm_grid_cl>: the grid for ctm field. 
    3. is_xy_data:<T/F>: True when data is in form of  [lon, lat, ...]
    
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
    
    5. tau_st:<float>: start time (Seconds since 1985.01.01 00:00)
    6. tau_end:<float>:end time (Seconds since 1985.01.01 00:00)
    
    
    
    Functions:
    -----------------------------------------------------------------------------

    1. __init__:Initialization
    2. get_grid: get CTM grid 
    
    # data retieval and maintainence. 
    3. update_data: update the data set and its attributes
    4. copy_data: get the copy of the data set.
    5. get_data: return the class data member.
    5. __getitem__: (over-riding) index function 
    6. get_profile: get profiles at the horizontal locations
    
    # attributes
    7. set_attr: set attribute
    8. get_attr: get attribute
    9. is_matched: compare attributes
    10. copy_attr_dict:  get attriute dictionary 
    11. copy: make a copy of the class
    



    # attribute functions
    
    
    """ 
    
    def __init__(self, name, data,  \
                     unit="",\
                     ctm_grd=None, \
                     is_xy_data=True,\
                     extra_axis=None,\
                     **keywords):
        
        """ initialization
        
        Input:
        --------------------------------------
        1. name:<str>:name of the object 
        2. data:<array>:  (2D, or 3D) data over ctm_grid
        Data is assumed to be in the form of
        --->(a). (lon, lat) for 2D spatial distributions  or
        --->(b). (lon, lat, ntime) for 2D time serie, 
        --->(c). (lon, lat, ne) 2D data ensemble
        
        3. ctm_grd:<ctm_grid_cl>: grid 
        4. is_xy_data:<T/F>:I
        
        5. keywords:<dict>: attributes
        reserved keywords:
        --->tracer_name:<str>: tracer name 
        --->tracer_id:<integer>: tracer number
        --->category:<str>: tracer category
        --->group:<str/integer>:tracer group 
        --->tau0:<float>: start time (hours since 1985.01.01 00:00)
        --->tau1:<float>:end time (hours since 1985.01.01 00:00)
        
        """
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_field})
        
        # S1 'Name and unit 
        name=name.strip()
        self.name=name
        self.attr_dict.update({'name':self.name})
        
        self.unit=unit
        
        self.attr_dict.update({'unit':self.unit})
        
        # S2 'Data'

        self.data=None
        self.dims=None 
        self.ndim=None
        
        self.update_data(data)
        
        self.is_xy_data=is_xy_data
        
        # S3 ctm_grd
        self.ctm_grd=None
        
        self.set_grid(ctm_grd)
        
        # S4 assign extral axis (such as time, or ensemble index etc)
        self.extra_axis=None
        
        if (oob.check_type(extra_axis, oob.ot_axis)):
            self.extra_axis=extra_axis
        else:
            if (extra_axis<>None):
                msm.show_err_msg('extra axis', msm.msm_wrong_type)
        
        # S5     
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
        
            
        
    def set_grid(self, ctm_grd):
        """ define the grid 
        Inputs: 
        ---------------------------------
        1.ctm_grid:<ctm_grid_cl>: grid
        """ 
        
        
        if (ctm_grd==None):
            if (self.ctm_grd<>None):
                del self.ctm_grd
            
            self.ctm_grd=None
        else:
            if (oob.check_type(ctm_grd,  oob.ot_grid)):
                self.ctm_grd=ctm_grd.copy()
                
            else:
                msm.show_err_msg('ctm grd', msm_wrong_type)
        
        
        
    def get_grid(self):

        """ retrieve CTM grid
        
        Returns: 
        -----------------------------------------------------
        1. ctm_grd:<ctm_grd_cl>: class
        """
        
        return self.ctm_grd
    
    
    def update_data(self, data, \
                        **keywords):
        """
        change data member
        
        Inputs: 
        ----------------------------------------
        1.data:<array, ([nlon], [nlat], [nz], ,,,)>: 
        data array with shape assumed to be consistent with ctm_grid
        2. keywords:<dict>: extra inputs
        
        
       
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
            msm.show_err_msg(name, msm.no_attr)
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
                    prof=prof[:]
                
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
        
        # attributes
        attrs=self.copy_attr_dict()
        
        # remove possible duplicated attributes
        
        
        if ('name' in attrs):
            del attrs['name']
        elif ('ot_type' in attrs):
            del attrs['ot_type']
            
        
        # S3 make new profiles
        
            
        # T1 duplicate the class
        new_profile=ctm__cl(self.name, fld,  \
                                is_xz_data=self.is_xz_data,\
                                extra_axis=self.extra_axis, **attrs)
        return new_profile
    

if (__name__=="__main__"):
    print '> form a grid 4x5 with 47 levels:'
    cm=gcgrd.gc_grid_2d(0, 0, mod_res='4x5')
    nx=cm.nlon
    ny=cm.nlat
    ps=npy.zeros([nx, ny], float)
    ps[:,:]=1000.0
    print 'nx, ny:', nx, ny
    
    cf=ctm_field_2d('mod_pres', ps, tra_name='pressure', tracer_id=3, 
                    unit='hPa')
    print '> construct the field'
    print 'name:', cf.name
    
    cf.set_grid(cm)

    print '> after set_grid'
    olat=npy.arange(-60, 60, 10.0)
    is_match=cf.is_matched(tra_name='pressure', tracer_id=3)
    print 'is_macth', is_match
    
    print cf[3,4]
    
    print '>interpolate'
    
    olon=npy.zeros(npy.shape(olat), float)
    prof=cf.get_profile(olon, olat)
    print 'profile:', prof
    
    
    
    ctm2=cf.get_grid()
    print '> make a copy of the grid'
    print 'nlon, nlat:', ctm2.nlon, ctm2.nlat
    print 'nlon, nlat:', ctm2.ix, ctm2.jx
    
    print '>regrid to 12 (lon) x 12 (lat) '
    
    regrid_mtx=ctm2.get_regrid_mtx(olon, olat)
    ops=regrid_mtx.get_profile(ps)
    
    print 'shape of origin ps:', npy.shape(ps)

    print 'shape of new ops:', npy.shape(ops)
    
    print ops[0,0]
    
 

    
    
    
               
