""" Class ) for model profile
    Authors: L. Feng, Edinburgh University
    History: v0.5, 2012.06.28
    History: v0.95, 2012.11.03
    
    Classes:
    -------------------------------------------------------
    1. ctm_profile_cl: general container for profiles for geo-physical variables at given locations. 
    

""" 
import numpy as npy

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm

import ctm_grid_3d as cgrd
import gc_grid_3d as gcgrd
import ctm_field_m as cfm
import ESA.util.vertical_profile as vpf


             
class ctm_profile_cl:
    
    """ general container for profiles for geo-physical variables at given locations. 
    
    Members:
    ------------------------------------------------------------------------------
    1. name: <str>: name (container name etc)
    2. olon:<array, (nob,)>: longitudes
    3. olat:<array, (nob,)>: latitudes 
    4. tau:<array, (nob,)>: time 
    
    5. dims:<array, (ndim,)>: shape of data
    6. ndim: <integer>: number of dimensions.
    7. extra_axis:<gp_axis_cl>: extra axis.
    8. attr_dict:<dict>: dictionary for attributes.
    9. data: <array, (nobs, nz, )>: data
    10. unit:<str/float>: unit of the data
    11. is_xz_data:<T/F>: True, if it is in shape of (nob, nz).
    12. mask_val:<float>: filling for missing or bad data
    


    Attributes (expected)
    ---------------------------------------------------
    1. tracer_name:<str>: tracer name 
    2. tracer_id:<integer>: tracer number
    3. category:<str>: tracer category
    4. group:<str/integer>: group id for the tracer
    
    5. tau_st:<float>: start time (Seconds since 1985.01.01 00:00)
    6. tau_end:<float>:end time (Seconds since 1985.01.01 00:00)
    
    
    
    Functions:
    ------------------------------------------------------------------------------

    1. __init__:Initialization
    2. get_location: retunr olon, olat
    
    
    #  data retieval and maintainence. 
    3. update_data: update the data set and its attributes
    4. copy_data: get the data set.
    5. __getitem__: (over-riding) index function 
    
    # attributes
    6.  set_attr: set attribute
    7.  get_attr: get attribute
    8. is_matched: compare attributes
    9. regrid_profile: vertically regrid profile. 
    10. copy:   duplicate itself
    
      
    
    """ 
    
    def __init__(self, name, data,  \
                     unit="",\
                     lon=None, \
                     lat=None, \
                     tau=None,\
                     mask_val=oob.fill_val,\
                     is_xz_data=True,\
                     extra_axis=None,\
                     
                     **keywords):
        
        """ initialization
        
        Input:
        --------------------------------------
        1. name:<str>:name of the object
        2. data:<array>:  (2D, 3D) data over location
        Data is assumed to be in the form of
        
        ---(a). (nob, lz) for cross-section
        ---(b). (nob, lz, ntime) for 2D time serie, 
        ---(c). (nob, lz, ne):2D data ensemble
        ---(d). (nob, lz, ne, ...): high-dimension data 
        
        3. unit:<str>: the unit of the profile 
        4. lon:<array, (nob,)>:longitudes
        5. lat:<array, (nob,)>:latitudes.
        6. tau:<array, (nob,)>: tau 
        7. is_xz_data:<T/F>: Ture when data is in the shape of  (nob, lz, ..)        
        8. extra_axis:<gp_axis>: extra axis (for example, time, etc)
        
        
        keywords (expected):
        
        1. tracer_name:<str>: tracer name 
        2. tracer_id:<integer>: tracer number
        3. category:<str>: tracer category
        4. group:<str/integer>:tracer group ID. 
        
        5. tau_st:<float>: start time (Seconds since 1985.01.01 00:00)
        6. tau_end:<float>:end time (Seconds since 1985.01.01 00:00)
        
        """
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_profile})
        
        # S1 'Name and unit 
        name=name.strip()
        self.name=name
        self.attr_dict.update({'name':self.name})
        
        self.unit=unit
        self.attr_dict.update({'unit':self.unit})
        
        # S2 'Data'
        self.lon=lon
        self.lat=lat
        self.tau=tau
        
        self.data=None
        self.dims=[]
        self.ndim=0
        self.mask_val=mask_val
        
        
        self.set_lon(lon)
        self.set_lat(lat)
        self.set_tau(tau)
        
        
        self.update_data(data)
        
        self.is_xz_data=is_xz_data
        
        
        # S3 assign extral axis (such as time, or ensemble index etc if interpolation is needed)
        
        self.extra_axis=None
        
        if (oob.check_type(extra_axis, oob.ot_axis)):
            self.extra_axis=extra_axis
        else:
            if (extra_axis<>None):
                msm.show_err_msg('extra axis', msm.msm_wrong_type)
            
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
        
       
    
    def update_data(self, data, \
                        **keywords):
        
        if (oob.get_ot_type(data)==oob.ot_array):
            self.data=npy.array(data)
        
        elif (oob.get_ot_type(data)==oob.ot_list):
            self.data=list(data)
        
        else:
            self.data=data
            
        self.dims=npy.shape(data)
        self.ndim=npy.size(self.dims)
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)
    
    
    def set_lon(self, olon):
        """ set longitude 
        Inputs:
        -------------------------------------------
        1. olon:<array>: longitude
        """
        
        self.lon=olon
        
    def set_lat(self, olat):
        
        """ set longitude 
        Inputs:
        -------------------------------------------
        1. olat:<array>: latitudes
        """
        
        self.lat=olat
        # self.attr_dict.update({'lat':self.lat})
        
    def set_tau(self, tau):

        """ set tau 
        
        Inputs:
        -------------------------------------------
        1. tau:<array>: time  
        """

        self.tau=tau
        
        
    def set_attr(self, name, val):
        """
        set attribute
        Inputs:
        -----------------------------
        1. name:<str>: attribute name
        2. val: <obj>: attribute value
        
        """

        
        
        if (name=='name'):
            self.name=val
            self.attr_dict.update({'name':self.name})
        
        elif (name=='unit'):
            self.unit=val
            self.attr_dict.update({'unit':self.unit})
        else:
            self.attr_dict.update({name:val})
            
        
    
    def get_attr(self, name):
        """ set attribute
        Inputs:
        -----------------------------
        1. name:<str>: attribute name

        Returns:
        1. val:<obj>: attribute value
        
        """

        
        if (name in self.attr_dict):
            return self.attr_dict[name]
        else:
            msm.show_err_msg(name, msm.msm_no_attr)
            return None
    
    def get_lon(self):
        
        """
        get lon  
        
        """
        
        if (oob.get_ot_type(self.lon)==oob.ot_array):
            return npy.array(self.lon)
        else: 
            return self.lon
    
    def get_lat(self):

        """
        get lat  
        
        """

        if (oob.get_ot_type(self.lat)==oob.ot_array):
            return npy.array(self.lat)
        else: 
            return self.lat
        
    
    def copy_attr_dict(self):
        """
        get the attribute dictionary 
        
        """
        return dict(self.attr_dict)
    
    
    def copy_data(self):
        if (oob.get_ot_type(self.data)==oob.ot_array):
            data=npy.array(self.data)
        elif (oob.get_ot_type(self.data)==oob.ot_list):
            data=list(self.data)
        else:
            data=self.data
            
        return self.data
    
    
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
    
    
    def is_matched(self, **keywords):
        """
        Check whether attribute is matched what is required
        Inputs:
        -------------------------------------------------
        1. keywords:<dictionary>: Attributes to be compared
        
        """
        
        matched=oob.is_attr_matched(self, **keywords)
        return matched
    
    
    def regrid_profile(self, vintpl):
        """ create a new profile at level specified in vintp
        
        Inputs:
        ------------------------------------------
        1. vintpl:<None, or initialized vinterpl_cl>: interpolation function 
        --- Pressure grid of vintpl should be in the shape of (nob, nz), 
        consistent with the first two dimensions of profile.
        
        --- ob pressure of the vintpl represent the grid to be projected to. 
        should be in the shape of (nob, nz_new)
        
        Outputs:
        --------------------------------------------
        1. new_profile:<ctm_profil_cl>: the profile class at the new vertical grid 
        
        """
        
        
        # S1 data 
        fld=self.data
        
        # attributes
        attrs=self.copy_attr_dict()

        # remove possible duplicated
        
        if ('name' in attrs):
            del attrs['name']
        elif ('lon' in attrs):
            del attrs['lon']
            
        elif ('lat' in attrs):
            del attrs['lat']
        
        elif ('ot_type' in attrs):
            del attrs['ot_type']
            
        elif ('tau' in attrs):
            del attrs['tau']
        
        # S3 make new profiles
        
        # T2 interpolation
        
        new_data=vintpl.interpolate_mod_prof(self, fld)
            
        # constructure the profile with new data 
        new_profile=ctm_profile_cl(self.name, new_data,  \
                                      lon=self.lon, \
                                      lat=self.lat, \
                                      tau=self.tau,\
                                      mask_val=self.mask_val,\
                                      is_xz_data=self.is_xz_data,\
                                      extra_axis=self.extra_axis, **attrs)
        
        
        return new_profile
    

    def shape(self):
        """ get the shape of the data
        """
        return self.dims
    
    def size(self):
        """ get the size of the data
        """

        nsize=1
        for ix in self.dims:
            nsize=nsize*ix
        

        return nsize
    
    
    def grid_profile(self, lons=None, lats=None, use_intpl=0):
        """
        gridded profiles to a gridded defined by lons, lats 

        Inputs:
        -------------------------------------------------
        1. lons:<array, (nlon,)>: longitude grid 
        2. lats:<array, (nlat,)>: latitude grid 
        3. use_intpl:<int>: 0 means profiles will be allocated the closest grid points. 
        --- others, meaning interpolation will be done. 
        
        Outputs:
        ------------------------------------------------
        1. counts:<array, ([nlon], [nlat])>: number of profile at each grid boxes
        2. sum_val:<array, ([nlon], [nlat], nz)>: summary of profile at each grid boxes
        3. sum_square:<array, ([nlon], [nlat], nz)>: summary of profile^2 at each grid boxes
        

        Notes:
        ------------------------------------------------------
        if data have a more than 2 dims, it will be reshape to (nob, nz) with nz=size(data)/nob
        
        """
        
        if (self.ndim>2):
            # reshape to 2D data
            
            profile=npy.reshape(self.data, [self.nob, -1])
        else:
            profile=self.data
        

        if (lons==None):
            # only lats is given 
            counts, sum_val, sum_square=vpf.grid_profile_1d(self.lat, profile, \
                                                                lats,  use_intpl, self.mask_val)
        elif (lats==None):
            # only lons is given 
            counts, sum_val, sum_square=vpf.grid_profile_1d(self.lon, profile, \
                                                                lons,  use_intpl, self.mask_val)
        else:
            # both 
            
            counts, sum_val, sum_square=vpf.grid_profile_2d(self.lon, self.lat, profile, \
                                                                lons,  lats, use_intpl, self.mask_val)
        return counts, sum_val, sum_square
    
     
    def copy(self):
        """ create duplicate itself
        
        Outputs:
        --------------------------------------------
        1. new_profile:<ctm_profil_cl>: the profile class at the new vertical grid 
        
        """
        
        
        # S1 data 
        
        fld=self.data
        
        # attributes
        attrs=self.copy_attr_dict()
        
        # remove possible duplicated
        
        if ('name' in attrs):
            del attrs['name']
        elif ('lon' in attrs):
            del attrs['lon']
            
        elif ('lat' in attrs):
            del attrs['lat']
        
        elif ('tau' in attrs):
            del attrs['tau']

        elif ('ot_type' in attrs):
            del attrs['ot_type']
            
        
        # S3 make copy of the profile
        
        new_profile=ctm_profile_cl(self.name, fld,  \
                                       lon=self.lon, \
                                       lat=self.lat, \
                                       tau=self.tau,\
                                       mask_val=self.mask_val,\
                                       is_xz_data=self.is_xz_data,\
                                       extra_axis=self.extra_axis, **attrs)
        return new_profile

    

if (__name__=="__main__"):
    print '>1. form a grid 4x5 with 47 levels:'
    cm=gcgrd.gc_grid_cl(0, 0, 47, mod_res='4x5')
    nx=cm.nlon
    ny=cm.nlat
    ps=npy.zeros([nx, ny], float)
    ps[:,:]=1000.0
    mod_pres=cm.compute_mod_pres(ps)
    
    cl_cf=cfm.ctm_field_cl('mod_pres', mod_pres, tra_name='pressure', tracer_id=3, 
                    unit='hPa')
    cl_cf.set_grid(cm)
    olat=npy.arange(-60, 60, 10.0)
    is_match=cl_cf.is_matched(tra_name='pressure', tracer_id=3)
    print is_match
    
    
    olon=npy.zeros(npy.shape(olat), float)
    print '>2 get profile at locations' 
    
    prof=cl_cf.get_profile(olon, olat)
    
    print '>3 construct profile class'
    
    cl_prof= ctm_profile_cl('Pressure', prof,  \
                                unit="hPa",\
                                lon=olon, \
                                lat=olat, \
                                tau=0.0,\
                                is_xz_data=True,\
                                extra_axis=None)
    
    print '  ndim:', cl_prof.ndim
    print '  dims:', cl_prof.dims
    print '  name:', cl_prof.get_attr('name')
    
    print '>4 set unit'
    
    
    cl_prof.set_attr('unit', 'hPa')
    
    
    print '>5 make a copy'
    
    new_prof=cl_prof.copy()
    # test the copy function
    
    print ' reset the unit' 

    new_prof.set_attr('unit', 'tx')
    
    
    print '  new unit:', new_prof.unit
    print '  old unit:', cl_prof.unit
    

    print '>6 grid profiles to model grid ' 
    
    lons=cm.get_lon()
    lats=cm.get_lat()
    
    counts, sum_val, sum_square=new_prof.grid_profile(lons=lons, lats=lats)
    
    print npy.shape(counts)
    print 'sum counts', npy.sum(counts.flat)/47
    
    print npy.shape(sum_val)
    
    


    
    

    
    
    
               
