""" Class for horizontal (lon-lat) interpolations of model fields


    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.09.30



    1. Functions:
    ===============================================
    1. is_the_same_array: whether two arrays are identical.
    

    2. Classes:
    ==============================================
    1. hinterp_cl: class for horizontal interpolation
    
    
"""


import numpy as npy
import gp_axis_m as axis_m
# import flib as flb
import sample_model_field as grd_pf
import message_m as msm
import otool_obj as oob


def is_the_same_array(a, b):
    """ 
    check whether two arrays are identical. 
    Inputs:
    1. a: <array>: array of a
    2. b: <array>: array of b
    Outputs: 
    1. is_same:<T/F>: True if a==b
    """
    
    is_the_same=False
    if (npy.size(a)<>npy.size(b)):
        return is_the_same
    
    dx=npy.sum(npy.abs(a-b))
    
    if (dx==0):
        is_the_same=True
        
    return is_the_same

    
class hinterp_cl:
    
    """horizontal interpolation class
    
    Members:
    ----------------------------------------
    1.  ax_lon:<gp_axis>:longitude axis
    2.  ax_lat:<gp_axis>:latitude axis
    3.  dims:<tuple>: shape of [ax_lon, ax_lat]
    4.  lon:<array>: x value
    5.  lat:<array>: y value
    6.  nlon:<integer>: number of longitude
    7.  nlat:<integer>: number of latitude
    8.  mask_val:<float>: filling value for missing or bad data
    9. attr_dict:<dict>: attribute dictionary
    

    # ob (target) grid 
    
    10.  olon:<array>:  longitudes to be interploated to
    11.  olat:<array>: latitudes to be interploated to
    12.  ob_only_1d:<T/F>: if true, olon and olat are along observation tracks,  
    ---otherwise olon and olat form another (target) grid 
    
    # variables for saving interpolation coefficients.
    
    13.  lonp1:<array>: longitude index for left corners of the encircling boxes.
    14.  lonp2:<array>: longitude index for right corners of the encircling boxes.
    15.  latp1:<array>: latitude index for lower corners of the encircling boxes.
    16.  latp2:<array>: latitude index for upper corners of the encircling boxes.
    17.  lonwgt:<array>: weighting factor for left corner
    18.  latwgt:<array>: weighting factor for lower corner
    
    19  w1, w2, w3, w4:<array> aggregated weighting factor for the 4 corners. 
    20. is_ready: <T/F>: True if interpolation has been initilzed. 
    
    Functions:
    --------------------------------------------------------
    1.  __init__: initialization
    2. set_attr: set attribute
    3. get_attr: get_attribute
    4. set_mask: set mask values
   
    5. init_interp: set parameters for interoplation 
    6. get_profile: get the profiles at given locations (or regridding field when  ob_only_1d=False
    7. get_wgt: Interpolation weight along lon or lat 
    
    
    """
    
    def __init__(self, ax_lon, ax_lat, mask_val=oob.fill_val):
        """  initialize horizontal interpolation class 
        Inputs:
        
        1. ax_lon:<gp_axis>: longitude axis 
        2. ax_lat:<gp_axis>: latitude  axis 
        3. mask_val: <float>: filling value for missing or bad data
        
        """

        self.attr_dict={}
        self.attr_dict.update({'ot_type': oob.ot_hinterp})
        
        self.mask_val=-999.0
        self.set_mask(mask_val)
        
        
        self.ax_lon=ax_lon
        self.ax_lat=ax_lat
        
        # make a copy of the lon/lat axis values 
        
        self.lon=ax_lon[:]
        self.lat=ax_lat[:]
        
        self.nlon=npy.size(self.lon)
        self.nlat=npy.size(self.lat)
        
        self.dims=tuple([self.nlon, self.nlat])
        self.ndim=2
        
        


        # initialize ob (target) grid to None
        self.ob_lon=None
        self.ob_lat=None
        self.ob_only_1d=None

        self.lonp1=None
        self.lonp2=None
        self.latp1=None
        self.latp2=None
        
        self.w1=None
        self.w2=None
        self.w3=None
        self.w4=None
        
        # Not  ready for interpolations
        self.is_ready=False
        
    def set_mask(self, mask):

        """ set values for data mask 
        
        Inputs:
        ------------------------------------------
        1. mask:<float>: new mask values
        
        """

        self.mask_val=mask
        self.attr_dict.update({'mask_val':mask})

    
    def set_attr(self, name, val):
        
        """
        add or replace attributes
        
        Inputs
        ---------------------------------
        1. name:<str>: attribute name
        2. val: <obj>: value of the attribute
        
        
        """
        if (name=='mask_val'):
            self.set_mask(val)
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


    def init_interp(self, olon, olat, ob_only_1d=True):
        
        """  initialize interpolation coefficients 
        
        Inputs:
        --------------------------------------------------
        1. olon:<array>: longitudes to be projected to 
        2. olat:<array>: latitudes to be projected to 
        3. ob_only_1d:<T/F>: True if points are along tracks
        
        """
        # S1 check the size
        
        if (ob_only_1d):
        
            if (npy.size(olon)<>npy.size(olat)):
                msg='input olon and olat have different size'
                msm.show_err_msg(msg)
                return None
        
        
        self.ob_only_1d=ob_only_1d
        
        self.ob_lon=olon
        self.ob_lat=olat
        
        self.ob_nlon=npy.size(olon)
        self.ob_nlat=npy.size(olat)
        
        if (ob_only_1d):
            
            self.ob_dims=tuple([self.ob_nlon])
        else:
            self.ob_dims=tuple([self.ob_nlon, self.ob_nlat])
        
        # get the boundaries (corners) and weights 
        
        lonp1, lonp2, lonwgt=self.ax_lon.getwgt(olon, mask_val=self.mask_val)
        latp1, latp2, latwgt=self.ax_lat.getwgt(olat, mask_val=self.mask_val)
        
        self.lonp1=lonp1
        self.lonp2=lonp2
        
        self.latp1=latp1
        self.latp2=latp2
        
        
        
        if (ob_only_1d):
            # sample along track
            
            # set the weight for 4 corners 
            self.w1=latwgt*lonwgt
            self.w2=(1.0-latwgt)*lonwgt
            self.w3=latwgt*(1.-lonwgt)
            self.w4=(1.0-latwgt)*(1.0-lonwgt)
        
        else:
            #  regrid 
            
            self.w1=lonwgt  # for left longitude point  
            self.w2=(1.0-lonwgt)  
            self.w3=latwgt # for lower latitude point 
            self.w4=(1.0-latwgt)
        
        # ready to use
        self.is_ready=True
        
        return self.is_ready
    
    
    def get_profile(self, gp_fld):
        
        
        """
        Get profiles at the given horizontal locations
        Inputs:
        ------------------------------
        1. gp_fld:<array, (nlon, nlat, nl, ...>: 2D (or more)  gridded data. 
        
        Outputs:
        1. prof:<array, (ob_nlon, nl, ...): if ob_only_1d is true. or
        --- prof:<array, (ob_nlon, ob_nlat, ...): if ob_only_1d is False
        
        
        """
        
        # S1. check the data size
        
        dims=npy.shape(gp_fld)
        ndim_fld=npy.size(dims)
        
        if (not self.is_ready):
            msg='Not initialized'
            msm.show_err_msg(msg)
            return None
        
        if (not oob.check_ot_shape(gp_fld, self.dims)):
            msg=dims
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
        
        if (self.ob_only_1d):
            # do interpolation
            if (ndim_fld>4): 
                # high dimension
                # T1 reshape 
                
                new_dims=dims[0:3]+tuple([-1])
                gp_fld=npy.reshape(gp_fld, new_dims)
                
                prof=grd_pf.sample_field_em(gp_fld, self.lonp1,self.lonp2, \
                                                self.latp1, self.latp2, \
                                                self.w1, self.w2, self.w3, self.w4, \
                                                mask_val=self.mask_val)
                # reshape back to [np, ...]
                
                np=self.ob_nlon
                
                new_dims= tuple([np])+dims[2:]
                prof=npy.reshape(prof, new_dims)
                
                               
            elif (ndim_fld==4): 
                # multiple tracers
                
                prof=grd_pf.sample_field_em(gp_fld, self.lonp1,self.lonp2, \
                                                self.latp1, self.latp2, \
                                                self.w1, self.w2, self.w3, self.w4, \
                                                mask_val=self.mask_val)
                
            else:
                # 3D data
                
                # print 'max', max(gp_fld.flat)
                # print gp_fld[20, 20]
                
                prof=grd_pf.sample_field(gp_fld, self.lonp1,self.lonp2, \
                                             self.latp1, self.latp2, \
                                             self.w1, self.w2, self.w3, self.w4,  \
                                             mask_val=self.mask_val)
                if (ndim_fld==2):
                    prof=prof[:,0]
                
        else:
            # do 'regridded
            
            if (ndim_fld>3):
                # reshape to [nlat, nlon, -1]
                
                new_dims=dims[0:2]+tuple([-1])
                gp_fld=npy.reshape(gp_fld, new_dims)
                
                prof=grd_pf.regrid_field(gp_fld,self.lonp1,self.lonp2, self.w1, \
                                             self.latp1, self.latp2, self.w3, \
                                             mask_val=self.mask_val)

                new_dims=tuple([self.ob_nlot, ob_nlat])+dims[2:]
                
                prof=npy.reshape(prof, new_dims)
                
                
            else:
                prof=grd_pf.regrid_field(gp_fld,self.lonp1,self.lonp2, self.w1, \
                                             self.latp1, self.latp2, self.w3, \
                                             mask_val=self.mask_val)
                
                if (ndim_fld==2):
                    prof=prof[:,:,0]
        
                
        return prof
    

    def get_wgt(self, xval, ax_id=0):
        """ get interpolation weight
        A wrap for getwgt for each axis
        
        Inputs:
        -----------------------------------------------------
        1.xval:<array, (nob,)>: values 
        2. ax_id:<integer/string>: axid =0 or 1 or ax_id='lon', 'lat'
        

        Returns:
        -----------------------------------------------------
        1. pl:<array, (nob,)>: left boundaries
        2. pr:<array, (nob,)>: right boundaries
        3. wgt:<array, (nob,)>: weighting factor for left boundaries
        
        """
        if (ax_id==0):
            gax=self.ax_lon
        elif (ax_id==1):
            gax=self.ax_lat
        elif (ax_id=='lon'):
            gax=self.ax_lon
        elif (ax_id=='lat'):
            gax=self.ax_lat
        else:
            msg='Wrong axis id: '
            msm.show_err_msg(msg)
            return None, None, None
        
        pl, pr, wgt=gax.getwgt(xval,mask_val=self.mask_val)
        return pl, pr, wgt
    
    
    
    
    
        
        

        
                
        
#<<< TESTS >>>

    
if (__name__=='__main__'):
     
    print '>>> test 1: form a grid 10(lon)x10(lat) with 47 levels <<<'
    
    rlon=npy.arange(-180, 180.0, 10)
    axis_lon=axis_m.gp_axis_cl('lon', rlon, ax_unit='Deg')
    
    
    rlat=npy.arange(-90, 90.0, 10.0)
    axis_lat=axis_m.gp_axis_cl('lat', rlat, ax_unit='Deg')
    
    
    h_itpl=hinterp_cl(axis_lon, axis_lat)
    olon=npy.array([10.0, 20.0, 30.0])
    olat=npy.array([10.0, 20.0, 30.0])
    
    h_itpl.init_interp(olon, olat)

    print 'Interpolation'
    
    fx=npy.zeros([npy.size(rlon), npy.size(rlat)], float)
    
    for ix in range(npy.size(rlon)):
        for iy in range(npy.size(rlat)):
            fx[ix, iy]=rlon[ix]+rlat[iy]
    
    pf0=range(npy.size(olon))
    
    for ip in range(npy.size(olon)):
        pf0[ip]=olon[ip]+olat[ip]
    
     
    pf=h_itpl.get_profile(fx)
    print npy.shape(pf)
    print pf
    print pf0
    
    print '********* regrid'
    
    h_itpl.init_interp(olon, olat, ob_only_1d=False)
    
    pf=h_itpl.get_profile(fx)
    print npy.shape(pf)
    print pf
    
    
    
    
