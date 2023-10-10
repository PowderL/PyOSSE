""" Class for 2D (lon, lat) grid used to store surface fields

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.09.30
    
    Classes:
    ---------------------------------------------------
    
    1. ctm_grid_2d:  derived class for 2D model grid.
    

"""

import ESA.util.gp_axis_m as axis_m
import ESA.util.horizontal_interp_m as hintpl_m
import ESA.util.vertical_interp_m as vintpl_m
import ESA.util.gp_grid_m as grid_m
import ESA.util.gp_axis_m as axis_m
import ESA.util.geo_constant as gc
import ESA.util.otool_obj    as oob
import ESA.util.pres_m       as pm
import compute_gc_grid       as cmg

import numpy as npy


class ctm_grid_2d(grid_m.gp_grid_cl): 
    
    """ derived class for model grid 
    It has grid of two axis named as 
    ['lon', 'lat']
    
    
    Members:
    ---------------------------------------
    # Added or over-ridden members
    
    1. attr_dict:<dict>: dictionary for attributes
    2. htype: Type of horizontal grid
    3. hresolution: resolution of the horizontal grid
    4. nlon:<integer>: size of longitude grid
    5. nlat:<integer>: size of latitude grid
    6. mask_val:<float>: filling for bad or missing data
    
    
    Functions (overridden and the new ones):
    -------------------------------------------
    
    1. __init__: initialize the class
    2. set_mask: set mask values
   
    # horizontal centres
    
    3. get_lon: return longitude grid 
    4. set_lon: set longitude grid 
    5. get_lat: return latitude grid 
    6. set_lat: set longitude grid 
    
    
    # edges 
    
    
    7. set_lon_edge: set longitude edge
    8. set_lat_edge: set latitude edge
    9. compute_lat_edge: compute_lat_edge
    10. compute_lon_edge: compute_lon_edge 
    11. get_area: calculate areas of grid boxes 
    
    # interploation
    12. get_hintpl: construct class for horizontal interpolation
    13. get_regrid_mtx: construct interpolation class for regridding fields to a new (lon, lat) grid
    
    
    Attributes (reserved):
    -----------------------------------------------------
    
    # grid edge 
    1. lon_edge:<array>: longitude grid 
    2. lat_edge:<array>: latitude grid
        
    # filling values
    3. mask_val:<float>: filling for bad or missing values. 
    
    
    
    
    
    """
    
    def __init__(self, lon, lat,htype='Gaussian', \
                     hresolution='4x5', \
                     mask_val=oob.fill_val, **keywords):
        
        """initialize grid 
        Inputs:
        ------------------------------------------------------
        1. lon, lat:<array>: array for model grids.
        ---their sizes are nlon, nlat, nz, respectively

        2. vtype:<string>: type of the vertical grid 
        
        3. htype:<string>: type of the horizontal grid
        4. hresolution:<string>: type of the horizontal grid 
        
        5. mask_val=oob.fill_val,

        6. keywords:<dict>: attributes 
        Reserved keywords:
        
        --->lon_edgee:<array, (nlon+1,)>: array of model lon edges
        --->lat_edgee:<array, (nlat+1,)>: array of model lon edges
        --->ap, bp:<array, (nz+1,)>: coefficients for pressure level
        --->surf_pres:<array, (nlon, nlat)>: pressures at surface
         
        
        """
        # S1: construct gp axis 
        
        
        ax_lon=axis_m.gp_axis_cl('lon', lon)
        ax_lat=axis_m.gp_axis_cl('lat', lat)
        
        self.htype=htype

        self.hresolution=hresolution
        
        # S2: initialize parent class 
        grid_m.gp_grid_cl.__init__(self, [ax_lon, ax_lat])
        
        self.mask_val=-999.0
        self.set_mask(mask_val)
        
        
        self.nlon=npy.size(lon)
        self.nlat=npy.size(lat)
        
        
        # S3: set attributes
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
            
        
        
    
    def set_mask(self, mask):
        """ set filling value for missing or bad data
        
        Inputs:
        -----------------------------------------
        1. mask:<float>: mask for missing or bad data
        
        
        """

        self.mask_val=mask
        self.set_attr('mask_val', mask)

        
    def get_lat(self):
        
        """ retrieve latitude grid 
        
        Returns: 
        -------------------------------------
        1. self.lat:<array, (nlat)>: latitudes 
        """
        # retrieve lat  grid
        
        ax_lat=self.get_axis('lat')
        lat=ax_lat[:]
        
        
        return lat

    def get_lon(self):
        """ retrieve longitude grid 
        
        Returns: 
        -------------------------------------
        1. self.lon:<array, (nlon)>: longitudes 
        """
        
        
        ax_lon=self.get_axis('lon')
        lon=ax_lon[:]
        
        return lon
    
    def set_lon(self, lon):

        """ set longitude  grid 
        
        Inputs:
        -------------------------------------
        1. lon:<array, (nlon)>: longitude 
        
        Returns: 
        -------------------------------------
        """
        
        
        # self.lon=npy.array(lon)
        
        self.nlon=npy.size(self.lon)
        ax_lon=axis_m.gp_axis('lon', lon)
        
        self.update_axis(ax_lon)
        # remove possible conflicts
 
        
        if ('lon_edge' in self.attr_dict):
            del self.attr_dict['lon_edge']
        
        
    def set_lat(self, lat):

        """ set latitude  grid 
        Inputs:
        ------------------------------------
        1. lat:<array>: latitude 
        ------------------------------------
        """
        
        
        # self.lat=npy.array(lat)
        
        self.nlat=npy.size(lat)
        
        ax_lat=axis_m.gp_axis('lat', lat)
        self.update_aixs(ax_lat)
        
        # remove possible conflicts

        if ('lat_edge' in self.attr_dict):
            del self.attr_dict['lon_edge']
            
            
    

    
   
    def get_regrid_mtx(self, newlon, newlat):
    
        """ construct interpolation classes for regridding fields to a new (lon, lat) grid
        
        Inputs:
        ------------------------------------
        1. newlon:<array, new_nlon>: longitude of the track
        2. newlat:<array, new_nlat>: latitude of the track
        
        
        Returns: 
        ------------------------------------
        hintpl:<interp_cl>: interpolation matrix for regridding

        """
        
        # S1 get axis 
        
        ax_lon=self.get_axis('lon')
        
        ax_lat=self.get_axis('lat')

        
        # S2 construct the interpolation matrix 
        
        hintpl=hintpl_m.hinterp_cl(ax_lon, ax_lat, mask_val=self.mask_val)
        
        # S3 set coefficients for interpolation 
        
        hintpl.init_interp(newlon, newlat,  ob_only_1d=False)
        
        return hintpl
    
    

    def get_hintpl(self, olon, olat):
        
        
        """ get weightings and locations for interpolating fields along a track 
        
        Inputs:
        ------------------------------------
        1. olon:<array, nob>: longitude of the track
        2. olat:<array, nob>: latitude of the track

        
        Returns: 
        ------------------------------------
        1. hintpl:<interp_cl>: interpolation matrix along the track 
        
        """
        
        # S1 get axis 
        
        ax_lon=self.get_axis('lon')
        
        ax_lat=self.get_axis('lat')

        # S2 construct the interpolation matrix 
        hintpl=hintpl_m.hinterp_cl(ax_lon, ax_lat, mask_val=self.mask_val)
        
        # S3 set coefficients for interpolation 
        
        hintpl.init_interp(olon, olat)
        
        return hintpl
    

    
    def compute_lon_edge(self, lon_st=-180.0, lon_end=180.0):
    
        """computer lontitude edge, and save it to the dictionary

        This function shall be over-ridden
                
        Inputs:
        ------------------------------------------
        1. lon_st:<float>: Starting point  of the longitude
        2. lon_end:<float>: End point of the longitude

       
        Returns:
        ------------------------------------------------------
        1. lon_edge:<array, nlon+1>: longitudes at model grid edges
        """
        
        lon_edge=npy.zeros(self.nlon+1)
        lon_edge[0]=lon_st
        ax_lon=self.get_axis('lon')
        lon=ax_lon[:]
        lon_up=lon[1:]
        lon_low=lon[0:-1]
        lon_edge[1:self.nlon]=0.5*(lon_up+lon_low)
        lon_edge[self.nlon]=lon_end

        self.set_lon_edge(lon_edge)
        
        return lon_edge
    
        
    def set_lon_edge(self, lon_edge):
        """
        save longitude box edges to the dictionary 
        
        """

        nedge=npy.size(lon_edge)
        
        if (nedge<>self.nlon+1):
            msg='Warning! The edge size is inconsistent with longitude size'
            msm.show_err_msg(msg)
        
        self.attr_dict.update({'lon_edge':lon_edge})
        
    
    def get_lon_edge(self, lon_edge):
        
        """
        
        Retrieve longitude box edges from dictionary 
         
        Returns:
        ------------------------------------------------------
        1. lon_edge:<array, nlon+1>: longitudes at model grid edges
        """
        

        if ('lon_edge' in self.attr_dict):
            
            lon_edge=self.attr_dict['lon_edge']
            return lon_edge
        
        else:
            
            msg='Longitudes at box edges have not been defined' 
            msm.show_err_msg(msg)
            return None
        
        
    def compute_lat_edge(self, lat_st=-90.0, lat_end=90):
        
        """
        Compute latitude box edge, and save it to the dictionary
        
        Inputs:
        ------------------------------------------
        1. lat_st:<float>: Starting point  of the latitude
        2. lat_end:<float>: End point of the latitude
        
        
        Returns:
        ------------------------------------------------------
        1. lat_edge:<array, nlat+1>: latitudes of model grid edges
        """
        
        lat_edge=npy.zeros(self.nlat+1)
        lat_edge[0]=lat_st
        
        ax_lat=self.get_axis('lat')
        lat=ax_lat[:]
        
        lat_up=lat[1:]
        lat_low=lat[0:-1]
        lat_edge[1:self.nlat]=0.5*(lat_up+lat_low)
        lat_edge[self.nlat]=lat_end
            
        self.set_lat_edge(lat_edge)
        
        return lat_edge
    
        
    def set_lat_edge(self, lat_edge):
        
        """
        Save latitude box edge into dictionary 
        
        Inputs:
        ---------------------------------------------------
        1. lat_edge:<array, nlat+1>: edge of the latitude boxes
        
        """
    
        nedge=npy.size(lat_edge)
        
        if (nedge<>self.nlat+1):
        
            msg='Warning! The edge size is inconsistent with longitude size'
            msm.show_err_msg(msg)

        self.attr_dict.update({'lat_edge':lat_edge})
        
    
    def get_lat_edge(self, lat_st=-90.0, lat_end=90):
    
        """
        
        Retrieve latitude box edges from dictionary 
        
        Returns:
        ------------------------------------------------------
        1. lat_edge:<array, nlat+1>: latitudes at model grid edges
        """
        
        
        if ('lat_edge' in self.attr_dict):
            
            lat_edge=self.attr_dict['lat_edge']
            return lat_edge
        
        else:
            
            msg='latitudes at box  edges have not been defined' 
            msm.show_err_msg(msg)
            return None
        
        
        
        
    def get_area(self,lon_st=-180.0, lon_end=180.0, \
                     lat_st=-90.0, lat_end=90):
        
        """
        Calculate areas for grid box (in m2)
        
        Inputs:
        --------------------------------------
        1. lon_st:<integer>: left edge of the first X grid box
        2. lon_end:<integer>: left edge of the last X grid box
        3. lat_st:<integer>: left edge of the first Y grid box
        4. lat_end:<integer>: left edge of the last Y grid box
        
        
        Returns:
        ------------------------------
        1. area:<array, (nx, ny)>: area of grid box in m^2
      
        """
        deg2rad=npy.pi/180.0
        # S1: get grid egdes 
        rlon=self.compute_lon_edge(lon_st=lon_st, lon_end=lon_end)
        rlat=self.compute_lat_edge(lat_st=lat_st, lat_end=lat_end)
        
        
        # S2: area
        
        area=cmg.get_area(rlon, rlat)
        
        return area

    
    
    def copy(self):

        """
        Make a copy of itself
        
        Returns:
        ---------------------------------------------------
        1. new_grid:<ctm_grid_cl>: new grid 
        
        """
        
        attrs=self.copy_attr_dict()
        
        if ('mask_val' in attrs):
            del attrs['mask_val']
        
        new_grid=ctm_grid_cl(self.lon, self.lat, \
                                 htype=self.htype, \
                                 hresolution=self.hresolution,\
                                 mask_val=self.mask_val, **keywords)
        return new_grid
    
    
            
if (__name__=='__main__'):
    
    print '>1. form a grid 4(lon)x5(lat) grid with 47 levels'
    
    rlon=npy.arange(-180, 180.0, 5)
    axis_lon=axis_m.gp_axis_cl('lon', rlon, ax_unit='Deg')
    
    
    rlat=npy.arange(-90.0, 90.0+1, 4.0)
    axis_lat=axis_m.gp_axis_cl('lat', rlat, ax_unit='Deg')
    

    axis_list=[axis_lon, axis_lat]
    
    mod_grid= ctm_grid_2d(rlon, rlat)
    
    lon_edge=mod_grid.compute_lon_edge()
    lon=mod_grid.get_lon()
    print '> lons:'
    
    print '0,-1:', rlon[0], rlon[-1]
    
    print 'lon edge [0,-1]:', lon_edge[0], lon_edge[-1]
    
    
    lat_edge=mod_grid.compute_lat_edge()
    
    print 'lat [0,-1]:', rlat[0], rlat[-1]
    print 'lat_edge [0,-1]:',  lat_edge[0], lat_edge[-1]
    
    print '>calculate area'
    areas=mod_grid.get_area()
    print 'shape of areas:', npy.shape(areas)
    
    sp=npy.zeros((mod_grid.nlon, mod_grid.nlat), float)
    sp[:,:]=1000.0
    
    
    olat=npy.arange(-90, 90, 30)
    nobs=npy.size(olat)
    olon=npy.zeros(nobs)
    olon[:]=-30.0
    
    print '>set hintpl'
    
    hintpl=mod_grid.get_hintpl(olon, olat)
    
            
     
    
