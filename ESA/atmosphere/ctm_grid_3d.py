""" Class for 3D (lon, lat, lz) grid used by CTM models.  

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.10.26

    
    Classes:
    ---------------------------------------------------
    1. ctm_grid_cl:  derived class for 3D model grid.
    
"""
import numpy as npy
# repath
# import gp_axis_m as axis_m
# import gp_grid_m as grid_m
# import ESA.util.flib as flb

#repath

import ESA.util.gp_axis_m as axis_m
import ESA.util.horizontal_interp_m as hintpl_m
import ESA.util.vertical_interp_m as vintpl_m
import ESA.util.gp_grid_m as grid_m
import ESA.util.gp_axis_m as axis_m
import ESA.util.geo_constant as gc
import ESA.util.otool_obj    as oob
import ESA.util.pres_m       as pm
import compute_gc_grid       as cmg


class ctm_grid_cl (grid_m.gp_grid_cl): 
    
    """ derived class for model grid 
    It has grid of three axis named as 
    ['lon', 'lat', 'lz']
    
    
    Members:
    ---------------------------------------
    Added or over-ridden members
    
    1. attr_dict:<dict>: dictionary for attributes
    2. vtype: Type of the vertical grid
    3. htype: Type of horizontal grid
    4. hresolution: resolution of the horizontal grid
    5. nlon:<integer>: size of longitude grid
    6. nlat:<integer>: size of latitude grid
    7. nz:<integer>: size of vertical grid
    8. mask_val:<float>: filling for bad or missing data
    
    
    Functions (overridden and the new ones):
    -------------------------------------------
    1. __init__: initialize the class
    2. set_mask: set mask values
    
    # horizontal 
    
    3. get_lon: return longitude grid 
    4. set_lon: set longitude grid 
    5. get_lat: return latitude grid 
    6. set_lat: set longitude grid 
    
    #
    7. set_lz: set vertical level 
    8. set_mod_pres: set model pressure
    9.  compute_mod_pres: compute model pressure and reset the attributes
    10. compute_mod_pres_edge:  compute model pressures at level edges and reset the attributes
    11. get_mod_pres: get model pressure 
    12. get_mod_pres_edge: get model pressure at level edges
    
    
    # edges 
    15. set_lon_edge: set longitude edge
    16. set_lat_edge: set latitude edge
    17. compute_lat_edge: compute_lat_edge
    18. compute_lon_edge: compute_lon_edge 
    19. get_area: calculate areas of grid boxes 
    
    # interploation
    20. get_hintpl: construct class for horizontal interpolation
    21. get_vintpl: construct class for vertical interpolation
    
    
    
    Attributes (reserved):
    -----------------------------------------------------
    
    # grid edge 
    1. lon_edge:<array>: longitude grid 
    2. lat_edge:<array>: latitude grid
    
    # pressure 
    3. surf_pres :<array, (nlon, nlat)>: surface pressure
    4. ap :<array,(nz+1)>: pressure coefficients
    5. bp:<array, (nz+1)>: pressure coefficients 
    6. mod_pres:<array, (nlon, nlat, nz)>: model pressure in hPa 
    7. mod_pres_edge:<array, (nlon, nlat, nz+1))>: pressure at model level edges 
    
    # filling values
    8. mask_val:<float>: filling for bad or missing values. 
    
    1. ap:  ap coefficients for pressure of hybrid grids
    2. bp:  bp coefficients for pressure of hybrid grids
    3. surf_pres: surface pressures in hPa
    
        
    """
    
    
    def __init__(self, lon, lat, lz, \
                     vtype='Hybrid', \
                     htype='Gaussian', \
                     hresolution='4x5', \
                     mask_val=oob.fill_val, \
                     **keywords):
        
        """initialize grid 
        Inputs:
        ---------------------------------------------------
        1. lon, lat, lz:<array>: array for model grids.
        ---their sizes are nlon, nlat, nz, respectively

        2. vtype:<string>: type of the vertical grid 
        
        3. htype:<string>: type of the horizontal grid
        4. hresolution:<string>: type of the horizontal grid 
        
        5.  mask_val=oob.fill_val,

        6. keywords: <dict>: attributes 
        
        Reserved words:
        
        --->lon_edgee:<array, (nlon+1,)>: array of model lon edges
        --->lat_edgee:<array, (nlat+1,)>: array of model lon edges
        --->ap, bp:<array, (nz+1,)>: coefficients for pressure level
        --->surf_pres:<array, (nlon, nlat)>: pressures at surface
        --->mod_pres:<array, 
        
        """
        # S1: construct gp axis 
        
        
        ax_lon=axis_m.gp_axis_cl('lon', lon)
        ax_lat=axis_m.gp_axis_cl('lat', lat)
        ax_lz=axis_m.gp_axis_cl('lz', lz)
        
        self.vtype=vtype
        self.htype=htype

        self.hresolution=hresolution
        
        # S2:  initialize parent class 
        grid_m.gp_grid_cl.__init__(self, [ax_lon, ax_lat, ax_lz])
        
        self.mask_val=-999.0
        self.set_mask(mask_val)

        # self.lon=npy.array(lon)
        # self.lat=npy.array(lat)
        # self.lz=npy.array(lz)
        
        self.nlon=npy.size(lon)
        self.nlat=npy.size(lat)
        self.nz=npy.size(lz)
        
       # S3: set attrib 

        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
            
        
        
    
    def set_mask(self, mask):
        """ set filling value for missing or bad data
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
            
            
    
    def set_lz(self, lz):

        """ set latitude  grid 
        Inputs:
        ------------------------------------
        1. lz:<array>: vertical coordinate
        returns: 
        ------------------------------------
        """
        
        # self.lz=npy.array(lz)
        self.nz=npy.size(self.lz)
        
        ax_lz=axis_m.gp_axis('lz', lz)
        self.update_aixs(ax_lz)
        
    
    def get_lz(self):
        """ retrieve lz (vertical) grid 
        returns: 
        --------------------------
        1. self.lz:<array>: longitudes 
        """
        ax_lz=self.get_axis('lz')
        lz=ax_lz[:]
        return lz

    
   
    def set_mod_pres(self, mod_pres):
        self.attr_dict.update({'mod_pres':mod_pres})
        
                
    def compute_mod_pres(self,sp, \
                             ap, bp, \
                             mod_pres_f=pm.get_mod_pres):
        
        """
        Compute model pressure
        
        Inputs: 
        -----------------------------------------------------------
        1. sp:<array, (nlon, nlat)>: surface pressure in hPa
        2. ap:<array, (nz+1)>: coefficient for model pressure level
        3. bp:<array, (nz+1)>: coefficient for model pressure levels

        Returns:
        1.mod_pres:<nlon, nlat, nz>: pressure in hPa

        
        Notes:
        ----------------------------------------
        1. In default routine, pressure at level edges are 
        caculated by 
        
        pres_edge=ap+bp*sp
                
        
        """
        
        mod_pres=mod_pres_f(sp, ap, bp)
        # update attributes 
        self.attr_dict.update({'mod_pres':mod_pres})
        self.attr_dict.update({'ap':ap})
        self.attr_dict.update({'bp':bp})
        self.attr_dict.update({'surf_pres':sp})
        
        return  mod_pres
    
    def compute_mod_pres_edge(self,sp, ap, bp, \
                                   mod_pres_edge_f=pm.get_mod_pres_edge):
        
        """

        Compute model pressures at level edges 
        
        
        Inputs: 
        -----------------------------------------------------------
        1. sp:<array, (nlon, nlat)>: surface pressure in hPa
        2. ap:<array, (nz+1)>: coefficient for model pressure level
        3. bp:<array, (nz+1)>: coefficient for model pressure levels

        Returns:
        1.mod_pres_edge:<nlon, nlat, nz+1>: pressures in hPa

        
        Notes:
        ----------------------------------------
        1. In this routine, pressures at level edges are 
        caculated by 
        
        pres_edge=ap+bp*sp
                
        
        """
        
        mod_pres_edge=mod_pres_edge_f(sp, ap, bp)
        self.attr_dict.update({'mod_pres_edge':mod_pres_edge})
        self.attr_dict.update({'ap':ap})
        self.attr_dict.update({'bp':bp})
        self.attr_dict.update({'surf_pres':sp})
        
        return  mod_pres_edge
    
    def get_mod_pres(self):
        """ 
        retrieve model pressure from dictionary 
        Returns
        -------------------------------
        
        1. mod_pres:<array, (nz,)>:pressure at model edge 

        """
        
        if ('mod_res' in self.attr_dict):
            mod_pres=self.attr_dict['mod_pres']
            return mod_pres
        else:
            msg="model pressures have not been calculated"
            msm.show_err_msg(msg)
            return None
        
    
    def get_mod_pres_edge(self):
        """ 
        retrieve model pressure at box edges from dictionary 
        
        Returns:
        ----------------------------------
        mod_pres_edge:<array, (nz+1,)>:pressure at model edge 

        """
        
        if ('mod_res_edge' in self.attr_dict):
            mod_pres_edge=self.attr_dict['mod_pres_edge']
            return mod_pres_edge
        else:
            msg="model pressure has not been calculated"
            msm.show_err_msg(msg)
            return None
        
    
        
    def get_regrid_mtx(self, newlon, newlat):
    
        """ get weightings and locations for interpolating fields along a track 
        
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
        hintpl:<interp_cl>: interpolation matrix along the track 

        """
        
        # S1 get axis 
        
        ax_lon=self.get_axis('lon')
        
        ax_lat=self.get_axis('lat')

        # S2 construct the interpolation matrix 
        hintpl=hintpl_m.hinterp_cl(ax_lon, ax_lat, mask_val=self.mask_val)
        
        # S3 set coefficients for interpolation 
        
        hintpl.init_interp(olon, olat)
        
        return hintpl
    

    
    def get_vintpl(self,mod_pres, opres, \
                       is_in_log, \
                       do_reverse, \
                       do_ob_reverse):

        """ get weightings and locations for vertically interpolating fields along a track 
        
        Inputs:
        ------------------------------------
        1. mod_pres:<array, (nob, nz, ..)>: model pressure along the track 
        2. opres:<array, (nob, ob_nz, ...)>: ob pressures to be interpolated to. 
        3. is_in_log:<T/F>: True, if mod_pres and  and opress are in log10
        4. do_reverse:<T/F>: True, if mod_pres is given in decending order 
        5. do_ob_reverse:<T/F>: True,  if opress is given in decending order 
        
        
        Returns: 
        ------------------------------------
        hintpl:<interp_cl>: vertical interpolation matrix  

        """
        
        vintpl=vintpl_m.vintpl_cl(mod_pres, \
                                      is_in_log, \
                                      do_reverse, \
                                      mask_val=self.mask_val)
        
        vintpl.init_interp(opres, is_in_log, do_ob_reverse)
        
        return vintpl



    
    
    
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
        #1 get grid egdes 
        rlon=self.compute_lon_edge(lon_st=lon_st, lon_end=lon_end)
        rlat=self.compute_lat_edge(lat_st=lat_st, lat_end=lat_end)
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
        
        new_grid=ctm_grid_cl(self.lon, self.lat, self.lz, vtype=self.vtype, \
                                 htype=self.htype, hresolution=self.hresolution,\
                                 mask_val=self.mask_val, **keywords)
        return new_grid
    
    
            
if (__name__=='__main__'):
     
    print '>1. form a grid 4(lon)x5(lat) grid with 47 levels'
    
    rlon=npy.arange(-180, 180.0, 5)
    axis_lon=axis_m.gp_axis_cl('lon', rlon, ax_unit='Deg')
    
    
    rlat=npy.arange(-90, 90.1, 4.0)
    axis_lat=axis_m.gp_axis_cl('lat', rlat, ax_unit='Deg')
    z=npy.arange(47)
    axis_z=axis_m.gp_axis_cl('z', z, ax_unit='None')
    

    axis_list=[axis_lon, axis_lat, axis_z]

    mod_grid= ctm_grid_cl(rlon, rlat, z)
    
    lon_edge=mod_grid.compute_lon_edge()
    lon=mod_grid.get_lon()
    
    print 'lon[0,-1]:', rlon[0], rlon[-1]
    
    print 'lon_edge[0,-1]:', lon_edge[0], lon_edge[-1]
    
    
    lat_edge=mod_grid.compute_lat_edge()
    
    print 'rlat[0,-1]:', rlat[0], rlat[-1]
    print 'lon_edge[0,-1]:', lat_edge[0], lat_edge[-1]
    
    print '>2. get area'
    
    areas=mod_grid.get_area()
    print 'shape of area:', npy.shape(areas)
    
    sp=npy.zeros((mod_grid.nlon, mod_grid.nlat), float)
    sp[:,:]=1000.0
    
    ap, bp=pm.get_geos5_ap_bp()
    
    mod_pres=mod_grid.compute_mod_pres(sp, ap, bp)
    print '>3. shape mod_pres', npy.shape(mod_pres)
    
    mod_pres_edge=mod_grid.compute_mod_pres_edge(sp, ap, bp)
    
    print '>4. shape mod_pres_edge', npy.shape(mod_pres_edge)
    
    # test interpolation 

    print '>5 test interpolation' 
    
    print 'a. test horizontal interpolation'

    olat=npy.arange(-90, 90, 30)
    nobs=npy.size(olat)
    olon=npy.zeros(nobs)
    olon[:]=-30.0
    opres_1d=npy.arange(100, 1000, 20)
    noz=npy.size(opres_1d)
    opres=npy.zeros([nobs, noz], float)
    opres[:,:]=opres_1d[npy.newaxis, :]
    opres[:,-1]=-999.0
    
    
    hintpl=mod_grid.get_hintpl(olon, olat)
    mod_pres_prof=hintpl.get_profile(mod_pres)
    print '  olon:', olon 
    print '  olat:', olat
    print '  model pressure profile [0,:]',  mod_pres_prof[0,:]
    
    
    print 'b. test veritcal interpolation'

    print '  shape--mode_pres_prof', npy.shape(mod_pres_prof)
    print '  shape-opres', npy.shape(opres)
    
    lg_opres=npy.array(opres)
    lg_mpres=npy.log10(mod_pres_prof)
    
    lg_opres[:,0:-1]=npy.log10(opres[:,0:-1])
    lg_mpres=npy.log10(mod_pres_prof)
    
    vintpl=mod_grid.get_vintpl(lg_mpres, lg_opres, True, True, False)
    
    
    print '  vpl ob:', vintpl.ob_vpl[0,:]
    print '  vpr ob:', vintpl.ob_vpr[0,:]
    print '  colwgt_ob:', vintpl.ob_colwgt[0,:]
    print '  sum colwgt:', sum(vintpl.ob_colwgt[0,:])
    
    print '  pres and lg_pres  for first profile'
    
    print 'grd pres[0, 0:2] ', vintpl.grd_pres[0, 0:2]
    print 'lg(grd)0, 0:2]:', vintpl.lg_grd_pres[0, 0:2]
    
    print 'opres[0, 0:2]:', vintpl.opres[0, 0:2]
    print 'lg(opres)[0, 0:2]:', vintpl.lg_opres[0, 0:2]
    
    
    
    
    mod_ob_prof=vintpl.interpolate_mod_prof(lg_mpres)
    print '  interoplate profile'
    print npy.shape(mod_ob_prof)
    print ' mod_ob_pres[0,:]:', mod_ob_prof[0,:]
    print ' ob_pre[0,:]s:', lg_opres[0, :]
    
    
    ob_mod_prof=vintpl.interpolate_ob_prof(lg_opres)
    
    print npy.shape(ob_mod_prof)
    print  '  ob_mod_prof[0,:]:', ob_mod_prof[0,:]
    print  '  mod pres[0,:]:', lg_mpres[0,:]
    
    
    prof_ob_ones=npy.ones(npy.shape(lg_opres), float)
    prof_md_ones=npy.ones(npy.shape(lg_mpres), float)
    print '>7. calculate column at model grid' 

    col_gp=vintpl.get_mod_column(prof_md_ones)
    
    print ' col_gp:', col_gp
    
    col_gp=vintpl.get_ob_column(prof_ob_ones)
    
    print '>8. calculate column at ob grid' 

    print '  col_op:', col_gp
    print '  reset ob profile values to 0.3 and 0.4'

    
    prof_md_ones[:]=0.3
    print '  shape', npy.shape(prof_md_ones)
    col_gp=vintpl.get_mod_column(prof_md_ones)
    print ' col_gp 1', col_gp
    
    
    print npy.shape(prof_ob_ones)
    prof_ob_ones[:,:]=0.3
    prof_ob_ones[3,:]=0.4
    
    print npy.shape(prof_ob_ones)
    
    col_gp=vintpl.get_ob_column(prof_ob_ones)
    print ' col_gp 2', col_gp
    
    
            
     
    
