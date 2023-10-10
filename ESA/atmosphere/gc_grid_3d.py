""" Class for grid used by GEOS-Chem CTM.  
   
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.10.12
    
    Classes:
    -------------------------------------------------------
    1. gc_grid_cl: derived class for GEOS-Chem model grid
    
    

"""
import numpy as npy
# repath
# import gp_axis_m as axis_m
# import gp_grid_m as grid_m
#repath
import compute_gc_grid as cmg
import ctm_grid_3d as ctm_grd
import ESA.util.pres_m       as pm
import ESA.util.otool_obj  as oob
import ESA.util.message_m as msm

class gc_grid_cl (ctm_grd.ctm_grid_cl): 
    
    """ derived class for model grid 
    
    Members (overrided and newly added):
    ---------------------------------------
    # window definitition
    
    1. ix: longitude  size 
    2. jx: latitude size  
    3. lx: vertical size 
    4. ifirst: longitude shift 
    5. jfirst: latitude shift 
    6. lfirst: vertical shift
    7. halfpolar: polar configuration/ 
    8. centre180: =1, if longitudes are represented in (0 to 360) degrees. 
    
    #   global H-grid (important for nested model)
    
    9.  glb_ix: global x size
    10. glb_jx: global y size
    11. glb_lx: global z size
    12. glb_lon: global longitude centers
    13. glb_lat: global latitude  centers
    14. glb_lon_edge:  global longitude edges 
    15. glb_lat_edge:  global latitude edges 

    
    # global V-vertical 
    
    21. is_nested: whether it is nested 
    
    22. lonres: longitude resolutions
    23. latres: latitude resolutions 
    
    24. geos_ver: geos-chem version 
    

    # initialization control
    25. resize_use_mod_res
    

    

    
    Functions: (new or overridden functions)

    ---------------------------------------------------------------------
        
    1. __init__: initialize the class
    2. get_lon: return longitude grid 
    3. set_lon: set longitude grid 
    4. get_lat: return latitude grid 
    5. set_lat: set longitude grid 
    
    #
    
    6. compute_mod_pres: compute model pressure and reset the attributes
    7. compute_mod_pres_edge:  compute model pressures at level edges and reset the attributes
    
    # edges 
    8. set_lon_edge: set longitude edge
    9. set_lat_edge: set latitude edge
    10. compute_lat_edge: compute_lat_edge
    11. compute_lon_edge: compute_lon_edge 
    12. get_area: calculate areas of grid boxes 

    # global grid 
    13.get_glb_lon: return global longitude grid 
    14.get_glb_lon_edge: return global longitude grid at box edges
  
    14.get_glb_lat: return global latitude grid 
    15.get_glb_lat_egde: return global latitude grid 
    



    Attributes (Added for resevered keywords):
    -------------------------------------------------------------
    1. glb_mod_pres: model pressure at level centers
    2. glb_mod_pres_edge: model pressures at level edges

    """
    
    
    def __init__(self, ix, jx, lx, \
		     ifirst=1, jfirst=1, lfirst=1,\
		     halfpolar=0, centre180=0, \
                     lonres=None, latres=None,\
                     is_nested=0,\
                     use_reduced=1, \
                     geos_ver=5, \
                     mod_res=None,\
                     resize_use_mod_res=True,\
                     mask_val=oob.fill_val, \
                     **keywords\
                     ):
        
        """ Initialize the class 
        Input:
        ----------------------------------
        1. ix, jx, lx:<nteger> : the sizes of lon, lat, z 
        2. ifirst, jfirst, lfirst:<integer> the start point for the data to store. 
        3. halfpolar, centre180:<integer>: 0 or 1 for the type of lon 
        4. lonres, latres:  <float>: the resolution of the lat and lon 
        
        """
        # S1 dictionary
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_grid})
        
        
        # S2 initialize values
        # T1 size
        self.ix=ix
        self.jx=jx
        self.lx=lx
        
        
        # T2 shift: 1 : No shift
        
        self.ifirst=ifirst
        self.jfirst=jfirst
        self.lfirst=lfirst
        self.halfpolar=halfpolar
        self.centre180=centre180
        
        self.is_nested=is_nested
        
        self.lonres=lonres
        self.latres=latres

        self.geos_ver=geos_ver
        self.use_reduced=use_reduced
        self.mod_res=mod_res
        self.resize_use_mod_res=resize_use_mod_res
        
          
        # S3   global grid (important for nested model)
                
        self.glb_ix=0
        self.glb_jx=0
        self.glb_lx=0

        self.glb_lon=None
        self.glb_lat=None

        self.glb_lon_edge=None
        self.glb_lat_edge=None
        
        # global V-vertical 
        
        
        
        
        ## T1 assign model resolution
        
        if (mod_res==None):
            if (self.is_nested==0):
                self.mod_res=cmg.get_model_resolution(self.ix, self.jx)
            else:
                if (geos_ver==4):
                    self.mod_res='1x1'
                elif (geos_ver==5):
                    self.mod_res='05x0666'
                else:
                    msg='wrong configuration for nested grid'
                    msm.show_err_msg(msg)
            

            
        else: # if the model resolution is given, ix, jx will be overrided
            
            self.mod_res=mod_res
            if ((self.is_nested==0) & (resize_use_mod_res)):
                lons=cmg.get_model_lon(model_res=self.mod_res)
                # call one
        
                lats=cmg.get_model_lat(model_res=self.mod_res)
                
                self.ix=npy.size(lons)
                self.jx=npy.size(lats)
        
                 
             
        
        # S5 initialised  longitude and latitude
        
        # N1 ifirst and jfirst will also been used

        # member lon, glb_lon, and glb_lon_edge will be set 
        
        
        
        lons=self.init_lon(self.ix)
        lats=self.init_lat(self.jx)
        
        # member lat, glb_lat, and glb_lat_edge  will be set 
        
        lats=self.init_lat(self.jx)
        
        #S4 initialized vertical grid 
      
        rlz=self.init_lz(self.lx)
        
        ctm_grd.ctm_grid_cl.__init__(self, lons, lats, rlz, \
                                         hresolution=self.mod_res, \
                                         mask_val=mask_val)
        
        # print self.attr_dict
        # print oob.ot_grid
        
        
        if (self.geos_ver==5):
            ap, bp=pm.get_geos5_ap_bp(self.use_reduced)
        
        if (self.geos_ver==4):
            ap, bp=pm.get_geos4_ap_bp(self.use_reduced)
        
        self.glb_lx=npy.size(ap)-1
        
        self.attr_dict.update({'ap':ap})
        self.attr_dict.update({'bp':bp})
        
        

        
        
        #S4 initializing parent members
        
        
      
        # assign ab bp
        
        
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
    
    
                 
    def init_lat(self, jx=None, jfirst=None, \
                     latres=None, is_nested=None, \
                     resize_use_mod_res=True):
        """ initialize latitude grid 
        Inputs:
        -------------------------------------
        1. jx:<integer>: size of longitude
        2. jfirst:<integer>: poistion of first longitude in global grid
        3. latres:<float>: longitude resolution
        4. is_nested:<integer>: 0==not nested grid 
        5. resize_use_mod_res:<T/F>: if Ture, jx will be calculated from latres when it is not 
        provided.
        

        Returns:
        -------------------------------------
        lats:<array, (ix,)>: longitudes
        
        """
        # 1 check inputs

            
        if (is_nested==None):
            is_nested=self.is_nested
        
        self.is_nested=is_nested
        
        if (jfirst==None):
            jfirst=self.jfirst
        
        self.jfirst=jfirst
                   
        
        # find latres if it is not given
        
        if (latres==None):
            if (is_nested==1):
                if (self.geos_ver==4):
                    latres=180.0/180.0 # 1 
                else:
                    latres=180.0/360.0  # 0.5
                    
            else:
                
                if (jx==None):
                    latres=self.latres
                else:
                    latres=180.0/(jx-1)
                    
        
                
        self.latres=latres
        
        if (jx==None):
            if (resize_use_mod_res):
                jx=int(180.0/latres)+1
            else:
                jx=self.jx
        
        self.jx=jx
        self.glb_jx=jx
        
        ##  re-do global_jx if nested modol is used
        if (is_nested==1):
            self.global_jx=int(180.0/latres)+1
        
        self.glb_lat=cmg.get_model_lat(nlat=self.glb_jx)
        self.glb_lat_edge=cmg.get_model_lat_edge(nlat=self.glb_jx)
        
        lats=cmg.get_wnd_lats(self.glb_lat,  self.jx, self.jfirst)
        
        return lats

    def init_lz(self, lx=None, lfirst=None, geos_ver=None, use_reduced=None, \
                    resize_use_glb=True):
        
        """ initialize vertical levels  
        
        Inputs:
        -------------------------------------
        1. lx:<integer>: size of vertical levels 
        2. lfirst:<integer>: poistion of first level
        3. geos_ver:<integer>: version for GEOS data 
        4. use_reduced:<integer>: 1==reduced vertical will be use. 
        5. resize_use_glb:<T/F>: if ture, lx will be resized to glb_lx when it is not given 
        
        
        

        Returns:
        -------------------------------------
        lz:<array, (lx,)>: longitudes
        """
        
        
        if (geos_ver==None):
            geos_ver=self.geos_ver
        
        self.geos_ver=geos_ver

        if (use_reduced==None):
            use_reduced=self.use_reduced
            
        self.use_reduced=use_reduced
        
        if (lx==None):
            if (resize_use_glb):
                lx=self.glb_lx
            else:
                lx=self.lx
        
        
        self.lx=lx
        
        if (lfirst==None):
            
            lfirst=self.lfirst
        
        self.lfirst=lfirst

        
        
        # print self.lfirst, self.lx
        # print self.glb_lx
        # print self.geos_ver
        
        lz=npy.arange(self.lfirst-1, self.lfirst-1+self.lx)
        return lz
    
    
    def init_lon(self, ix=None, ifirst=None, lonres=None, \
                     centre180=None, is_nested=None, \
                     resize_use_mod_res=True):
        
        """ initialize longitude grid 
    
        Inputs:
        -------------------------------------
        1. ix:<integer>: size of longitude
        2. ifirst:<integer>: poistion of first longitude in global grid
        3. lonres:<float>L: longitude resolution
        4. centre180:<integer>: 0==centerred around 0
        5. is_nested:<integer>: 0==not using nested grid 
        6. resize_use_mod_res:<T/F>: True, ix will be calculated from lonres, if it is not 
        provided. 
        
 
        Returns:
        ------------------------------------------
        1. lons:<array, (ix,)>: longitudes
        
        -------------------------------------
        """
        #S1 check inputs, if they are none, they will be filled by self values
        

        if (is_nested==None):
            is_nested=self.is_nested
        
        self.is_nested=is_nested
        
        if (ifirst==None):
            ifirst=self.ifirst
        
            
        self.ifirst=ifirst
        
        
        ##T1 lonres will be calculated, if it is not provided 

        if (lonres==None):
            if (is_nested==1): # if it is a nested grid 
                if (self.geos_ver==4):
                    lonres=1.0
                else:
                    lonres=360.0/540

            else:  # get it from ix if it is not a nested grid 
                if (ix==None):
                    lonres=self.lonres
                else:
                    lonres=360.0/ix
        
        self.lonres=lonres
        
        if (ix==None):
            if (resize_use_mod_res):
                ix=int(360.0/lonres)
            else:
                ix=self.ix
        
        self.lonres=lonres
        self.ix=ix
        
        ##T2 global longitude size will be reset if the nested grid is used
        self.glb_ix=ix
        if (is_nested==1):
            self.glb_ix=int(360.0/self.lonres)            
            
        if (centre180==None):
            centre180=self.centre180
            
        self.centre180=centre180
        
        #S2 calculate longitude grid 
        self.glb_lon=cmg.get_model_lon(nlon=self.glb_ix)
        self.glb_lon_edge=cmg.get_model_lon_edge(nlon=self.glb_ix)
        
        #S3 check lons are really in use. 
        
        lons=cmg.get_wnd_lons(self.glb_lon, self.ix, self.ifirst, self.centre180)
        
        return lons
    
    
        

    def set_lon(self, ix=None, ifirst=None, lonres=None, centre180=None, \
                    is_nested=None):

        
        
        
        """ Re-set longitude  grid (over-riding function)
 
        Inputs:
        -------------------------------------
        1. ix:<integer>: size of longitude
        2. ifirst:<integer>: poistion of first longitude in global grid
        3. lonres:<float>L: longitude resolution
        4. centre180:<integer>: 0==centerred around 0
        5. is_nested:<integer>: 0==not using nested grid 
        
        
        Returns:
        --------------------------------------
        1. lons:<array, (ix,)>: longitudes
        
        -------------------------------------
        """
        
        lons=self.init_lon(ix=ix, ifirst=ifirst, lonres=lonres, centre180=centre180, \
                               is_nested=is_nested)
        
        # updating parent
        
        ctm_grd.ctm_grid_cl.set_lon(self, lons)
        return lons
    
        
    def set_lat(self, jx=None, jfirst=None, latres=None, is_nested=None):
        
        """ set latitude  grid (over-riding function 
        
        Inputs:
        -------------------------------------
        1. jx:<integer>: size of latitude
        2. jfirst:<integer>: poistion of first latitude in global grid
        3. latres:<float>: latitude resolution
        
 
        Returns:
        -----------------------------------------
        1. lats:<array, (ix,)>: longitudes
        ------------------------------------
        """

        lats=self.init_lat(jx=jx, jfirst=jfirst, latres=latres, is_nested=is_nested)
        
        # updating parent
        ctm_grd.ctm_grid_cl.set_lat(self, lats)
        
        return lats
    
    def get_mod_pres_edge(self):
    
        """ 
        
        Retrieve model pressure at box edges from dictionary 
        
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
    
    def get_mod_pres(self):
        
        """ 
        Retrieve model pressure from dictionary 
        
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
    
    def compute_mod_pres(self,sp):
        """
        set pressure at model level centres 
        Inputs:
        -------------------------------------
        1. sp: <array, (ix,jx)> surface pressure
        
        Returns:
        ----------------------------------
        1. mod_res: <array, (ix, jx, lx)>: pressure at model level
        """
        ap=self.attr_dict['ap']
        bp=self.attr_dict['bp']
        nlon, nlat=npy.shape(sp)
        
        if ((nlon==self.glb_ix) & (nlat==self.glb_jx)):
            # sp is given as global ones
            glb_mod_pres=pm.get_mod_pres(sp, ap, bp)
            self.attr_dict.update({'glb_mod_pres':glb_mod_pres})
            self.attr_dict.update({'glb_surf_pres':sp})
            
            
            mod_pres=glb_mod_pres[self.ifirst-1:self.ifirst+self.ix-1, \
                                      self.jfirst-1:self.jfirst+self.jx-1,self.lfirst-1:self.lx+self.lfirst-1]
            wnd_sp=sp[self.ifirst-1:self.ifirst+self.ix-1, \
                                      self.jfirst-1:self.jfirst+self.jx-1]
            
            self.attr_dict.update({'surf_pres':wnd_sp})
            
        else:
            
            mod_pres=pm.get_mod_pres(sp, ap, bp)
            mod_pres=mod_pres[:,:,self.lfirst-1:self.lx+self.lfirst-1]
            # only window values given 
            
            self.attr_dict.update({'surf_pres':sp})
            
        
        self.attr_dict.update({'mod_pres':glb_mod_pres})
            
        return  mod_pres
    
    def compute_mod_pres_edge(self,sp):
        """
        Set model pressure at model level edges
        
        Inputs:
        -------------------------------------
        1. sp: <array, (ix,jx)> surface pressure
        
        Returns:
        1. mod_res: <array, (ix, jx, lx)>: pressure at model level
        """

        ap=self.attr_dict['ap']
        bp=self.attr_dict['bp']
        nlon, nlat=npy.shape(sp)
        
        if ((nlon==self.glb_ix) & (nlat==self.glb_jx)):
        
            # sp is given as global ones
            glb_mod_pres_edge=pm.get_mod_pres_edge(sp, ap, bp)
            self.attr_dict.update({'glb_mod_pres_edge':glb_mod_pres_edge})
            self.attr_dict.update({'surf_pres':sp})
            
            mod_pres_edge=glb_mod_pres_edge[self.ifirst-1:self.ifirst+self.ix-1,\
                                                self.jfirst-1:self.jfirst+self.jx-1,\
                                                self.lfirst-1:self.lx+self.lfirst]
        else:
            
            self.attr_dict.update({'surf_pres':sp})
            mod_pres_edge=pm.get_mod_pres_edge(sp, ap, bp)
            mod_pres_edge=mod_pres_edge[:,:,self.lfirst-1:self.lx+self.lfirst]
            
        
        self.attr_dict.update({'mod_pres_edge':mod_pres_edge})
            
        
        return  mod_pres_edge

    
    def get_lon_edge(self):
        """
        Returns:
        ----------------------------------------
        1. lons:<array, (ix+1,)>: longitude at grid edges
        
        """
        
         
        if ('lon_edge' in self.attr_dict):
            lons=self.attr_dict['lon_edge']
         
        else:
            lons=cmg.get_wnd_lons(self.glb_lon_edge, \
                                      self.ix+1, self.ifirst, self.centre180)
            self.set_lon_edge(lons)
        
        return lons

    def compute_lon_edge(self):
        
        """
        Returns:
        1. lons:<array, (ix+1,)>: longitude at grid edges
        
        """
        
        lons=cmg.get_wnd_lons(self.glb_lon_edge, \
                                  self.ix+1, self.ifirst, self.centre180)
        self.set_lon_edge(lons)
        
        return lons

     

    def get_lat_edge(self):
        """
        Returns:
        ------------------------------------------
        1. lats:<array, (jx+1,)>: longitude at grid edges
        
        """
        if ('lat_edge' in self.attr_dict):
            lats=self.attr_dict['lat_edge']
        else:
            lats=cmg.get_wnd_lats(self.glb_lat_edge, self.jx+1, self.jfirst)
            self.set_lat_edge(lats)
        
        return lats

    
    def compute_lat_edge(self):
        
        """
        Returns:
        -----------------------------------------
        1. lats:<array, (jx+1,)>: longitude at grid edges
        
        """
        lats=cmg.get_wnd_lats(self.glb_lat_edge, self.jx+1, self.jfirst)
        self.set_lat_edge(lats)
        
        return lats

    
    def get_area(self):
        """
        calculate area for grid box
        Returns:
        1. area:<array, (nx, ny)>: areas for grid boxes in m^2
        
        """

        lons=self.get_lon_edge()
        lats=self.get_lat_edge()
        
        area=cmg.get_area(lons, lats)
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
        
        # make a copy
        
        new_grid=gc_grid_cl(self.ix, self.jx, self.lx, ifirst=self.ifirst, \
                                jfirst=self.jfirst, lfirst=self.lfirst,\
                                halfpolar=self.halfpolar, \
                                centre180=self.centre180, \
                                lonres=self.lonres, latres=self.latres,\
                                is_nested=self.is_nested,\
                                use_reduced=self.use_reduced, \
                                geos_ver=self.geos_ver, \
                                mod_res=self.mod_res,\
                                resize_use_mod_res=self.resize_use_mod_res,\
                                mask_val=self.mask_val, **attrs)
        return new_grid
    
    def get_glb_lon(self):
        """
        get global longitude grid
        Returns
        -----------------------------------------

        1. glb_lon:<array, (glb_ix)>: global longitudes

        """
        if (self.glb_lon==None):
            return None
        else:
            # one copy
            
            return npy.array(self.glb_lon)
        

    def get_glb_lon_edge(self):
    
        """
        get global longitude grid
        Returns
        -----------------------------------------

        1. glb_lon_edge:<array, (glb_ix+1)>: global longitudes

        """
        
        if (self.glb_lon_edge==None):
            return None
        else:
            return npy.array(self.glb_lon_edge)
        


    def get_glb_lat(self):
        
        """
        get global latitude grid
        Returns
        -----------------------------------------

        1. glb_lat:<array, (glb_jx)>: global latitudes
        
        """
        
        if (self.glb_lat==None):
            return None
        else:
            # one copy
            
            return npy.array(self.glb_lat)
        

    def get_glb_lat_edge(self):
    
        """
        get global latitude grid
        
        Returns
        -----------------------------------------
        1. glb_lat_edge:<array, (glb_jx+1)>: global latitudes

        """
        
        if (self.glb_lat_edge==None):
            return None
        else:
            return npy.array(self.glb_lat_edge)
    


if (__name__=='__main__'):
    
    print '1. form a grid 4x5 with 47 levels:'
    cm=gc_grid_cl(0, 0, 47, mod_res='4x5')
    nx=cm.nlon
    ny=cm.nlat
    print '===a. longitude' 
    lon=cm.get_lon()
    lon_edge=cm.get_lon_edge()
    print '--->size:', npy.size(lon)
    print '--->first, last:', lon[0], lon[-1]
    print '--->edge size:', npy.size(lon_edge)
    print '--->edge first, last:', lon_edge[0], lon_edge[-1]
    print ''
    
    
    print '=====b. latitude' 
    lat=cm.get_lat()
    lat_edge=cm.get_lat_edge()
    print '--->size:', npy.size(lat)
    print '--->first, last:', lat[0], lat[-1]
    print '--->edge size:', npy.size(lat_edge)
    print '--- edge first, last:', lat_edge[0], lat_edge[-1]
    print ''
    
    
    
    print '2. set model level presssure'
    ps=npy.zeros([nx, ny], float)
    ps[:,:]=1000.0
    mod_pres=cm.compute_mod_pres(ps)
    print npy.shape(mod_pres)
    print mod_pres[0,0,:]
    
    print '3. set model level edge presssure'
    
    mod_pres=cm.compute_mod_pres_edge(ps)
    print npy.shape(mod_pres)
    print mod_pres[0,0,:]
    
    print '4. get_area'
    areas=cm.get_area()
    
    print npy.shape(areas)

    new_cm=cm.copy()
    lons=new_cm.get_lon()
    
    print lons[:]
    
    dims=new_cm.shape()
    print dims

    print new_cm.size('lon')
    
    # global grid 

    glb_lon=new_cm.get_glb_lon()
    glb_lon_edge=new_cm.get_glb_lon_edge()


    glb_lat=new_cm.get_glb_lat()
    glb_lat_edge=new_cm.get_glb_lat_edge()
    

    print npy.size(glb_lon)
    print npy.size(glb_lon_edge)
    
    
    print npy.size(glb_lat)
    print npy.size(glb_lat_edge)
    
    
    
    
    
    
            
     
