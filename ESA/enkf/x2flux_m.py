"""
    Class for projecting coefficient to flux
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.11.12
    History: v0.95, 2013.01.07
    
    Classes:
    ===============================================
    1. x2flux_cl: Class for projecting coefficients vector 
    

    Functions:
    ===================================================
    1. set_localization_wgt: 
"""


import ESA.util.otool_gcfile_io as gcfio 
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.surface.bf_file_m as bf_ncio
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.great_circle_distance as gdist
import ESA.util.flux_regrid as fgrd
import ESA.util.gen_plots as gpl

import numpy as npy 
import pylab as plt 



def set_localization_wgt(grd_lon, grd_lat, olon, olat, \
                             wlen, threshold, \
                             by_relative,\
                             **keywords):
    """
    
    Set up weighting function for localization according to distance from 
    observations. 
    
    Usually, localization is necessary when a partial representation is used. 
    
    Inputs:
    ---------------------------------------------
    1. grd_lon:<array, (nlon)>: map longititude
    2. grd_lat:<array, (nlat)>: map latititude
    3. olon, olat:<array, (nlat)>:
    4. wlen:<float>: length 
    5. threshold:<float>: threshold above which wgt will be set to 1.0 (see Note 1)
    6. by_relative:<T/F>: if true, threshold will be interprested as percentage
    7. keyword:<dict>: dictionary extra inputs
    
    Returns:
    -------------------------------------------
    1. wmap:<array, (nlon, nlat)>: wmap 
    

    Notes:
    -----------------------------------------
    1. threshold normally should be less than 1.0
    

    """
    

    # S1: set distance map for points at (0, grd_lat)
    
    dist_table=gdist.set_circle_distance_table(grd_lon, grd_lat)

    # S2: get weighting factor
    
    plt.figure(6)
    gpl.plot_map(dist_table[:,:,10], grd_lon, grd_lat, use_pcolor=1)
    plt.figure(1)
    
    lon_period=360.0
    
    wmap=fgrd.get_obs_wgt_map(grd_lon, grd_lat, olon, olat, dist_table, wlen,lon_period)
    
    
    if (by_relative):
        max_wgt=npy.maxval(wmap.flat)
        wmap=wmap/max_wgt
        
    # S3: set values above threshold to be 1
    
    wmap=npy.where(wmap>=threshold, 1, wmap)
    
    
    return wmap

    
   
def read_bf_to_dict(datapath, \
                        flnm, \
                        yyyy, mm, dd, \
                        varname_lst=['lon', 'lat', 'flux', 'map', 'area'],\
                        **keywords):
    
    """
    
    read basis function (ensemble function from a netcdf file
    
    Inputs:
    

    Inputs:
    --------------------------------------------------------
    1. datapath:<str>: file path for aod data 
    2. flnm: <flnm>: file name for aod file
    3. yyyy: <int>:  year 
    4. mm:<]int>:     month 
    5. dd:<int>:     day 
    6. varname_lst:<list, t:str>: list of vars to be read
    
    Returns:
    ----------------------------------------
    var_dict:<dict>: dictionary for vars read from BF file
    
    """

    # S1: get file name 
    
    doy=tm.day_of_year(yyyy, mm, dd)
    
    ncflnm=ovar.construct_filename(datapath, flnm, \
                                       XYYYYX=yyyy,\
                                       XMMX=mm,\
                                       XDDX=dd,\
                                       XDOYX=doy)
    
    print ncflnm
    
    # S2 read in data
    
    var_lst=ncfio.ncf_read(ncflnm, varnames)
    var_dict={}
    nvars=len(varnames)
    for ivar in range(nvars):
        var_dict.update({varname[ivar]:var_lst[ivar]})
        
    return var_dict

    

  
                

def read_bf(datapath, flnm, yyyy, mm, dd, \
                lon_nm='lon', lat_nm='lat', \
                flux_nm='flux', \
                **keywords):
    """
    
    read basis function (ensemble function from a netcdf file
    Inputs:
    
    
    Inputs:
    --------------------------------------------------------
    1. datapath:<str>: file path for aod data 
    2. flnm: <flnm>: file name for aod file
    3. yyyy: <int>:  year 
    4. mm:<]int>:     month 
    5. dd:<int>:     day 
    6. lon_nm:<str>: name of longitude
    7. lat_nm:<str>: name of latitude
    8. flux_nm:<str>: name of flux
    
    
    """

    # S1: get file name 
    
    doy=tm.day_of_year(yyyy, mm, dd)
    
    ncflnm=ovar.construct_filename(datapath, flnm, \
                                       XYYYYX=yyyy,\
                                       XMMX=mm,\
                                       XDDX=dd,\
                                       XDOYX=doy)
    print ncflnm
    
    # S2 read in data
    varnames=[lon_nm, lat_nm, flux_nm]
    lon, lat, regflux=ncfio.ncf_read(ncflnm, varnames)
    return lon, lat, regflux
    

class x2flux_cl: 
    
    def __init__(self,\
                     datapath, flnm, \
                     yyyy, mm, dd, \
                     bf_get=read_bf, \
                     bfcl=None, \
                     lon_nm='lon', lat_nm='lat', \
                     flux_nm='flux', \
                     fio_keywords={}, \
                     **keywords):
        
        
        
        """ initialization 
        
        Inputs:
        ------------------------------------------------------
        1. datapath:<str>: file path for aod data 
        2. flnm: <flnm>:  file name for aod file
        3. yyyy: <int>:   year 
        4. mm:<]int>:     month 
        5. dd:<int>:      day 
        6. bf_get:<func>: function for reading basis function (or ensemble function)
        7. bfcl:<bf_cl>:  basis function class
        8. lon_nm:<str>: name of longitude
        9. lat_nm:<str>: name of latitude
        10. flux_nm:<str>: name of flux
        11. fio_keywords: <dict>:  extra inputs for file accessing
        12. keywords:<dict>: attribute 
        

        
        """
        
        # S1: Set the values
        
        self.flnm=flnm
        self.datapath=datapath
        
         
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        self.fio_keywords=fio_keywords
        self.bf_get=bf_get
        self.bfcl=bfcl
        self.lon=None
        self.lat=None
        self.flux=None
                
        self.lon_nm=lon_nm
        self.lat_nm=lat_nm
        self.flux_nm=flux_nm
        
        
        # S2: set attributes
        
        self.attr_dict={}
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
        
        
    def load_flux_from_file(self, **keywords):
        """
        Inputs:
        ---------------------------------------------
        1. datapath:<str>: file path for aod data 
        2. flnm: <flnm>:  file name for aod file
        3. keywords:<dict>: extra inputs
        ---reserved keys:
        --->lon_nm:<str>: name of longitude
        --->lat_nm:<str>: name of latitude
        --->flux_nm:<str>: name of flux
        --->flnm:<str>: name of flux file
        --->path:<str>: name of data path
        --->yyyy, mm, dd:<int>: year, month, day 
        
        


        Notes:
        --------------------------------------------
        1. if existing self.bfcl, self.bfck.load_bf will be called. 
        otherwise, self.bf_get will be invoked. 
        
        
        """
        # S1: check the keyword 
        
        if ('lon_nm' in keywords):
        
            self.lon_nm=keyword['lon_nm']
        
        lon_nm=self.lon_nm
        
        if ('lat_nm' in keywords):
            self.lat_nm=keyword['lat_nm']
            
        lat_nm=self.lat_nm
        
        if ('flux_nm' in keywords):
            self.flux_nm=keyword['flux_nm']
        
        flux_nm=self.flux_nm


        
        if ('path' in keywords):
            self.datapath=keywords['path']
            
        datapath=self.datapath
        
        if ('flnm' in keywords):
            self.flnm=keywords['flnm']

        flnm=self.flnm
        
        if ('yyyy' in keywords):
            self.yyyy=keywords['yyyy']

        if ('mm' in keywords):
            self.mm=keywords['mm']

        if ('dd' in keywords):
            self.dd=keywords['dd']

            
        
        
        # S2: read in basis (ensemble) flux function 

        if (self.bfcl==None):
            # T1: read in basis (ensemble) via bf_get function 
            
            lon, lat, flux=self.bf_get(self.datapath, self.flnm, \
                                           self.yyyy, self.mm, self.dd, \
                                           lon_nm=lon_nm, lat_nm=lat_nm, flux_nm=flux_nm,\
                                           **self.fio_keywords)
            
        else:
            # T2: read in basis (ensemble) via bfcl class

            # #c: read 
            self.bfcl.load_bf(self.datapath, self.flnm, self.yyyy, \
                                  self.mm, self.dd,  **self.fio_keywords)
            
            lon=self.bfcl.fdesc.get_data(lon_nm)
            lat=self.bfcl.fdesc.get_data(lat_nm)
            flux=self.bfcl.fdesc.get_data(flux_nm)
                
         
        
        self.lon=lon
        self.lat=lat
        self.flux=flux
        
        
        

    def x2flux(self, x, \
                   flux=None, \
                   locz_wgt=None,\
                   keep_flux=False, **keywords):
        
        """
        Inputs:
        ------------------------------------------
        1. x:<array, (ne, (ne))>: coefficient matrix 
        2. flux:<array, (nlon, nlat, ne)>: flux basis function or ensemble 
        3. locz_wgt:<array, (nlon, nlat)>: localization scaling factor
        4. keep_flux:<T/F>: if true, flux will be saved 
        

        """
        
        
        
        # S1: read flux if necessary     
        
        if (flux==None):
            self.load_flux(**keywords)
            flux=self.flux
        
        # S2: calculate x*flux or flux*wgt*x

        dims=npy.shape(x)
        ndim=npy.size(dims)
        
        if (ndim==1):
            
            x=npy.reshape(x,[-1,1])
            
        if (locz_wgt==None):
            flux_out=fgrd.do_flux_by_coef(flux, x)

        else:
            flux_out=fgrd.do_flux_by_coef_wgt(flux, locz_wgt, x)
        
        # S3: keep flux if required 
        
        if (keep_flux):
            self.flux=flux

        else:
            self.flux=None
        
        flux_out=npy.squeeze(flux_out)
        
        return flux_out
    
    def copy(self, yyyy, mm, dd, **keywords):
        """
        Make a copy of itself, and replace the keywords 
        
        Inputs:
        ---------------------------------------------
        1. yyyy, mm, dd:<int>: date
        2. keywords:<dict>: extra inputs 

        Returns
        -------------------------------------------
        cl_new:<x2flux_cl>: a copy of itself
        

        """
        all_keywords=dict(self.attr_dict)
        all_keywords.update(keywords)
        
        
        cl_new=x2flux_cl(self.datapath, self.flnm, \
                             yyyy, mm, dd, \
                             bf_get=self.bf_get, \
                             bfcl=self.bfcl, \
                             lon_nm=self.lon_nm, lat_nm=self.lat_nm, \
                             flux_nm=self.flux_nm, \
                             fio_keywords=self.fio_keywords, \
                             **all_keywords)
        
        
        return cl_new
    
    
# <<< TEST >>> 

if (__name__=='__main__'):
    
    # file and data

    datapath='/scratch/local/otool_data/surface_flux/'
    flnm='CO2_EMISSION.XYYYYXDXDOYX.nc'
    yyyy, mm, dd=2009, 1,1
    
    # create class
    
    xflx=x2flux_cl(datapath, flnm, \
                       yyyy, mm, dd, \
                       bf_get=read_bf, \
                       bf_cl=None, \
                       lon_nm='longitude', \
                       lat_nm='latitude', \
                       flux_nm='flux')

    
    
    # load flux 

    xflx.load_flux_from_file()
    
    grd_lon=xflx.lon
    grd_lat=xflx.lat
    
    # set observation 
    
    olon=npy.arange(-60, 60, 30.0)
    olat=npy.arange(-60, 60, 30.0)
   
    
    # set localization parameters for 
    
    wlen=1000.0 # 800 km 
    threshold=0.5 #
    by_relative=False
    
    wmap=set_localization_wgt(grd_lon, grd_lat, \
                                  olon, olat, \
                                  wlen, threshold, \
                                  by_relative)
    
    
    flux=xflx.flux
    nx, ny, nz=npy.shape(flux)
    x=npy.ones(nz, float)
    x=x
    
    
    flux2=xflx.x2flux(x, \
                          flux, \
                          locz_wgt=wmap)
    
    
    plt.figure(2)
    gpl.plot_map(flux2, grd_lon, grd_lat, use_pcolor=1)
    
    
    flux3=xflx.x2flux(x, \
                          flux, \
                          locz_wgt=None)
    
    
    print npy.shape(flux)
    
    plt.figure(3)
    gpl.plot_map(flux3, \
                     grd_lon, grd_lat, \
                     use_pcolor=1)
    
    plt.figure(4)
    gpl.plot_map(npy.sum(flux, axis=2), grd_lon, grd_lat, use_pcolor=1)
    
    plt.show()
    
    
