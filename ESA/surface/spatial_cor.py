"""
Functions for generating error correlation matrix

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.03.01
    
    
    Functionso:
    ===============================================
    1. construct_region_spatial_err_cov: Generate  error correlations between regions
    
    

"""

import numpy as npy
import ESA.util.great_circle_distance as gcd
import pylab as plb
import numpy.linalg as nlg
import ESA.util.otool_ncfile_io as nfio
import ESA.util.gen_plots as gpl
    

def construct_grid_spatial_err_cov(lon, lat, \
                                       err, \
                                       cor_length, \
                                       scaling_factor, \
                                       dist_threshold=4):
    
    """ 
    
    Generate  error correlations between select_grid boxes
    
    Inputs:
    ------------------------------------------------------
    1. olon:<array, (npt,)>: longitudes for selected points
    2. olat:<array, (npt,)>: latitudes for selected points
    3. err:<arry, (npt,)>: uncerntity 
    4. cor_len:<float>: correlation length 
    5. scaling_factor:<float>: scaling factor for off-diagonl terms
    6. threshold:<float>: threshold for ratio distance/cor_length where correlations to be ignored
    
    Return:
    ---------------------------------------------------------
    1. err_cov:<array, (npt, npt)>: error spatial correlation matrix 
    
    """

    # S1:   initializing 

    
    npt=npy.size(lon)
    err_cov=npy.zeros([npt, npt], float)
    # #c: covert err from scalar to matrix 
    
    err_wy=npy.diag(err) 
    
    for i in range(npt):
        err_cov[i, i]=1.0
    
    
    # S1: calculate distance

    dist=gcd.get_circle_distance_xy(lon, lat, lon, lat)

    # S2: calculate schur_factor 
    schur_factor=npy.zeros(npy.shape(dist), float)
    dist=dist/cor_length
    
    usd_idx=npy.where(dist<dist_threshold) 
    
    schur_factor[usd_idx]=npy.exp(-dist[usd_idx])
   
    # S3: calculate error correlation matrix 

    
    # err_cor= err*schur_factor*err
    
    
    schur_factor=npy.dot(schur_factor, err_wy)
    # print 'shape--schur_factor 1', shape(schur_factor)
    
    err_cov=npy.dot(err_wy, schur_factor)
            
    return err_cov


    

def construct_region_spatial_err_cov(lon, lat, \
                                         reg_idx_lst, \
                                         reg_err_lst, \
                                         cor_length, \
                                         scaling_factor, \
                                         dist_threshold=4):
    
    """ 
    
    Generate  error correlations between regions
    
    Inputs:
    ------------------------------------------------------
    1. lon:<array, (nlon,)>: longitudes
    2. lat:<array, (nlat,)>: latitudes
    3. reg_idx_lst:<list>: list of cells in each regions
    4. reg_err_lst:<list>: list of flux error in each regions
    5. cor_len:<float>: correlation length 
    6. scaling_factor:<float>: scaling factor for off-diagonl terms
    7. threshold:<float>: threshold for ratio distance/cor_length where correlations to be ignored
    
    Return:
    ---------------------------------------------------------
    1. err_cov:<array>: error spatial correlation matrix 
    
    """

    #  S1: get the distance
    lon_m, lat_m=plb.meshgrid(lon, lat)
    lon_m=npy.transpose(lon_m)
    lat_m=npy.transpose(lat_m)
    
    nreg=len(reg_idx_lst)
    err_cov=npy.zeros([nreg, nreg], float)
    
    for i in range(nreg):
        # ##c: index for each cell in region i
        px=reg_idx_lst[i]
        # ##c: flux error at each cell in region i
        wx=reg_err_lst[i]
        
        lonx=lon_m[px]
        latx=lat_m[px]

        err_cov[i, i]=1.0
            
        
        for j in range(i+1, nreg):
            # ##c: index for each cell in region j
            py=reg_idx_lst[j]
            # ##c: flux error at each cell in region j
            wy=reg_err_lst[j]
            
            lony=lon_m[py]
            laty=lat_m[py]
    
	    # print shape(laty), shape(lony), shape(lonx), shape(latx)
            
            # distance (km) between cells in range i and j
            
	    dist=gcd.get_circle_distance_xy(lonx, latx, lony, laty)
            # print 'max--dist', max(dist.flat)
            # tttt=raw_input()
            
            schur_factor=npy.zeros(npy.shape(dist), float)
            dist=dist/cor_length
            
            usd_idx=npy.where(dist<dist_threshold) # only close points are used
            
            schur_factor[usd_idx]=npy.exp(-dist[usd_idx])
            # corrected_error=wx*schur_factor*wy

            schur_factor=npy.dot(schur_factor, wy)
            # print 'shape--schur_factor 1', shape(schur_factor)
            
            schur_factor=npy.dot(wx, schur_factor)
            # print schur_factor
            
            # print  'shape--schur_factor 2',shape(schur_factor)
            
            errx=npy.sum(wx)
            erry=npy.sum(wy)
            # print errx, erry
            
            
            if (abs(errx*erry)>0):
                err_cov[i, j]=scaling_factor*sqrt(abs(schur_factor/(errx*erry)))
                
                err_cov[j, i]=err_cov[i, j]
    
    
    return err_cov



def get_mean_location(lon, lat, flux, area):
    """
    get flux-mean location of regions

    Inputs:
    ---------------------------------------------------
    
    1. lon:<array, (nlon,)>: longitudes
    2. lat:<array, (nlat,)>: latitudes
    3. flux:<array, (nlon, nlat, nlayer)>: region flux
    4. area:<array, (nlon, nlat)>: area
    
    Returns
    
    ----------------------------------------------------
    1. clon:<array, (nz)>: mean longitude
    2. clat:<array, (nz)>: mean latitude
    
    """
    
    clon=list()
    clat=list()
    nx, ny, nz=npy.shape(flux)
    
    for iz in range(nz):
        
        cflux=flux[:,:,iz]
        cflux=area*cflux
        cx_flux=npy.sum(abs(cflux), axis=1)
        cy_flux=npy.sum(abs(cflux), axis=0)
        total_sum=npy.sum(cx_flux)
        
        x=npy.sum(cx_flux*lon)/total_sum
        y=npy.sum(cy_flux*lat)/total_sum
        
        clon.append(x)
        clat.append(y)

    clon=npy.array(clon)
    clat=npy.array(clat)
    
    return clon, clat

    
        
        
def construct_region_spatial_cor_factor(lon, lat, \
                                            flux,\
                                            area,\
                                            cor_length, \
                                            dist_threshold=4):
    
    """ 
    
    Generate  factor for spatial error correlations between regions
    
    Inputs:
    ------------------------------------------------------
    1. lon:<array, (nlon,)>: longitudes
    2. lat:<array, (nlat,)>: latitudes
    3. flux:<array, (nlon, nlat, nz)>: region flux
    4. area:<array, (nlon, nlat)>: area
    5. cor_len:<float>: correlation length 
    6. threshold:<float>: threshold for ratio distance/cor_length where correlations to be ignored
    
    Return:
    ---------------------------------------------------------
    1. cor_factor:<array, (nz, nz)>: spatial correlation factor
    
    """

    # S1: get the locations
    
    lon_m, lat_m=get_mean_location(lon, lat, flux, area)
    
    # S2: Distance
    
    dist=gcd.get_circle_distance(lon_m, lat_m)
    
    
    schur_factor=npy.zeros(npy.shape(dist), float)
    dist=dist/cor_length
    
    usd_idx=npy.where(dist<dist_threshold) # only close points are used
    
    # S3: schur factor
    
    schur_factor[usd_idx]=npy.exp(-dist[usd_idx])
    

    # S4: calculate the factor matrix=sqrt(schur_factor)
    
    
    u, w, v=nlg.svd(schur_factor)
    w=npy.sqrt(w)
    w=npy.diag(w)
    
    cor_factor=npy.dot(u, w)
    cor_factor=npy.dot(cor_factor, npy.transpose(u))
    
    
    
    
    return cor_factor


# <<< test >>>

if (__name__=='__main__'):
    
    
    varnames=['longitude', 'latitude', 'flux', 'area', 'map']
    flnm='reg_144_ml.4x5.nc'
    lon, lat, flux, area, regmap=nfio.ncf_read(flnm, varnames)
    
    cor_length=1000.0
    cor_factor_land=construct_region_spatial_cor_factor(lon, lat, \
                                                            flux[:,:, 0:100],\
                                                            area,\
                                                            cor_length, \
                                                            dist_threshold=6)

    plb.figure(3)
    
    plb.pcolor(cor_factor_land)
    
    plb.show()
    

    
    
    cor_len=2000.0

    cor_factor_ocean=construct_region_spatial_cor_factor(lon, lat, \
                                                             flux[:,:, 100:],\
                                                             area,\
                                                             cor_length, \
                                                             dist_threshold=6)
    
    plb.figure(4)
    
    plb.pcolor(cor_factor_ocean)
    
    plb.show()
    
