"""
    Functions for doing composition of regional Basis functions (i.e, error)
    

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.24
    
    
    Functions:
    ===========================================
    1. do_region_err_svd: decomposition of regional error 
    
"""

import numpy as npy
import numpy.random as rnd
import numpy.linalg as nlg
import ESA.util.flux_regrid as fgrd
import spatial_cor as spcor
import pylab as plb

def do_region_err_svd(grd_lon, grd_lat, \
                          reg_map, reg_err, \
                          cor_length, \
                          threshold=0.1, \
                          scale=1.0, \
                          dist_threshold=4.0, \
                          maxsel=20,\
                          **keywords):
    
    """
    
    SVD error map over region to several  sub-regions defined by (nsb_lon, nsb_lat)
    
    Inputs:
   
    --------------------------------------
    1. grd_lon:<array, (nlon)>: longitude
    2. grd_lat:<array, (nlat)>: latitude
    3. reg_map:<array, (nlon, nlat)>: ratio  map defining one region. 
    4. reg_err:<array, (nlon, nlat)>: err (or BF) defined over given region
    5. cor_length <float>: length of the correlation 
    6. threshold:<float>: all grid boxes with flux (err) below this threshold will be treated as one single region 
    7. scale:<float>: scaling factor applied to this region 
    8. dist_threshold:<float>: Max distance/cor_length where correlation will be considered 
    9. maxsel:<int>: the maximum number of eigen vector to be kept 
    
    Returns:
    --------------------------------------------
    1. flux:<array, (nlon, nlat, neig)>: new regions masked with sub-region number 
    2. new_map_lst:<list, t:array>: the layer map with each layer for each sub-regions
        
    
    """
        
    print '--->svd_region_flux:<---'
    
    # S1:  array of points inside the region 
    
    
    nx, ny=npy.shape(reg_map)
    usd_id=npy.where(reg_map>0)
    
    # #T: lon, lat
    idx=usd_id[0]
    idy=usd_id[1]
    npt=npy.size(idx)
    
    

    err=reg_err[idx, idy]  # from map to err array 

    mean_err=npy.mean(abs(err))  # the mean value
    total_err=npy.sum(err*err)
    
    print '1. mean_err:', mean_err
    print '2. total_err:', total_err
    
    # #T2:separate region into a small flux (err) one and another with larger fluxes
    
    sflux_id=npy.where(err<threshold*mean_err)
    
    if (npy.size(sflux_id)>0):
        sflux_idx=idx[sflux_id]
        sflux_idy=idy[sflux_id]
        sflux_err=err[sflux_id]
        sflux_map=npy.where(reg_err<threshold*mean_err, reg_map, 0)
    
    bflux_id=npy.where(err>=threshold*mean_err)
    bflux_idx=idx[bflux_id]
    bflux_idy=idy[bflux_id]
    bflux_err=err[bflux_id]
    
    bflux_map=npy.where(reg_err>=threshold*mean_err, reg_map, 0)
    
        
            
        
    # S2:  calculate correlation 
    
    clon=grd_lon[bflux_idx] # lons of  region cells
    clat=grd_lat[bflux_idy] # lats of  region cells
    
    err_cor=spcor.construct_grid_spatial_err_cov(clon, clat, \
                                                     bflux_err, \
                                                     cor_length, \
                                                     scale, \
                                                     dist_threshold=dist_threshold)
    
    
    
    sum_org=npy.sum(bflux_err*bflux_err)
    print '3. total error^2 for svd:', sum_org
    
    # S3: SVD 
    
    u, w, v=nlg.svd(err_cor)
    
    
    # plb.plot(w)
    # plb.show()
    
    
    
    
    
    if (maxsel==None):
        maxsel=npy.size(w)
    else:
      if (maxsel>npy.size(w)):
          maxsel=npy.size(w)
    
      
    sum_w30=npy.sum(w[0:maxsel])
    sum_w=npy.sum(w)
    dz=sum_w30/sum_w
    
    print '4. ratio of retained error after SVD', dz
    
    w=npy.sqrt(w)
    
    w=npy.diag(w)
    
    
    # S4:  select the vector to be kept 
    
    sel_err_vect=npy.dot(u, w)
    sel_err_vect=sel_err_vect[:,0:maxsel]
    ## c: applying transpose(u) will localize vector_err
    
    #  vector_err=npy.dot(vector_err, npy.transpose(u))
    
    new_err_cor=npy.dot(sel_err_vect, npy.transpose(sel_err_vect))
    new_total_err=npy.sum(new_err_cor)
    


    print '5. new total err^2 after select:', new_total_err
    
    # #T1: rescale so that total of (err_cor)==total of (bflux_err)
   
    
    scaling_factor=npy.sqrt(total_err/new_total_err)
    print '6. scaling factor back to original total:', scaling_factor
    
    sel_err_vect=scaling_factor*sel_err_vect
    
    

    
    if (npy.size(sflux_id)>0):
        
        sflux_at_grid=npy.zeros(npy.shape(sflux_map), float)
        sflux_at_grid=fgrd.add_val_to_map(sflux_at_grid, sflux_err, sflux_idx, sflux_idy)
        sum_sflux=npy.sum(sflux_at_grid)
        # #c: as now sflux and sel_err_vect independent 
        
        new_total_err=total_err+sum_sflux*sum_sflux
        
        # rescaling 
        scaling_factor=npy.sqrt(total_err/new_total_err)
        sflux_at_grid=scaling_factor*sflux_at_grid
        sel_err_vect=scaling_factor*sel_err_vect
        print '7. new total_err inc.small:', new_total_err
        print '8. scaling factor inc small:', scaling_factor
    
        
    # S5:  fill svd_map_lst, svd_flux_lst
    svd_map_lst=list()
    svd_flux_lst=list()
    
    
    
    for isel in range(maxsel):
        # #c: still use the regional map 
        svd_map_lst.append(bflux_map)
        sel_u=sel_err_vect[:,isel]
        svd_flux=npy.zeros(npy.shape(bflux_map), float)
        
        svd_flux2=fgrd.add_val_to_map(svd_flux, sel_u, bflux_idx, bflux_idy)
        
        svd_flux_lst.append(svd_flux2)
        
    
    if (npy.size(sflux_id)>0):
        svd_flux_lst.append(sflux_at_grid)
        svd_map_lst.append(sflux_map)
        
    
        
    
    out_dict={'lon':grd_lon, \
                  'lat':grd_lat, \
                  'map':svd_map_lst,\
                  'flux':svd_flux_lst}
    
    print '8. size of map; flux:', len(svd_map_lst), len(svd_flux_lst)
    print '---<svd_region_flux >----'
    
    return out_dict

    


        
