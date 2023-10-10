"""
    Functions for sampling flux perturbations to generable flux ensemble 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.24
    
    
    Functions:
    ===========================================
    1. sample_by_reg: Sample every regional bf
    2. random_sample_reg: Randomly sample regional bf 
    3. svd_sample_reg: Sample regional bf by using SVD of error covariance.
    

"""

import numpy as npy
import numpy.random as rnd
import numpy.linalg as nlg
import ESA.util.flux_regrid as fgrd
import spatial_cor as spcor
import pylab as plb
import ESA.util.otool_ncfile_io as ncio



def sample_by_reg(grd_lon,\
                      grd_lat,\
                      reg_flux, \
                      reg_map, \
                      pid, \
                      **keywords):
    """
    
    Sample every regional bf
    
    Inputs:
    ---------------------------------------------------------------------
    1. grd_lon:<array, (nlon)>: longitude
    2. grd_lat:<array, (nlat)>: latitude
    
    3. reg_flux:<array, (nlon, nlat, nlayer)>: regional flux. 
       It could be single layer 

    4. reg_map:<array, (nlon, nlat, nlayer)>: regional map 
    5. pid:   <array, (nlayer/nreg)>: the parent (T3) region ID
    
    6. keywords:<dict>: the extra inputs for sampling 
    ---reserved keywords:
    --->sel_reg_lst:<list, t:int>: List of region ID to be included 
    
    Outputs:
    --------------------------------------------------------------------
    1. out_dict:<dict>: It has 6 entries
    
    ---lon:<array, (nlon)>: longitude
    ---lat:<array, (nlon)>: latitude
    ---flux:<array, (nlon, nlat, nsel)>: selected regional flux perturbations 
    ---map:<array, (nlon, nlat, nsel)>:  regional map 
    ---pid:<array, (nsel)>: Parent ID of the selected regions
    ---sel_reg_lst:<array, t:int>: selected layers 
    ---coef:<array, (nsel, nsel)>: projection from seleced region to the out put
    

    """
    # S1: get list for regions to be included
    
    is_ml_map=False
    
    dims=npy.shape(reg_map)
    if (len(dims)==3):
        nlayer=dims[2]
        is_ml_map=True
    else:
        nlayer=int(npy.max(reg_map.flat))
        if (nlayer==0):
            nlayer=1


    sel_reg_lst=None
    
    
    if ('sel_reg_lst' in keywords):
        sel_reg_lst=keywords['sel_reg_lst']
        
    
    if (sel_reg_lst==None):
        sel_reg_lst=range(1, nlayer+1)
    

    # S2: convert reg_id from (1 to nlayer) to (0 to nlayer-1)
    
    id_reg_lst=npy.array(sel_reg_lst)

    id_reg_lst=id_reg_lst-1

    id_reg_lst=id_reg_lst.tolist()
    
    # print id_reg_lst
    # print type(id_reg_lst)
    
    
    if (is_ml_map):
        # ##c: multi-layer flux map
        
        
        sel_flux=reg_flux[:,:,id_reg_lst]
        sel_map=reg_map[:,:,id_reg_lst]
        sel_pid=pid[id_reg_lst]
        nreg=len(id_reg_lst)
        coef=npy.identity(nreg)
        
        out_dict={'lon':grd_lon,\
                  'lat':grd_lat,\
                      'flux':sel_flux, 'map':sel_map, \
                      'pid':sel_pid, \
                      'sel_reg_lst':sel_reg_lst, \
                      'coef':coef}

        
        return out_dict
    
    else:
        
        # ##c: single-layer flux map
        
        sel_map=[]
        sel_flux=[]
        sel_pid=[]
        
        for ireg in sel_reg_lst:
            
            layer_flux=npy.where(sel_map==ireg, reg_flux, 0)
            sel_flux.append(layer_flux)
            layer_map=npy.where(sel_map==ireg, 1, 0)
            sel_map.append(layer_map)
            sel_pid.append(pid[ireg-1])
        
        nreg=len(id_reg_lst)
        coef=npy.identity(nreg)
        
        sel_flux=npy.array(sel_flux)
        sel_map=npy.array(sel_map)
        sel_pid=npy.array(sel_pid)
        
        sel_flux=npy.transpose(sel_flux, [1,2,0])
        sel_map=npy.transpose(sel_map, [1,2,0])
        sel_reg_lst=npy.array(sel_reg_lst)
        
        out_dict={'lon':grd_lon,\
                  'lat':grd_lat,\
                      'flux':sel_flux, 'map':sel_map, \
                      'pid':sel_pid,\
                      'sel_reg_lst':sel_reg_lst, \
                      'coef':coef}
        
        
        return out_dict
 
        



def random_sample_reg(grd_lon,\
                          grd_lat,\
                          reg_flux, \
                          reg_map, \
                          pid, \
                          **keywords):
    
    
    """
    
    Randomly sample regional bf 
    
    Inputs:
    ---------------------------------------------------------------------
    1. grd_lon:<array, (nlon)>: longitude
    2. grd_lat:<array, (nlat)>: latitude
    
    3. reg_flux:<array, (nlon, nlat, nlayer)>: regional flux. 
    It could be single layer 
    
    4. reg_map:<array, (nlon, nlat, nlayer)>: regional map 
    5. pid:   <array, (nlayer/nreg)>: the parent (T3) region ID
    
    6. keywords:<dict>: the extra inputs for sampling 
    
    ---reserved keywords:
    --->sel_reg_lst:<list, t:int>: List of region ID to be included 
    --->nsample:<int>: the upper-limit of the set of the sampled fluxes
    --->scale:<array, (nsample)>: values to scaling BF functions
        
    Outputs:
    --------------------------------------------------------------------
    1. out_dict:<dict>: It has seven entries
    
    ---grd_lon:<array, (nlon)>: longitude
    ---grd_lat:<array, (nlat)>: latitude
    ---flux:<array, (nlon, nlat, nsample)>: combined regional flux perturbations 
    ---map:<array, (nlon, nlat, nsel_reg>: individual regional map 
    ---sel_reg_lst:<list/array, t:int(nsel_reg)>: selected layers 
    ---coef:<array, (nsel_reg, nsample>: random coefficient 
    ---pid:<array, nsel_reg>: parent ID for selected region 
    
        
    """
    
    # #S1: check parameter for random sampling 
    
    # #c: check whether it is a multi-later map
    
    is_ml_map=False
    
    dims=npy.shape(reg_map)
    if (len(dims)==3):
        nlayer=dims[2]
        is_ml_map=True
    else:
        nlayer=int(npy.max(reg_map.flat))
        if (nlayer==0):
            nlayer=1

    sel_reg_lst=None
    
    # #c: check sel_reg_lst
    
    if ('sel_reg_lst' in keywords):
        sel_reg_lst=keywords['sel_reg_lst']
        
    if (sel_reg_lst==None):
        sel_reg_lst=range(1, nlayer+1)
    
    # #c: check number of sample to be kept. 
    
    if ('nsample' in keywords):
        nsample=keywords['nsample']
    else:
        nsample=nlayer


    
    
    # #c: check the scale 
    
    scale=npy.ones(nsample, float)
    
    if ('scale' in keywords):
        scale=keywords['scale']
    
    # S2: generate the random coefficients
    
    # #c: sel_coef in shape  [nreg, nsample]
    
    sel_coef=[]
    nsel_reg=len(sel_reg_lst)
    
    for i in range(nsel_reg):
        
        pb_scale=scale[i]
        coef_rnd=rnd.normal(scale=pb_scale, size=nsample)
        # #c: renormalize
        
        coef_rnd=coef_rnd/npy.sqrt(nsample-1)
        
        sel_coef.append(coef_rnd)
    
    sel_coef=npy.array(sel_coef)
    
    
    if (is_ml_map):
        # ##c: multi-layer map
        
        # #T1: convert reg_id from (1 to nlayer) to (0 to nlayer-1)
        
        id_reg_lst=npy.array(sel_reg_lst)
        
        # ##c: to python index
        
        id_reg_lst=id_reg_lst-1
        
        sel_pid=pid[id_reg_lst]
        
        usd_flux=reg_flux[:,:,id_reg_lst]
        
        # #T2: apply random coefficient to generate combined fluxes. 
        
        sel_flux=fgrd.do_flux_by_coef(usd_flux, sel_coef)
        sel_map=reg_map[:,:,id_reg_lst]
        
        out_dict={'lon':grd_lon, 'lat':grd_lat, \
                      'flux':sel_flux, 'pid':sel_pid, \
                      'coef':sel_coef, \
                      'map':sel_map,'sel_reg_lst':sel_reg_lst}
        return out_dict
    
    
    else:
        # ##c: single layer map 
        # #T1: construct multiple-layer  regional flux and map 
        
        sel_flux=[]
        usd_map=[]
        sel_pid=[]
        for ireg in sel_reg_lst:
            cur_flux=npy.where(sel_map==ireg, reg_flux, 0)
            cur_map=npy.where(sel_map==ireg, 1, 0)
            usd_flux.append(cur_flux)
            sel_map.append(cur_map)
            sel_pid.append(pid[ireg-1])
        
        usd_flux=npy.array(usd_flux)
        usd_flux=npy.transpose(usd_flux, [1,2,0])
    
        sel_map=npy.array(sel_map)
        sel_map=npy.transpose(sel_map, [1,2,0])
    
        # #T2 apply random coefficient 
        
        sel_flux=grd.do_flux_by_coef(usd_flux, sel_coef)
        
        out_dict={'lon':grd_lon, 'lat':grd_lat, \
                      'pid':sel_pid,\
                      'flux':sel_flux, 'coef':sel_coef, \
                      'map':sel_map,'sel_reg_lst':sel_reg_lst}
        
        return our_dict
    
    

def svd_sample_reg(grd_lon,\
                       grd_lat,\
                       reg_flux, \
                       reg_map, \
                       pid, \
                       **keywords):
    
    
    """
    
    SVD sample regional bf 
    
    Inputs:
    ---------------------------------------------------------------------
    1. grd_lon:<array, (nlon)>: longitude
    2. grd_lat:<array, (nlat)>: latitude
  
    3. reg_flux:<array, (nlon, nlat, nlayer)>: regional flux. 
    It could be single layer 
    
    4. reg_map:<array, (nlon, nlat, nlayer)>: regional map 
    3. pid:   <array, (nlayer/nreg)>: the parent (T3) region ID
    
    5. keywords:<dict>: the extra inputs for sampling 
    
    ---reserved keywords:
    --->sample_numb:<int>: the upper-limit of the size of sampled flux
    --->sel_reg_lst:<list, t:int>: List of region ID to be included 
    --->err_cor:<array>: error correlation matrix to be SVD
    
    
        
    Outputs:
    --------------------------------------------------------------------
    1. out_dict:<dict>: It has three entries
    ---flux:<array, (nlon, nlat, nsel_reg)>: combined regional flux perturbations 
    ---map:<array, (nlon, nlat, nsel_reg>: 
    ---sel_reg_lst:<array, t:int>: selected layers 
    ---sel_coef:<array, t:int>: selected 
    
        
    """
    
    # S1: check parameter for sampling 
    
    is_ml_map=False
    
    dims=npy.shape(reg_map)
    if (len(dims)==3):
        nlayer=dims[2]
        is_ml_map=True
    else:
        nlayer=int(npy.max(reg_map.flat))
        if (nlayer==0):
            nlayer=1

    sel_reg_lst=None
    
    # #c: check sel_reg_lst
    
    if ('sel_reg_lst' in keywords):
        sel_reg_lst=keywords['sel_reg_lst']
        
    if (sel_reg_lst==None):
        sel_reg_lst=range(1, nlayer+1)
    
    # #c: check number of sample to be kept. 
    
    if ('nsample' in keywords):
        nsample=keywords['nsample']
    else:
        nsample=nlayer


    # S4: SVD err_cor
    
    # ##c: generate the coefficients by SVD
    
    err_cov=keywords('err_cov')
    u,w, v=nlg.svd(err_cov)
    wd=npy.sqrt(w)
    wd=npy.diag(wd)
    u=npy.dot(u, wd)
    sel_coef=u[:, 0:sample_number]
    
    
    if (is_ml_map):
        # ##c multi-layer map
        
        # ##T1:: convert reg_id from (1 to nlayer) to (0 to nlayer-1)
        
        id_reg_lst=npy.array(sel_reg_lst)
        
        id_reg_lst=id_reg_lst-1
    
        usd_flux=reg_flux[:,:,id_reg_lst]
        sel_flux=grd.do_flux_by_coef(usd_flux, sel_coef)
        sel_map=reg_map[:,:,sel_reg_lst]
        out_dict={'flux':sel_flux, 'sel_coef':sel_coef, \
                      'map':sel_map,'sel_reg_lst':sel_reg_lst, 'coef':sel_coef}
        
        return out_dict


    
    
    else:
        # ##c: single layer map 
        # #T1: construct multiple-layer  regional flux and map 
        
        sel_flux=[]
        usd_map=[]
        for ireg in sel_reg_lst:
            cur_flux=npy.where(sel_map==ireg, reg_flux, 0)
            cur_map=npy.where(sel_map==ireg, 1, 0)
            usd_flux.append(cur_flux)
            sel_map.append(cur_map)
        
        usd_flux=npy.array(usd_flux)
        usd_flux=npy.transpose(usd_flux, [1,2,0])
    
        sel_map=npy.array(sel_map)
        sel_map=npy.transpose(sel_map, [1,2,0])
    
        # #T2 apply random coefficient 
        
        sel_flux=grd.do_flux_by_coef(usd_flux, sel_coef)
        
        out_dict={'lon':grd_lon, 'lat':grd_lat, 'flux':sel_flux, 'sel_coef':sel_coef, \
                      'map':sel_map,'sel_reg_lst':sel_reg_lst, 'coef':sel_coef}
        
        
        return out_dict



    
# <<< TEST >>>

if (__name__=='__main__'):
                  
    import bf_file_m as bfio
    varname_lst=bfio.bf_varname_lst
    
    varname_lst[0]='longitude'
    varname_lst[1]='latitude'
    
    flnm='reg_144_ml.4x5.nc'
    lon, lat,layer, reg_map, reg_flux,area, pid=ncio.ncf_read(flnm, varname_lst)
    var_dict=random_sample_reg(lon,\
                                   lat,\
                                   reg_flux, \
                                   reg_map, \
                                   pid)
    
    
