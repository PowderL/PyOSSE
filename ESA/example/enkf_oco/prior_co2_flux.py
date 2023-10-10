""" read co2 emission at model grid box and then convert them into regional fluxes

"""
import numpy as npy

import ESA.util.bpch2_rw_smp as brw # need to handle the data
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.geo_constant as gc
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.geo_constant as gc

XNUMOL_CO2=gc.An/(1.0e-3*gc.mco2)
flux_unit_factor=1.0/XNUMOL_CO2
flux_unit_c_co2=gc.mco2/gc.mc

# #c: notes: the unit is in kgCO2/s, which is different 
# #c: from outputs of functions in ESA.surface






def make_ctm_flnm(datapath, fctm, ctm_date0, nmext, \
                      ftracerinfo='tracerinfo', \
                      fdiaginfo='diaginfo'):
    
    """ 
    make ctm file name to read in flux
    CTM is saved at a time step different from GEOS-Chem  ts files
    
    Inputs:
    --------------------------------------------------------------
    1.datapath:<str>: directory of the model outputs
    2.fctm:<str>: name (format of name) of the ctm file 
    3.nmext: <str>: extension of the ctm file 
    4.ctm_date0:<str>: date tag for ctm file  
    5.ftracerinfo:<str>: name of the tracer info file
    6.fdiaginfo:<str>:name of diaginfo file 
    
    Returns:
    -------------------------------------------------------------
    1. ctmflnm, ftraceinfo, fdiaginfo:<str>: name of CTM, tracerinfo and diaginfo file 
    
    """
    fctm=fctm.replace('XEXTX', nmext)
    fctm=fctm.replace('XYYYYMMDDX', ctm_date0)
    
    ctmflnm=datapath+"/"+fctm
    
    print '1', ftracerinfo
    
    ftracerinfo=ftracerinfo.replace('XEXTX', nmext)
    ftracerinfo=ftracerinfo.replace('XYYYYMMDDX', ctm_date0)
    print '2', ftracerinfo
    
    fdiaginfo=fdiaginfo.replace('XEXTX', nmext)
    fdiaginfo==fdiaginfo.replace('XYYYYMMDDX', ctm_date0)

    ftracerinfo=datapath+"/"+ftracerinfo
    fdiaginfo=datapath+"/"+fdiaginfo
    return ctmflnm, ftracerinfo, fdiaginfo


def read_surface_emission_from_ctm(ctmflnm,\
                                       ftracerinfo,\
                                       fdiaginfo,\
                                       tau0,\
                                       src_cat= ['CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE'],\
                                       src_name=['FF', 'OC', 'BAL',  'BB',  'BF'],\
                                       src_id=[1,  2,  3, 4, 5]):
    
    
    
    """ read surface emission rate from CTM file 
    
    INPUTS: 
    ---------------------------------------------------------------
    1. ctmflnm:<str>: name of CTM file 
    2. ftraceinfo:<str>: name of tracerinfo file 
    3. fdiaginfo:<str>: name of diaginfo file 
    4. tau0: <float>: starting time for emissions
    5. src_cat:<list, t:str>: surface emission category
    6. src_name:<list, t:str>: surface emission name 
    7. src_id:<list, t:int>: surface emission id (check tracer info files)
    

    OUTPUTS: 
    -----------------------------------------------------------
    1. em_tau0, em_tau1:<float>: emission starting and end time 
    2. all_src:<dict>:surface emission in each grid box in unit of mole /s
    
    """
    
    
    # S1: form the selection criteria
    
    sel_categorys=src_cat
    
    nsrc=len(src_cat)
    
    if (src_id==None):
        src_id=[-1]*(nsrc) # set to be  -1
    
    sel_tracers=src_id
    
    sel_taus=[tau0]*(nsrc)
    
    
    # S2: load the emission data
    
    bpch2_ctm=brw.bpch2_file_rw(ctmflnm, "r", do_read=2,\
                                    sel_tracers=sel_tracers,\
                                    sel_taus=sel_taus,\
                                    sel_categorys=sel_categorys,\
                                    ftracerinfo=ftracerinfo,\
                                    fdiaginfo=fdiaginfo)
    
    # S3: summary over different emission type
    
    usd_cat=list()
    all_src={}
    nsrc=len(src_id)
            
    for isrc in range(nsrc):
    
        cat=src_cat[isrc]
        categorys=[cat]

        itra=src_id[isrc]
        
        if (itra>0):
            tracers=[itra]
        else:
            tracers=None
        
        sname=src_name[isrc]
        
        
        tranames=None
        taus=None
        
        data_list, founded=bpch2_ctm.get_data(categorys, tracers, taus, tranames)
        
        # print isrc, sname, len(data_list)
        
        bpdata=data_list[0]
        
        if (isrc==0):
            em_tau0, em_tau1=bpdata.get_attr(['tau0', 'tau1'])
        
        surface_src=npy.array(npy.squeeze(bpdata.data))
        print 'name, max, min', sname, npy.max(surface_src.flat), npy.min(surface_src.flat)
        
        all_src.update({sname:surface_src})
        
    # loop:  range(nsrc) end
    
    return em_tau0, em_tau1, all_src


def sum_reg_emission(region, sf_area, \
                     surface_src, src_key, nocean=44):

    """ calculate regional emissions from emissions at
    grid box and regional map for Transcom-3 23 regions. 
    INPUTS:
     region  -----in----- regional maps
     surface_src  ----in emission at grid box
    OUTPUTS:
      reg_src   ----- out---- the regional source
      
    """
    XNUMOL_CO2=gc.An/(1.0e-3*gc.mco2)
    DTSRCE=1.0e4
    factor=DTSRCE/XNUMOL_CO2
    
    nreg=npy.size(region[0,0,:])
    
    reg_src=npy.zeros(nreg, float)
    
    nland=nreg-nocean

    
    if (src_key=='OC'):
        nst, nend=nland, nreg

    elif (src_key=='FF'):
        # nst, nend=0, nreg
        nst, nend=0, nland
    
    else:
        nst, nend=0, nland
    
    
    
    tx=npy.sum(region[:,:, nst:nend], axis=2)
    usd_idx=npy.where(tx>0.0)
    
    mask=npy.ones(npy.shape(region[:,:,0]), float)
    mask[usd_idx]=0.0
    residual_src=surface_src*mask
    
    
    
    for ireg in range(nst, nend):
    
        reg_map=region[:,:,ireg]
        
        # reg_map=where(reg_map>0.1, 1, 0.0)
        # /tx
        
        flux_map1=reg_map[usd_idx]*surface_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]

        # flux_map1=flux_map1*sf_area

        flux_map1=npy.sum(flux_map1.flat)
        flux_map1=factor*flux_map1
        
        reg_src[ireg]=flux_map1
    

    
    
    if (nst>0):

        tx=npy.sum(region[:,:,0:nst], axis=2)
        usd_idx=npy.where(tx>0.0)
        
        for ireg in range(0, nst):
            reg_map=region[:,:,ireg]
            
            # reg_map=where(reg_map>0.1, 1, 0.0)
            # /tx
            
            flux_map1=reg_map[usd_idx]*residual_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]
            
            # flux_map1=flux_map1*sf_area
            
            flux_map1=npy.sum(flux_map1.flat)
            flux_map1=factor*flux_map1
            
            reg_src[ireg]=flux_map1
            
    if (nend<nreg):
        tx=npy.sum(region[:,:,nend:nreg], axis=2)
        usd_idx=npy.where(tx>0.0)
        
        for ireg in range(nend, nreg):
            reg_map=region[:,:,ireg]
            
            # reg_map=where(reg_map>0.1, 1, 0.0)
            # /tx
            
            flux_map1=reg_map[usd_idx]*residual_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]
            
            # flux_map1=flux_map1*sf_area
            
            flux_map1=npy.sum(flux_map1.flat)
            flux_map1=factor*flux_map1
            
            reg_src[ireg]=flux_map1
    
    
    reg_src=npy.array(reg_src)
    total_flux=factor*(sum(sf_area*surface_src))
    
    print 'sum_flux ', src_key, npy.sum(reg_src), total_flux
    total_flux=npy.sum(reg_src[55:])
    total_flux=npy.sum(reg_src)
    
    return reg_src, total_flux



def get_reg_emission_sum(reg_map, sf_area, src):
    """
    get sum for flux in regions defined by reg_map 
    Inputs:
    ----------------------------------------
    1. reg_map:<array, (nlon, nlat, nlayer)>: regional maps
    2. sf_area:<array, (nlon, nlat)>: grid area
    3. src: <array, (nlon, nlat)>: emission
    
    Returns:
    -------------------------------------------
    1. reg_src:<array, (nlayer)>: regional source

    """
    
    # S1: get dimensions
    nlon, nlat, nlayer=npy.shape(reg_map)
    reg_src=list()
    # S2: calculate flux over each region
    
    for ilayer in range(nlayer):
        reg_flux=sf_area[:,:]*reg_map[:,:,ilayer]*src[:,:]
        reg_src.append(npy.sum(reg_flux))
    
    reg_src=npy.array(reg_src)
    
    return reg_src





        



def get_gc_prior_flux( yyyy, mm, dd,\
                           datapath, fctm,\
                           nmext,\
                           ctm_date0, \
                           ftracerinfo,\
                           fdiaginfo,\
                           tau0,\
                           reg_map,\
                           sf_area,\
                           do_save=True,\
                           out_flnm='gc_co2',\
                           src_cat= ['CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE'],\
                           src_name=['FF', 'OC', 'BAL',  'BB',  'BF'],\
                           src_id=[1,  2,  3, 4, 5]):
    
    
    """ calculate the regional fluxes

    Inputs:
    -----------------------------------------
    1. yyyy, mm, dd:<int>: date
    2. datapath:<str>: directory of the model outputs
    2.fctm:<str>: name (format of name) of the ctm file 
    3.nmext: <str>: extension of the ctm file 
    4.ctm_date0:<str>: date tag for ctm file  
    5.ftracerinfo:<str>: name of the tracer info file
    6.fdiaginfo:<str>:name of diaginfo file 
    7.tau0:<float>: starting time for emissions
    8. reg_map:<array, (nlon, nlat, nem)>: regional map 
    9. sf_area:<array, (nlon, nlat)>: area
    10. do_save:<T/F>: save teh dtat 

    
    OUTPUTS:
    ----------------------------------------------
    
    1. total_reg_src:<array, (nlayer)>: total_regional emission for nreg in kgCO2/s
    
    """
    # S1: ensemble emissions
    
    
    # read in the emissions 
    
    print '-----get_gc_prior_flux-----' 
    full_flnm, full_ftracerinfo, full_fdiaginfo=make_ctm_flnm(datapath, fctm, ctm_date0, nmext, \
                                                               ftracerinfo=ftracerinfo, \
                                                               fdiaginfo=fdiaginfo)
    
    print 'Read flux from :', full_flnm
    print 'tracerinfo:' , full_ftracerinfo
    print 'diaginfo  :' , full_fdiaginfo
    
    
    
    em_tau0, em_tau1, all_src=read_surface_emission_from_ctm(full_flnm,\
                                                                 full_ftracerinfo,\
                                                                 full_fdiaginfo,\
                                                                 tau0,\
                                                                 src_cat= src_cat,\
                                                                 src_name=src_name,\
                                                                 src_id=src_id)
    


    print 'Time: (tau0, em_tau0, em_tau1)', em_tau0, em_tau1
    
    src_info_lst=list()
    isrc=0
    unit='kgCO2/s'
    
    
    for src_key in all_src.keys():
        
        surface_src=all_src[src_key]
        
        # reg_src, total_sf_src=sum_reg_emission(reg_map, sf_area, \
        #                                           surface_src, src_key, nocean=44)
        
        reg_src= get_reg_emission_sum(reg_map, sf_area, surface_src)
        
        reg_src=npy.squeeze(reg_src)
        
        reg_src=1.0e4*flux_unit_factor*reg_src
        
        # outputs
        
        src_info=ovar.io_var_cl(src_key, 'f', ['layer'], reg_src, varattr={"unit":unit, "desc":src_key})
        src_info_lst.append(src_info)
        if (isrc==0):
            total_reg_src=npy.array(reg_src)
        else:
            total_reg_src=total_reg_src+reg_src
        
        
        
        isrc=isrc+1
    
    
    nlayer=npy.size(total_reg_src)
    
    
    if (do_save):
        sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm,dd)
        xlayer=npy.arange(nlayer)
        layer_info=ovar.io_var_cl('layer', 'i', ['layer'], xlayer, \
                                      varattr=None)
        out_ncflnm=out_flnm+"."+sdate+".nc"
        ncfio.ncf_save_var(out_ncflnm,  src_info_lst, [layer_info], \
            create_new=True)   
        
    
    return  em_tau0, em_tau1, total_reg_src



def create_bf_read_cl(yyyy, mm, dd, cfg):
    
    """
    Inputs:
    -----------------------------------------
    1. yyyy, mm, dd:<int>: date
    2. cfg:<menu_cl>: configuration for bf
    """
    
    cur_bf=bfcl.bf_cl(cfg['bf.path'], \
                          cfg['bf.flnm'], \
                          yyyy, mm, dd,  \
                          varname_lst=cfg['bf.varname_lst'],\
                          varname_dict=cfg['bf.varname_dict'],\
                          fopen=cfg['bf.fopen'], \
                          fread=cfg['bf.fread'], \
                          fget=cfg['bf.fget'], \
                          fclose=cfg['bf.fclose'], \
                          fio_keywords=fio_keys)
    return cur_bf


def read_reg_bf(cur_bf, yyyy, mm, dd, \
                    varname_lst=['lon', 'lat', 'map', 'flux', 'area']):
    
    """
    read in regional flux

    Inputs:
    -----------------------------------------
    1. yyyy, mm, dd:<int>: date
    2. varname_lst:name of variables 

    Returns:
    --------------------------------
    1. 
    

    """
    # #c: longitude
    
    cur_bf.load_bf("", "", yyyy, mm, dd)
    lon, lat, reg_map, area, reg_flux=cur_bf.get_bf(0, varname_lst=varname_lst)
    reg_err=list()
    nlon, nlat, nlayer=npy.shape(reg_map)
    dims=npy.shape(area)
    print dims
    print npy.shape(reg_flux)
    
    for ilayer in range(nlayer):
        pb_flux=reg_flux[:,:,ilayer]*area[:,:]
        pb_flux=npy.sum(pb_flux)
        pb_flux=1.0e4*flux_unit_factor*pb_flux
        
        reg_err.append(pb_flux)
        
    return lon, lat, reg_map,  area, reg_err



if (__name__=="__main__"):
    print 'check oco_assim_step'
    
    
    
    
       
    
