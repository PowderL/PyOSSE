""" read co2 emission at model grid box and then convert them into regional fluxes

"""
from pylab import *
import ESA.util.bpch2_rw_smp as brw # need to handle the data
import ESA.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.geo_constant as gc


def make_ctm_flnm(datapath, fctm, ctm_date0, nmext, \
                      ftracerinfo='tracerinfo', \
                      fdiaginfo='diaginfo'):
    
    """ ------------------------------------------------------------------
    ctm_date0 ----in ------ the starting date for ctm file  
    datapath -----in ------- the directory of the model outputs
    em_st, em_end  -------in------- the tagged tracer numbers.
               
    the name of CTM file is supposed to ctm.EN[em_st]--EN[em_end].[ctm_date0].bpch
    
    """
    
    ctmflnm=datapath+"/"+fctm+"."+nmext+"."+ctm_date0+".bpch"
    ftraceinfo=datapath+"/"+ftracerinfo+"."+nmext+".dat"
    fdiaginfo=datapath+"/"+fdiaginfo+"."+nmext+".dat"
    return ctmflnm, ftraceinfo, fdiaginfo


def read_surface_emission_from_ctm(ctmflnm,\
                                       ftraceinfo,\
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
                                    ftracerinfo=ftraceinfo,\
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
        
        surface_src=npy.array(squeeze(bpdata.data))
        all_src.update({sname:surface_src})
        
    # loop:  range(nsrc) end
    
    return em_tau0, em_tau1, all_src






def sum_reg_emission(region, reg_pid, sf_area, \
                     surface_src, src_key, nocean, factor):

    """ calculate regional total emissions 
    
    INPUTS:
    -------------------------------------------------- 
    1. region:<npy.array, (nlon, nlat, nlayer)>: regional maps
    2. reg_pid:<npy.array, (nlayer)>: parent id of the regions 
    3. sf_area:<npy.array, (nlon, nlat)>: area (in m2 by default)
    4. src_key:<list, t:str>: name of emission sources
    5. nocean:<int>: name of ocean in the regional maps
    
    
    OUTPUTS:
    ------------------------
    reg_src:<list>: the regional source
      
    """
    
    nreg=npy.size(reg_pid)
    reg_src=zeros(nreg, float)
    nland=nreg-nocean
    
    surface_src=factor*surface_src
    
    # S1: choose which land or ocean regions  will be calculated 

    
    if (src_key=='OC'):
        
        nst, nend=nland, nreg
        
    elif (src_key=='FF'):
        # nst, nend=0, nreg
        
        nst, nend=0, nland
    
    else:
        
        nst, nend=0, nland
        
        
    
    tx=sum(region[:,:, nst:nend], axis=2)
    usd_idx=where(tx>0.0)
    
    mask=ones(shape(region[:,:,0]), float)
    mask[usd_idx]=0.0
   
    # #c: residual_src is for emissions have been allocated to different surface type
    
    residual_src=surface_src*mask
    
    
    
    for ireg in range(nst, nend):
    
        reg_map=region[:,:,ireg]
        
        flux_map1=reg_map[usd_idx]*surface_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]

        # flux_map1=flux_map1*sf_area

        flux_map1=sum(flux_map1.flat)
        flux_map1=flux_map1
        
        reg_src[ireg]=flux_map1
    

    
    # S3: use residual to fill other regions corresponding to different types
    
    if (nst>0):

        tx=sum(region[:,:,0:nst], axis=2)
        usd_idx=where(tx>0.0)
        
        for ireg in range(0, nst):
            reg_map=region[:,:,ireg]
            
            # reg_map=where(reg_map>0.1, 1, 0.0)
            # /tx
            
            flux_map1=reg_map[usd_idx]*residual_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]
            
            # flux_map1=flux_map1*sf_area
            
            flux_map1=sum(flux_map1.flat)
            flux_map1=flux_map1
            
            reg_src[ireg]=flux_map1
            
    if (nend<nreg):
        tx=sum(region[:,:,nend:nreg], axis=2)
        usd_idx=where(tx>0.0)
        
        for ireg in range(nend, nreg):
            reg_map=region[:,:,ireg]
            
            # reg_map=where(reg_map>0.1, 1, 0.0)
            # /tx
            
            flux_map1=reg_map[usd_idx]*residual_src[usd_idx]*sf_area[usd_idx]/tx[usd_idx]
            
            # flux_map1=flux_map1*sf_area
            
            flux_map1=sum(flux_map1.flat)
            flux_map1=flux_map1
            
            reg_src[ireg]=flux_map1
    
    
    reg_src=npy.array(reg_src)
    total_flux=(sum(sf_area*surface_src))
    
    total_flux2=sum(reg_src)
    
    return reg_src, total_flux





def read_region_map(sdate, \
                        varnames=['map', 'area', 'parent_id'],\
                        datapath='./surface_flux/',\
                        prefix="CO2_EMISSION"):
    
    """ read the predefined regional masks  """
    
    

    map_flnm=datapath+'/'+prefix+'.'+sdoy+'.nc'
    region, sf_area, reg_pid=ofb.ncf_read(map_flnm, varnames)
    nlayer=npy.size(reg_pid)
    sf_area=1.0e4*sf_area
    
    tx=sum(region[:,:,:], axis=2)
    print max(tx.flat), min(tx.flat)
    region=region[:,:,:]/tx[:,:,newaxis]
    print 'end-divid'
    
    return region, sf_area, reg_pid


        

def read_reg_flux_pb(sdoy,\
                         region, sf_area,\
                         datapath='./surface_flux/',\
                         prefix="CO2_EMISSION"):
    
    """ read region flux and constructur error correlation
    """
    
    XNUMOL_CO2=gc.An/(1.0e-3*gc.mco2)
    DTSRCE=1.0
    factor=DTSRCE/XNUMOL_CO2
    
    
    varnames=['flux']
    map_flnm=datapath+'/'+prefix+'.'+sdoy+'.nc'
    flux=ofb.ncf_read(map_flnm, varnames)
    flux=squeeze(flux)
    print shape(flux)
    
    nreg=npy.size(flux[0,0,:])
    err=list()
    # error listed
    
    
    # idx for cells in one region
    
    reg_idx_lst=list()
    # cell flux in one region
    
    
    reg_err_lst=list()

    
    for ireg in range(nreg):
        reg_pb_flux=flux[:,:,ireg]
        reg_map=region[:,:,ireg]
        used_idx=where(reg_map>0.1)
        # used_idx=squeeze(used_idx)
        # 
        reg_err=factor*reg_pb_flux[used_idx]*sf_area[used_idx]
        
        reg_err_lst.append(reg_err)
        reg_idx_lst.append(used_idx)
        
        
        uflux=reg_pb_flux # [used_idx]
        uarea=sf_area #[used_idx]
        reg_pb=sum(uflux*uarea)
        reg_pb=factor*reg_pb
        
        err.append(reg_pb)
        
    err=npy.array(err)
    
    return err,reg_idx_lst, reg_err_lst


def get_reg_emissions_info(yyyy, mm, dd,\
                           em_st=1, em_end=2, \
                           err_flux_path='./surface_flux/',\
                           prior_flux_path=gcdf.data_path,\
                           do_write=True, \
                           nland=100,\
                           prefix="co2_emission",\
                           unit='kgCO2/s'\
                           ):
    
    """ calculate the regional fluxes

    OUTPUTS:
    
    total_reg_src ---out--- regional emission for nreg in kgCO2/s
    reg_err ----- out------- region emission for error in kgCO2/s
    total_glb_src ---- total globle emissions kgCO2/s
      
      
    
    """
    
    doy=tmdl.day_of_year(yyyy, mm, dd)
    sdate=r'%4.4dD%3.3d' % (yyyy, doy)
    
    region_m, sf_area, reg_pid=read_region_map(sdate, \
                                                   datapath=err_flux_path)
    
    
    utc=tmdl.time_array_to_utc(yyyy,mm, dd)
    ctm_date0=r'%4.4d%2.2d%2.2d' % (yyyy, mm,dd)
    tau0=tmdl.utc_to_tai85(utc)
    tau0=tau0/(3600.0)
    
    src_cat= ['CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE', 'CO2-SRCE']
    src_name=['FF', 'OC', 'BAL',  'BB',  'BF']
    src_id=[1,  2,  3, 4, 5]
    # read in the emissions 
    
    em_tau0, em_tau1, all_src=read_surface_emission(yyyy,\
                                                        tau0, \
                                                        ctm_date0,\
                                                        prior_flux_path,\
                                                        em_st=em_st, em_end=em_end, \
                                                        src_cat= src_cat,\
                                                        src_name=src_name,\
                                                        src_id=src_id \
                                                        )
    
    src_info_list=list()
    reg_src_list=list()
    reg_err_list=list()
    
    isrc=0
    total_glb_src=list()
    
    for src_key in all_src.keys():
        
        surface_src=all_src[src_key]
        reg_src, total_sf_src=sum_reg_emission(region_m, reg_pid, \
                                                   sf_area, surface_src,src_key)
        # reg_src=factor*reg_src
        reg_src=squeeze(reg_src)
        
        all_name="Sub T3 Regions" 
        
        total_glb_src.append(total_sf_src)
        
        src_info=ofb.geos_varinfo(src_key, 'f', ['nreg'], reg_src, varattr={"unit":unit, "desc":all_name})
        src_info_list.append(src_info)
            
        if (isrc==0):
            total_reg_src=npy.array(reg_src)
        else:
            total_reg_src=total_reg_src+npy.array(reg_src)

        if (src_key=='BAL'):
            src_bio=npy.array(reg_src)
        
        if (src_key=='OC'):
            src_oc=npy.array(reg_src)
            
        
        isrc=isrc+1
    
    # get the regional emission errors
    
    
    reg_err,reg_idx_lst,reg_err_lst  =read_reg_flux_errors(sdate, \
                                                               region_m, sf_area,\
                                                               datapath=err_flux_path)
    
    # reg_err=diag(reg_err)
    
    print npy.size(reg_err), npy.size(src_bio), npy.size(src_oc), max(region_m.flat)

    nreg=npy.size(reg_err)

    src_info=ofb.geos_varinfo('err', 'f', ['nreg'], \
                                  reg_err, varattr={"unit":unit, "desc":"total uncertainties in TransCom 3"})
    
    src_info_list.append(src_info)
    
    if (do_write):
        sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm,dd)
        xreg=arange(nreg)
        dimnames=['nreg']
        dimvars=[xreg]
        dimtypes=['i']
        all_name=all_name.strip()
        resflnm=prefix+"."+sdate+".nc"
        ofb.ncf_write_by_varinfo(resflnm, dimnames, dimtypes, dimvars, src_info_list)
        
    
    return total_reg_src, reg_err, total_glb_src


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
    
    cur_bf.load_bf(yyyy, mm, dd)
    lon, lat, ref_map, reg_flux, area=cur_bf.get_bf(0, varname_lst=varname_lst)
    reg_err=list()
    
    for ilayer in range(nlayer):
        pb_flux=reg_flux[:,:,ilayer]*area
        pb_flux=npy.sum(pb_flux)
        reg_err.append(pb_flux)
        
    return lon, lat, reg_map,  area, reg_err



    



if (__name__=="__main__"):
    from pylab import *
    import oco_units as ount
    yyyy, mm, dd=2009, 1, 1
    for mm in range(1,13):
        x1, x2, x3=get_reg_emissions_info(yyyy, mm, dd,\
                                          em_st=1, em_end=2, \
                                          err_flux_path='./surface_flux/', \
                                          prior_flux_path='./enkf_rerun/',\
                                          do_write=True, \
                                          prefix="co2_emission",\
                                          unit='kgCO2/s'\
                                          )
    
        print shape(x1), shape(x2)
        x3=array(x3)
        if (mm==1):
            total_x3=array(x3)
        else:
            total_x3=total_x3+x3

        print ount.kg_s_to_GtC_Y*x3
    
    
    print ount.kg_s_to_GtC_Y*total_x3/12
    
    
    
    
       
    
