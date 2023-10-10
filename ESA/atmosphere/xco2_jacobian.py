""" set up jacobian for monthly co inversion
"""
from pylab import *
import bpch2_rw_v2 as brw # need to handle the data
import gp_axis as gax
import gp_field as gf
import geo_constant as gc
import geos_chem_def as gcdf
import time_module as tmdl
import pres_mod_py as pm
import gen_plots as gpl
import read_cloud_map as rcm
# import scipy.weave as weave
import field_read as frd
import oco_feedback as ofb
import geo_constant as gc
import grid_prof as grd_pf
import flib as flb
from numpy import *
import standard_atmosphere as stam


def setup_vert_intpl(grd_pres, mod_pres, obs_pres):

    """
    calculate
    
    1) vpl, vpr, vwgt

    location (vpl/vpr) and  coefficicents (vwgt) for vertical interpolations from model grid (mod_pres)
    to observation vertical grid (obs_pres)

    2) col_wgt  
    coefficent for calculating column integration (i.e., mass weight)

    """
    
    do_vert_reverse=False
    
    if (mod_pres[0,2]>mod_pres[0,6]):
        do_vert_reverse=True
        lg_mod_pres=log10(mod_pres[:, ::-1])
    else:
        lg_mod_pres=log10(mod_pres)
        
    usd_idx=where(grd_pres>-990.0)
    lg_grd_pres=array(grd_pres)
    lg_grd_pres[usd_idx]=log10(grd_pres[usd_idx])
    # lg_opres=obs_pres
    usd_idx=where(obs_pres>-990.0)
    lg_obs_pres=array(obs_pres)
    lg_obs_pres[usd_idx]=log10(obs_pres[usd_idx])
    
    
    
    vpl, vpr, vwgt=flb.get_vertical_wgt_1d(lg_mod_pres, lg_grd_pres)

    obs_vpl, obs_vpr, obs_vwgt=flb.get_vertical_wgt_1d(lg_obs_pres, lg_grd_pres)
    
    colwgt=flb.get_col_wgt_1d(lg_grd_pres)
    colwgt=squeeze(colwgt)
    
    
    # print 'do_vert_reverse', do_vert_reverse
    return vpl, vpr, vwgt, obs_vpl, obs_vpr, obs_vwgt, colwgt, do_vert_reverse



def get_col_wv(vpl, vpr, vwgt, grd_colwgt, prof_wv, do_vert_reverse=False):
    """ calculate water vapor column (col_wv)
    from water-vapor profile prof_wv """
    if (do_vert_reverse):
        prof_wv=prof_wv[:,::-1]
    
    grd_wv=flb.prof_vertical_intpl_1d(vpl, vpr, vwgt, prof_wv)
    col_wv=flb.col_int_1d(grd_wv, grd_colwgt)
    
    return col_wv

def get_col_dry_air(grd_colwgt, wv_at_grd):
    """ calculate water vapor column (col_wv)
    from water-vapor profile prof_wv """
    rmass=gc.mh2o/gc.mg
    prof_dry_air=ones(shape(wv_at_grd), float)
    used_idx=where(wv_at_grd<-990.0)
    prof_dry_air[used_idx]=-999.0
    
    
    
    prof_dry_air=flb.array_divide_array_2d(prof_dry_air, wv_at_grd, 1.0)
    col_dry_air=flb.col_int_1d(prof_dry_air, grd_colwgt)
    
    return col_dry_air



def get_xgp0(obs_apr, \
             obs_ak,\
             obs_vpl,\
             obs_vpr,\
             obs_vwgt,\
             wv_at_grd,\
             grd_colwgt):
    

    no_use_idx=where(obs_apr<>obs_apr)
    
    obs_ak[no_use_idx]=-999.0
    obs_apr[no_use_idx]=-999.0
    obs_apr=where(obs_apr<-990.0, -999.0, obs_apr)
    
    
    new_ak=flb.refill_bad_points_1d(obs_ak)
    new_apr=flb.refill_bad_points_1d(obs_apr)

    # project apr to cb grid
    
    apr_at_grd=flb.prof_vertical_intpl_wv(obs_vpl, obs_vpr, obs_vwgt, wv_at_grd, new_apr)
    ak_at_grd=flb.prof_vertical_intpl_1d(obs_vpl, obs_vpr, obs_vwgt, new_ak)
    res_ak_at_grd=1.0-ak_at_grd
    rmass=gc.mh2o/gc.mg
    wv_org=array(wv_at_grd)
    usd_idx=where(wv_org>-990.0)
    wv_org[usd_idx]=1.0-wv_org[usd_idx]/rmass
    
    apr_at_grd=flb.array_divide_array_2d(apr_at_grd, wv_org, 0.0)
    
    
    xgp0=flb.ak_col_int_1d(grd_colwgt, res_ak_at_grd, apr_at_grd)
    return  ak_at_grd, xgp0

def get_xgp_em(grd_pres,
            prof,\
            wv_at_grd,\
            ak_at_grd, \
            vpl, \
            vpr,\
            vwgt,\
            grd_colwgt,\
            mod_offset=0.0,\
            do_vert_reverse=False, debug=False):
    
    """ calculate xgp
    
    Notes:
    1) this version use column water vapor instead of the profile, which could be problematic
 
    2) xgp does not include a-priori
    
    """
    
    # calculate xgp at observation locations
        
    # co2 profile after removing model_offset
    
    prof=prof-mod_offset
    
    if (do_vert_reverse):
        print shape(prof)
        
        prof=prof[:,:, ::-1]
        
    
    prof_at_grd=flb.prof_vertical_intpl_em(vpl, vpr, vwgt, wv_at_grd, prof)
    xgp=flb.ak_col_int_em(grd_colwgt, ak_at_grd, prof_at_grd)
    
    return xgp



    

def get_xgp(grd_pres,
            prof,\
            wv_at_grd,\
            ak_at_grd, \
            vpl, \
            vpr,\
            vwgt,\
            grd_colwgt,\
            mod_offset=0.0,\
            do_vert_reverse=False, debug=False):
    
    """ calculate xgp
    
    Notes:
    1) this version use column water vapor instead of the profile, which could be problematic
 
    2) xgp does not include a-priori
    
    """
    
    # calculate xgp at observation locations
        
    # co2 profile after removing model_offset
    
    prof=prof-mod_offset
    
    if (do_vert_reverse):
        prof=prof[:,::-1]
        
    # print prof[42,:]
    # wv_at_grd[:,:]=0.0
    
    prof_at_grd=flb.prof_vertical_intpl_wv(vpl, vpr, vwgt, wv_at_grd, prof)
    # print prof_at_grd[42,:]
    
    xgp=flb.ak_col_int_1d(grd_colwgt, ak_at_grd, prof_at_grd)
    
    return xgp



def setup_daily_field_em(yyyy, mm, dd, \
                         olon, olat,\
                         gp_list,\
                         em_st_list, em_end_list,\
                         em_yst_list,\
                         em_step_list,\
                         datapath, \
                         gpname='CO2',\
                         prefix="ts_satellite"):
    
    
    """ set up  the model fields from time-series files and ctm files
    
    yyyy, mm, dd            ----in----    year, month, day
    olat, olon              ----in----    obsvervation location
    
    em_st_list, em_end_list ----in---  list of the output tagged tracers.
    The name for time series files is in the format of ts.EN[emst]-EN[emend].[date].bpch
    
    
    datapath  ------ in -----------   directory of the model outputs
    
    gp_list ----- in/out  ----- collection of the tagged tracers
    
    mod_pres  ---out--- model pressure
    mod_wv  ---out --- model water vapor
    
    nem ----- out -----  the number of tagged tracers
    
    """
    
    sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
    nset=len(em_st_list)

    do_sp_wv=False
    
        
    for iset  in range(nset):
        
        em_st=em_st_list[iset]
        em_end=em_end_list[iset]
        cur_yyyy=em_yst_list[iset]
        em_step=em_step_list[iset]

        
        syyyy=r'%4.4d' % cur_yyyy
        sem=r'ST%3.3d.EN%4.4d-EN%4.4d' % (em_step, em_st, em_end)
        # sem=r'EN%4.4d-EN%4.4d' % (em_st, em_end)
        full_datapath=datapath
        
        full_flnm=full_datapath+"/"+syyyy+"/"+prefix+"."+sem+"."+sdate+".bpch"
        
        ftraceinfo=full_datapath+"/"+syyyy+"/"+"tracerinfo"+"."+sem+".dat"
        fdiaginfo=full_datapath+"/"+syyyy+"/"+"diaginfo"+"."+sem+".dat"
        
        bpch2_ts=brw.bpch2_file_rw(full_flnm, "r", \
                                   do_read=1,  ftracerinfo=ftraceinfo,\
                                   fdiaginfo=fdiaginfo)
        
        # print bpch2_ts.tranames
        categorys=None
        if (1==-999):
            for bpdata_in in bpch2_ts.data: 
                print bpdata_in.category, bpdata_in.ntracer
        
        tracers=None
        taus=None
        
        tranames=['CO2']
        print full_flnm
        
        print 'READ CO2'
        
        data_list=None
        print categorys,tracers,taus,tranames, 'am'
        data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
        print 'found co2 data', founded, len(data_list)
        
            
        if (iset==0):

            bpdata=data_list[0]
            
            rlat=bpdata.grid.get_lat()
            rlon=bpdata.grid.get_lon()
            rz=bpdata.grid.get_z()
            
            ax_lon=gax.gp_axis('lon', rlon)
            ax_lat=gax.gp_axis('lat', rlat)
            ax_lz=gax.gp_axis('level', rz)

            lonp1, lonp2, lonwgt=ax_lon.getwgt(olon)
            latp1, latp2, latwgt=ax_lat.getwgt(olat)

            w1=latwgt*lonwgt
            w2=(1.0-latwgt)*lonwgt
            w3=latwgt*(1.-lonwgt)
            w4=(1-latwgt)*(1.0-lonwgt)
            # categorys=["PS-PTOP"]
            
            categorys=None
            tracers=None
            tranames=None
            
            tranames=['PSURF']
            print 'READ Pressure'
            
            ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
            
            # PS
            if (len(ax_data_list)==0):
                
                categorys=["PEDGE-$"]
                tracers=None
                tranames=None
                
                # tranames=['PSURF']
                print 'READ Pressure  again using category '
            
                ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
                # print 'problem'
                # zzzz=raw_input()
                
            gp=ax_data_list[0].data
            gp=squeeze(gp[:,:,0])
            
            
            ps=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            
            # mod pressure
              
            levels=ax_lz[:]
            levels=array(levels)
            levels=levels+1
            levels=levels.astype(int)
            
            pres=pm.get_pres(ps, levels)
            
            mod_pres=array(pres)

            ps=squeeze(ps)
            
            # WV
            tranames=['AVGW']
            
            categorys=None
            tracers=None
            # tranames=None
            ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
            gp=ax_data_list[0].data

            
            
            mod_wv=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            rmass=gc.mh2o/gc.mg
            mod_wv=rmass*mod_wv
            
            # print 'mod wv', mod_wv[0,0:10]
            
            
        # add co2 staff
        
        
        ibp=0
        bp=data_list[0]
        gp=bp.data
        prof0=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)


        icut=0
        if (em_st>1):
            # the last one is not useful
            icut=1
        ndata=len(data_list)
        
        for idata in range(1, ndata-icut):
            bp=data_list[idata]
            gp=bp.data
            prof=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            
            prof=prof-prof0
            gp_list.append(prof)
        
        
        
                
    # we now have 3-D model pressure (mod_pres), model water vapor profile (mod_wv)
    # 3D CO2 (gp_list)
    mod_wv[:,:]=1.0e-8
    
    return mod_pres, mod_wv, ps, gp_list


def setup_daily_field(yyyy, mm, dd, \
                      olon, olat,\
                      gp_list,\
                      em_st_list, em_end_list,\
                      datapath, \
                      gpname='CO2',\
                      prefix="ts_satellite"):
    
    
    """ set up  the model fields from time-series files and ctm files
    
    yyyy, mm, dd            ----in----    year, month, day
    olat, olon              ----in----    obsvervation location
    
    em_st_list, em_end_list ----in---  list of the output tagged tracers.
    The name for time series files is in the format of ts.EN[emst]-EN[emend].[date].bpch
    
    
    datapath  ------ in -----------   directory of the model outputs
    
    gp_list ----- in/out  ----- collection of the tagged tracers
    
    mod_pres  ---out--- model pressure
    mod_wv  ---out --- model water vapor
    
    nem ----- out -----  the number of tagged tracers
    
    """
    
    sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
    nset=len(em_st_list)

    do_sp_wv=False
    
        
    for iset  in range(nset):
        em_st=em_st_list[iset]
        em_end=em_end_list[iset]

        
        syyyy=r'%4.4d' % yyyy
        
        sem=r'ST%3.3d.EN%4.4d-EN%4.4d' % (1, em_st, em_end)
        full_datapath=datapath
        
        full_flnm=full_datapath+"/"+prefix+"."+sem+"."+sdate+".bpch"

        ftraceinfo=full_datapath+"/"+"tracerinfo"+"."+sem+".dat"
        fdiaginfo=full_datapath+"/"+"diaginfo"+"."+sem+".dat"
        
        bpch2_ts=brw.bpch2_file_rw(full_flnm, "r", \
                                   do_read=1,  ftracerinfo=ftraceinfo,\
                                   fdiaginfo=fdiaginfo)
        
        # print bpch2_ts.tranames
        categorys=None
        if (1==-999):
            for bpdata_in in bpch2_ts.data: 
                print bpdata_in.category, bpdata_in.ntracer
        
        tracers=None
        taus=None
        
        tranames=['CO2']
        print full_flnm
        
        print 'READ CO2'
        
        data_list=None
        print categorys,tracers,taus,tranames, 'am'
        data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
        print 'found co2 data', founded, len(data_list)
        
            
        if (iset==0):

            bpdata=data_list[0]
            
            rlat=bpdata.grid.get_lat()
            rlon=bpdata.grid.get_lon()
            rz=bpdata.grid.get_z()
            
            ax_lon=gax.gp_axis('lon', rlon)
            ax_lat=gax.gp_axis('lat', rlat)
            ax_lz=gax.gp_axis('level', rz)

            lonp1, lonp2, lonwgt=ax_lon.getwgt(olon)
            latp1, latp2, latwgt=ax_lat.getwgt(olat)

            w1=latwgt*lonwgt
            w2=(1.0-latwgt)*lonwgt
            w3=latwgt*(1.-lonwgt)
            w4=(1-latwgt)*(1.0-lonwgt)
            # categorys=["PS-PTOP"]
            
            categorys=None
            tracers=None
            tranames=None
            
            tranames=['PSURF']
            print 'READ Pressure'
            
            ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
            
            # PS
            if (len(ax_data_list)==0):
                
                categorys=["PEDGE-$"]
                tracers=None
                tranames=None
                
                # tranames=['PSURF']
                print 'READ Pressure  again using category '
            
                ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
                # print 'problem'
                # zzzz=raw_input()
                
            gp=ax_data_list[0].data
            gp=squeeze(gp[:,:,0])
            
            
            ps=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            
            # mod pressure
              
            levels=ax_lz[:]
            levels=array(levels)
            levels=levels+1
            levels=levels.astype(int)
            
            pres=pm.get_pres(ps, levels)
            
            mod_pres=array(pres)

            ps=squeeze(ps)
            
            # WV
            tranames=['AVGW']
            
            categorys=None
            tracers=None
            # tranames=None
            ax_data_list, founded=bpch2_ts.get_data(categorys, tracers, taus, tranames)
            gp=ax_data_list[0].data

            
            
            mod_wv=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            rmass=gc.mh2o/gc.mg
            mod_wv=rmass*mod_wv
            
            # print 'mod wv', mod_wv[0,0:10]
            
            
        # add co2 staff
        
        
        ibp=0
        for bp in data_list:
            gp=bp.data
            if (ibp==-1000):
                add_bp=data_list[1]
                add_data=add_bp.data -100.0e-6
                
                gp=gp+add_data
    
            prof=grd_pf.grid_prof(gp,lonp1,lonp2,latp1,latp2,w1,w2, w3, w4)
            # to the real dry-air mixing ratio
            
            # prof=(1.0+mod_wv)*prof
            gp_list.append(prof)
            ibp=ibp+1
    
        
                
    # we now have 3-D model pressure (mod_pres), model water vapor profile (mod_wv)
    # 3D CO2 (gp_list)
    mod_wv[:,:]=1.0e-8
    
    return mod_pres, mod_wv, ps, gp_list

    

def xco2jacob(yyyy, mm, dd, \
              olon, olat, obs_pres, \
              opsurf,\
              obs_apr, obs_ak,\
              em_st_list, em_end_list, \
              em_yst_list,\
              em_step_list,\
              mod_offset=0.0,\
              datapath=gcdf.data_path,\
              do_debug_jco2=False):
    
    """
    read jacobian at observation locations from model outputs
    """
       
    gp_list=list()
    
    mod_pres, mod_wv,  mod_sp, gp_list=setup_daily_field_em(yyyy, mm, dd, \
                                                            olon, olat,\
                                                            gp_list,\
                                                            em_st_list, em_end_list,\
                                                            em_yst_list,\
                                                            em_step_list,\
                                                            datapath)
    
    if (mod_pres[0,2]>mod_pres[0,6]):
        do_vert_reverse=True
        new_mod_pres=mod_pres[:, ::-1]
    else:
        new_mod_pres=mod_pres

     
    
    no_use_idx=where((obs_ak<>obs_ak) | (obs_apr<>obs_apr))
    
    obs_ak[no_use_idx]=-999.0
    obs_apr[no_use_idx]=-999.0
    # remove the problematic top level
    
    obs_ak[:,0]=obs_ak[:,1]
    # comment out by lf
    opsurf=mod_sp
    
    
   
    
    # obs_ak just used to indicate which levels are valid

    grd_pres, obs_pres_out, obs_ak_out=flb.combine_vertical_grid(new_mod_pres, obs_pres, obs_ak, opsurf)
    
    obs_pres=obs_pres_out
    obs_ak=array(obs_ak_out)
    
    
    illegal_idx=where(grd_pres<>grd_pres)
    grd_pres[illegal_idx]=-999.0
    
    
    # print 'obs 0', grd_pres[0,:]
    # print 'obs 1', grd_pres[1,:]
    
    # coefficients for vertical interpolation
    
    vpl, vpr, vwgt, \
         obs_vpl, obs_vpr, obs_vwgt, \
         grd_colwgt, do_vert_reverse=setup_vert_intpl(grd_pres, mod_pres, obs_pres)
    
    # water vapor to grid box 
    
    if (do_vert_reverse):
        prof_wv=mod_wv[:,::-1]
    else:
        prof_wv=mod_wv
        
    wv_at_grd=flb.prof_vertical_intpl_1d(vpl, vpr, vwgt, prof_wv)
    col_dry_air=get_col_dry_air(grd_colwgt, wv_at_grd)
    
    # set to be zeros 

    # wv_at_grd[:,:]=1.e-8
    
    # 
    # a-prior at grid box
    
    ak_at_grd, xgp0=get_xgp0(obs_apr, \
                             obs_ak,\
                             obs_vpl,\
                             obs_vpr,\
                             obs_vwgt,\
                             wv_at_grd,\
                             grd_colwgt)
    
    
    
    
    nem=len(gp_list)    
    
    prof=array(gp_list)
    
    
    xgp=get_xgp_em(grd_pres,\
                   prof,\
                   wv_at_grd,\
                   ak_at_grd, \
                   vpl, \
                   vpr,\
                   vwgt,\
                   grd_colwgt,\
                   mod_offset=0.0,\
                   do_vert_reverse=do_vert_reverse, debug=False)
    
    
    xgp=flb.array_divide_array_2d_1d(xgp, col_dry_air, 0.0)
    
    prof_h=transpose(xgp)
    prof_h=1.0e6*prof_h # changed to ppm
    
    return prof_h, mod_sp



def xco2jacob_sel(yyyy, mm, dd, \
                  olon, olat, obs_pres,
                  opsurf,\
                  obs_apr,\
                  obs_ak,\
                  em_st_list, em_end_list, \
                  sel_em=1,\
                  mod_offset=0.0,\
                  datapath=gcdf.data_path,do_debug_jco2=False):
    
    
    
    
    """
    read model xgp values at observation locations
    
    """


    gp_list=list()
    
    mod_pres, mod_wv, mod_sp, gp_list=setup_daily_field(yyyy, mm, dd, \
                                                olon, olat,\
                                                gp_list,\
                                                em_st_list, em_end_list,\
                                                datapath)
    
    if (mod_pres[0,2]>mod_pres[0,6]):
        do_vert_reverse=True
        new_mod_pres=mod_pres[:, ::-1]
    else:
        new_mod_pres=mod_pres

    
    # no_use_idx=    
    no_use_idx=where(obs_ak<-990.0)
    # print 'obs_ak'
    # print obs_ak[184,:]
    
    obs_ak[no_use_idx]=-999.0
    obs_apr[no_use_idx]=-999.0
    no_use_idx=where(obs_apr<-990.0)
    obs_ak[no_use_idx]=-999.0
    obs_apr[no_use_idx]=-999.0
    # obs_ak[:,:]=1.0
    obs_ak[:,0]=obs_ak[:,1]
    
    # obs_ak just used to indicate which levels are valid
    print shape(mod_sp)

    print shape(opsurf)
    
    tsurf=opsurf-mod_sp
    print max(tsurf), min(tsurf)
    # print mod_pres[0,:]
    # comment outby lf
    
    opsurf=mod_sp
    

    # print 'obs_ak', obs_ak[0,:]
    
    # print obs_pres[0, :]

    # print opsurf[0]
    
    
    grd_pres, obs_pres_out, obs_ak_out=flb.combine_vertical_grid(new_mod_pres, obs_pres, obs_ak, opsurf)
    obs_pres=array(obs_pres_out)
    obs_ak=array(obs_ak_out)
    
    # print 'obs_ak 2'
    # print obs_ak[184,:]
    # print obs_pres[184,:]
    # obs_ak[:,:]=1.0
    
    # print grd_pres[184,:]
    
    illegal_idx=where(grd_pres<>grd_pres)
    grd_pres[illegal_idx]=-999.0
    
    # print 'obs 0', grd_pres[0,:]
    # print 'obs 1', grd_pres[1,:]
    
    # coefficients for vertical interpolation
    
    vpl, vpr, vwgt, \
         obs_vpl, obs_vpr, obs_vwgt, \
         grd_colwgt, do_vert_reverse=setup_vert_intpl(grd_pres, mod_pres, obs_pres)
    
    # water vapor to grid box 
    
    if (do_vert_reverse):
        prof_wv=mod_wv[:,::-1]
    else:
        prof_wv=mod_wv
        
    wv_at_grd=flb.prof_vertical_intpl_1d(vpl, vpr, vwgt, prof_wv)
    wv_at_grd=squeeze(wv_at_grd)
    # set to be zeros 
    col_dry_air=get_col_dry_air(grd_colwgt, wv_at_grd)
    # wv_at_grd[:,:]=1.e-8
    
    # 
    illegal_idx=where(wv_at_grd<>wv_at_grd)
    illegal_idx=squeeze(illegal_idx)
    if (len(illegal_idx)>0):
        print 'illegal value'
        
        wv_at_grd[illegal_idx]=-999.0
    else:
        print 'wv_at_grd', max(wv_at_grd.flat), min(wv_at_grd.flat)
    
    print shape(wv_at_grd)
    
    # a-prior at grid box

    ak_at_grd, xgp0=get_xgp0(obs_apr, \
                             obs_ak,\
                             obs_vpl,\
                             obs_vpr,\
                             obs_vwgt,\
                             wv_at_grd,\
                             grd_colwgt)
    
    
    
    
    prof=gp_list[sel_em-1]
    # print 'prof', prof[0,:]
    # print 'prior', obs_apr[0,:]
    # print 'xgp0', xgp0[0:5]
    
    
    xgp=get_xgp(grd_pres,\
                prof,\
                wv_at_grd,\
                ak_at_grd, \
                vpl, \
                vpr,\
                vwgt,\
                grd_colwgt,\
                mod_offset=0.0,\
                do_vert_reverse=do_vert_reverse, debug=False)
    
    xgp=xgp+xgp0
    xgp=flb.array_divide_array(xgp, col_dry_air, 0.0)
    
    xgp=1.0e6*xgp
    
    inter_idx=where(xgp<315.0)
    inter_idx=squeeze(inter_idx)
    # inter_idx=[184,44]
    
    nidx=size(inter_idx)
    if (nidx>1):
        print 'small values'
        itx=inter_idx[0]
        itx=108
        
        print itx
        print 'prof'
        print  prof[itx,:]
        print 'pressure'
        print grd_pres[itx,:]
        print 'xgp0'
        
        print 1.0e6*xgp0[itx]
        print 'xgp'
        print xgp[itx]
        print 'new mod pres'
        print new_mod_pres[itx, :]
        print 'obs pres'
        print obs_pres[itx, :]
        print 'opsurf'
        print opsurf[itx]
        
        ddd=raw_input()
    return xgp, mod_sp





if (__name__=='__main__'):
    import os
    import satellite_operator as sop
    em_st_list=[1]
    em_end_list=[2]
    obs_datapath='./gosat_v28_obs/'
    viewtype='gosat_v28'
    viewmod='nadir'
    data_path='./enkf_rerun/'
    yyyy=2009
    doy_st=192
    doy_end=doy_st+8
    op=sop.satellite_xgp(viewtype, obs_datapath, 0, True, True)
    for doy in range(doy_st, doy_end):
        op.read_obs(yyyy, doy)
        yyyy, mm, dd=tmdl.doy_to_time_array(doy, yyyy)
        xgp=op.read_sel_xgp(yyyy, mm, dd,\
                            1, 2,\
                            1,\
                            data_path)
        
        
         
        
               
    
           

