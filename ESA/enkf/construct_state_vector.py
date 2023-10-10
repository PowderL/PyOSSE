"""
   this code is used to construct state vector 
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.06.28
   History: v0.95, 2013.02.15
   
   Functions:
   
   =================================================
   1. construct_step_tag: Construct tag to identify ensemble run steps
   2. add_step_prior: Adding apriori for one step into the state vector
   3. hm_subtract_m1: subtracting reference values for response function

   
   


"""
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_descfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio
import ESA.enkf.run_desc_file_m as rdescio
import ESA.enkf.run_desc_m as rdesc
import ESA.enkf.state_vector as stv_m
import ESA.enkf.x2flux_m as x2flux


def construct_step_tag(yyyy, mm, dd, step, fmt, **keywords):
    """
    Construct tag to identify ensemble run steps 
    
    Inputs:
    ---------------------------------------------------
    1. yyyy, mm, dd:<int>: date
    2. step: <int>: ensemble run step 
    3. fmt:<str>: format of the tag 
    4. keywords:<dict>: extra inputs
    

    Returns:
    --------------------------------------------
    1. tag:<str>: tag used to identify step 
    
    """
    path=""
    s_step=r'%3.3d' % (step)
    all_keywords={'XYYYYX':yyyy, 'XSTEPX':s_step}
    all_keywords.update(keywords)
    tag=ovio.construct_filename(path, fmt, **all_keywords)
    return tag


def add_step_prior(yyyy, mm, dd, \
                       cl_run_desc, cl_stv, \
                       cl_x2flux=None,\
                       tag_fmt='XYYYYX.STXSTEPX',\
                       ftag=construct_step_tag, \
                       fscale_m=None,\
                       **keywords):
    
    
    """
    Adding apriori for one step into the state vector
    
    Inputs:
    -------------------------------------------------

    1. yyyy, mm, dd:<int>: date
    2. cl_run_desc:<run_desc_cl>: class for ensemble run description 
    3. cl_stv: <state_vector>: class for state vector 
    4. cl_x2flux:<x2flux_cl>:  class for project x to fluxes 
    5. tag_fmt:<str>:  format for tag used to identify format 
    6. ftag   :<func>: function for generating tag 
    7. fscale_m:<func>: for generating correlation 
    8. keywords:<dict>: extra inputs
    

    Returns:
    --------------------------------------------------
    1. iadd_x:<int>: number of state vector added 
    
    """
    
    # S1: check coverage on given date
    
    nset=cl_run_desc.extract_coverage(yyyy, mm, dd)
    
    
    tag_lst=[]
    cur_tag=""
    ems_st=9999
    ems_end=0
    tshift=0
    
    iadd_x=0
    
    # S2: check each set to find if a new x (period) to be added 
    
    for iset in range(nset): 
        
        # #T1: starting time 
        
        istep=cl_run_desc.step_lst[iset]
        yyyy_st=cl_run_desc.st_year_lst[iset]
        yyyy_st=int(yyyy_st)
        doy_st=cl_run_desc.st_doy_lst[iset]
        doy_st=int(doy_st)
        
        
        yyyy_st, mm_st, dd_st=tm.doy_to_time_array(doy_st, yyyy_st)
        tau_st=tm.get_tau(yyyy_st, mm_st, dd_st)

        # #T2: mode for indexing ensemble members
                
        mod=cl_run_desc.idx_mod_lst[iset]
        
        # #T3:ensemble member in the run 
        
        enr_st=cl_run_desc.st_no_lst[iset]
        enr_end=cl_run_desc.end_no_lst[iset]
        
        
       
        # #T4: end time 

        yyyy_end=cl_run_desc.end_year_lst[iset]
        yyyy_end=int(yyyy_end)
        
        doy_end=cl_run_desc.end_doy_lst[iset]
        doy_end=int(doy_end)
        
        yyyy_end, mm_end, dd_end=tm.doy_to_time_array(doy_end, yyyy_end)
        tau_end=tm.get_tau(yyyy_end, mm_end, dd_end)

        

        
        # #T5: generate Tag
        
        new_tag=ftag(yyyy, mm, dd, istep, fmt=tag_fmt, **keywords)

        print 'set:',         iset
        print 'start time:', yyyy_st, mm_st, dd_st
        print 'end time:', yyyy_end, mm_end, dd_end
        print 'st end no:',  enr_st, enr_end
        print 'tag:',        new_tag
        print 'mod:',        mod
        
        
        
        # #T6: check whether tag is already included in the class cl_stv
        # ##: if not, a new member is added in
        
        
        if (cur_tag==""):
            # #c: just start 
            
            if (mod==0):
                # there is 1 redundant member in the run 
                tshift=tshift+1
        
            ems_st=enr_st     # ensemble emission file start 
            ems_end=enr_end   # emission file end 
            cur_tag=new_tag   
            cur_yyyy, cur_mm, cur_dd=yyyy_st, mm_st, dd_st

        elif (cur_tag==new_tag):
            
            # #c: the same set of ensemble emission
            
            if (mod==0):
                # there is 1 redundant member in the run 
                tshift=tshift+1
        
            if (ems_st>enr_st):
                ems_st=enr_st
                
            if (ems_end<enr_end):
                ems_end=enr_end

            # ##: starting time for emissions
            
            cur_yyyy, cur_mm, cur_dd=yyyy_st, mm_st, dd_st

        else:
            # #c:  move to a new tag 
            
            # ##d:  add aprior for previous sets if necessary
            
            
            if (not cur_tag in cl_stv.step_tag_lst):
                # ##:  check  size
                
                ems_end=ems_end-tshift
                nx=ems_end-ems_st+1
                ne=ems_end-ems_st+1
                
                # ##: make a copy of class x2flux_cl (and move dates)
                
                if (cl_x2flux<>None):
                    new_cl_x2flux=cl_x2flux.copy(cur_yyyy, cur_mm, cur_dd, \
                                                     ems_st=ems_st, ems_end=ems_end)
                else:
                    
                    new_cl_x2flux=cl_x2flux
                

                # ##: check scale_m is necessary 
                scale_m=None
                if (fscale_m<>None):
                    scale_m=fscale_m(x, new_cl_x2flux)
                    
                # ##: create default aprior for this time period
                
                
                x, dx=cl_stv.vf_construct_aprior(nx, ne)
                
                # ##: add 
                
                cl_stv.add_new_x_to_window(x, \
                                               dx, \
                                               tau_st,\
                                               tau_end,\
                                               step_tag=cur_tag,\
                                               enlarge_factor=1.0,\
                                               scale_m=scale_m, \
                                               lag_window=None,\
                                               cl_x2flux=cl_x2flux)
                
                
                iadd_x=iadd_x+1
                

                # #c: set data to iset 

                ems_st=9999
                ems_end=0
                tshift=0
                
                
                if (mod==0):
                    # there is 1 redundant member in the run 
                    tshift=tshift+1
                
                ems_st=enr_st
                ems_end=enr_end
                cur_tag=new_tag
                cur_yyyy, cur_mm, cur_dd=yyyy_st, mm_st, dd_st
                
    
        
    # ##c: 
    if (ems_end>0):
        ems_end=ems_end-tshift
        # ##c:
        if (not cur_tag in cl_stv.step_tag_lst):
            
            nx=ems_end-ems_st+1
            ne=ems_end-ems_st+1
            
            if (cl_x2flux<>None):
                new_cl_x2flux=cl_x2flux.copy(cur_yyyy, cur_mm, cur_dd, ems_st=ems_st, ems_end=ems_end)
            else:
                new_cl_x2flux=cl_x2flux
                 
                 
            scale_m=None
            if (fscale_m<>None):
                scale_m=fscale_m(new_cl_x2flux)

                
            x, dx=cl_stv.vf_construct_aprior(nx, ne)
                
            cl_stv.add_new_x_to_window(x, \
                                           dx, \
                                           tau_st,\
                                           tau_end,\
                                           step_tag=cur_tag,\
                                            enlarge_factor=1.0,\
                                           scale_m=scale_m, \
                                           lag_window=None,\
                                           cl_x2flux=cl_x2flux)
             
            iadd_x=iadd_x+1
             
    return iadd_x

def hm_subtract_m1( enr_yst_lst,\
                        enr_dst_lst,\
                        enr_step_lst,\
                        enr_em_st_lst,\
                        enr_em_end_lst, hm):
    
    new_hm=list()
    nset=len(enr_em_st_lst)
    ipos=0
    for iset  in range(nset):
        
        enr_em_st=enr_em_st_lst[iset]
        enr_em_end=enr_em_end_lst[iset]
        enr_yyyy=enr_yst_lst[iset]
        enr_step=enr_step_lst[iset]
        enr_doy=enr_dst_lst[iset]
        icut=0
        print 'iset, ipos', iset, ipos
        
        if (enr_em_st>1):
            # the last one is not useful
            icut=1
        else:
            print 'set vh0'
        
            vh0=hm[:, ipos]
        
        ipos=ipos+1
        nenr=enr_em_end-enr_em_st+1
        for usd_idx in range(1, nenr-icut):
            vh=hm[:,ipos]
            vh=vh-vh0
            new_hm.append(vh)
            ipos=ipos+1
        
        ipos=ipos+icut
    
    new_hm=npy.array(new_hm)
    new_hm=npy.transpose(new_hm)
    return new_hm



#<<< TEST >>> 

if (__name__=='__main__'):

    # build class  for flux
    
    datapath='/scratch/local/otool_data/surface_flux/'
    flnm='CO2_EMISSION.XYYYYXDXDOYX.nc'
    yyyy, mm, dd=2009, 3,1
    
    # create class
    
    
    cl_x2flux=x2flux.x2flux_cl(datapath, flnm, \
                                   yyyy, mm, dd, \
                                   bf_get=x2flux.read_bf, \
                                   bf_cl=None, \
                                   lon_nm='longitude', \
                                   lat_nm='latitude', \
                                   flux_nm='flux')
    
    
    

    
    
    # build class for run_desc
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    flnm='enr_cfg_XYYYYX.dat'

    
    
    cl_run_desc=rdesc.enr_desc_cl(flnm,\
                                      datapath, \
                                      yyyy, mm, dd,\
                                      fopen=rdescio.open_enr_desc_file,\
                                      fread=rdescio.read_enr_desc_file,\
                                      fclose=rdescio.close_enr_desc_file,\
                                      fget=rdescio.get_enr_desc_table)
    
    
    # build class for state vector 
    
    max_step=20
    lag_window=5
    
    cl_stv= stv_m.ens_stat_cl(yyyy,\
                                  mm, \
                                  dd,\
                                  max_step,\
                                  lag_window,\
                                  mod_bias=[],\
                                  mod_bias_err=[],\
                                  use_tcor=False,\
                                  tcor_len=30.0*24.00,\
                                  tcor_factor=1.0)
    
    
    add_step_prior(yyyy, mm, dd, \
                       cl_run_desc, cl_stv, \
                       cl_x2flux=cl_x2flux,\
                       tag_fmt='XYYYYX.STXSTEPX',\
                       ftag=construct_step_tag, \
                       fscale_m=None)
    
    
    print cl_stv.wnd_nx, cl_stv.wnd_ne
    
