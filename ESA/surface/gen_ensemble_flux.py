"""
Functions for generating pulse-like ensemble fluxes.  

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.21
    
    
    Functions:
    ===============================================
    1. flux_pb_normalization: normalize regional flux to a given unit
    
    2. gen_ensemble_flux_by_reg: generate ensemble flux 
    
"""
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.otool_var_io as ovio
import ESA.util.geo_constant as gc
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.flux_regrid as fgrd

import bf_m as bfcl

import numpy as npy
import pylab as plb

#<<<parameters>>>

mol_mass_factor=gc.mc/gc.An # in gram 

# convertion factor for the total flux from Mol/s*m^2 to PgC/yr




# default setting 

import ESA.util.otool_menu_m as menu_m


def process_configuration_file(cfg_flnm):
    """ process configuration menu 
    Inputs:
    ----------------------------------------------
    1. cfg_flnm:<str>: name of the configuration file 
    
    Returns:
    ------------------------------------------------
    1. root_menu: <menu_cl>: menu class for configuration 
    
    """
    
    if ('.xml' in cfg_flnm):
        # #c: if the menu is given in xml  file 
        root_menu=menu_m.xml_to_menu(cfg_flnm)
    else:
        # #c: if the menu is given in format text file 
        root_menu=menu_m.txt_to_menu(cfg_flnm)
    
    root_menu.print_menu()

    return root_menu

# #c: show definition 



def flux_pb_normalization(reg_flux, reg_map, sf_area, \
                              unit_flux):
    
    """
    normalize regional fluxes to a given unit flux 
      
    Inputs:
    ----------------------------------------
    1. reg_flux:<array, (nlon, nlat [,nlayer])>: regional flux
    2. reg_map:<array, (nlon, nlat [,nlayer])>: region map
    3. sf_area:<array, (nlon, nlat)>: area for grid boxes
    4. unit_flux:<float>: unit for regional fluxes to be normalized. 
    
    Returns:
    ------------------------------------
    1. reg_flux:<array>: normalized regional fluxes. 
    
    
    """
    # S1: check whether it is multiple layer map 
    
    dims=npy.shape(reg_flux)
    
    if (len(dims)>2):
        # #c: mult-layer map
        
        nlayer=npy.size(reg_flux)
        
        # S2: loop over layer 
        
        for ilayer in range(nlayer):
            
            sel_flux=reg_flux[:, :, ilayer]
            sum_flux=npy.sum(sel_flux*sf_area)
            if (abs(sum_flux)>0):
                ratio=unit_flux/sum_flux
                sel_flux=ratio*sel_flux
                reg_flux[:, :, ilayer]=sel_flux
            else:
                # no flux over this region 
                reg_flux[:, :, ilayer]=0.0
        # loop ilayer end
        
        return reg_flux
        
    else:
        
        # #c: single-layer map 
        
        nlayer=int(npy.max(reg_map.flat))
        
        if (nlayer==0):
            nlayer==1

        total_flux=npy.zeros(npy.shape(reg_flux), float)
        # #S3: loop over region
        
        for ireg in range(1, nlayer+1):
            sel_flux=npy.where(reg_map==ireg,reg_flux, 0)
            sum_flux=npy.sum(sel_flux*sf_area)
        
            if (abs(sum_flux)>0):
                ratio=unit_flux/sum_flux
                sel_flux=ratio*sel_flux
            
            total_flux=total_flux+sel_flux
        # loop ireg end
        
        return total_flux


def gen_ensemble_flux_by_reg(cfg):
    
    """ 
    generate basis functions by time 
    
    Inputs:
    ----------------------------
    1. cfg:<menu>: the configuration files 
    
    
    """
    
    
    
    # S1: time steps for ensemble flux 
    
    
    
    yyyy_st=cfg['sample.year']
    mm_st=cfg['sample.month']
    dd_st=cfg['sample.day']
    
    nstep=cfg['sample.nstep']
    time_step=cfg['sample.time']
    
    
    ts_keywords=cfg['sample.step_keywords']
    
    # #c: function to generate time 
    ens_step_f=cfg['sample.fstep']
    
    # #c: generate ensemble list for starting time
    
    
    ens_doy_lst, ens_yyyy_lst=ens_step_f(yyyy_st, mm_st, dd_st, time_step, nstep, **ts_keywords)
    
    print ens_doy_lst

    
    # S2: generate ensemble fluxes for each time step 
    
    # #c: file to be saved 
    
    desc_flnm=""  # description file output for ensemble fluxes  
    fl_desc=None
    
    
    for istep in range(0, nstep+1):
        # #c: current doy
        
        cur_doy=ens_doy_lst[istep]
        cur_yyyy=ens_yyyy_lst[istep]
    
        # #c: 
        
        next_doy=ens_doy_lst[istep+1]
        next_yyyy=ens_yyyy_lst[istep+1]
        
        doy=cur_doy
        
        yyyy, mm, dd=tm.doy_to_time_array(cur_doy, cur_yyyy)

        # S3: get Basis functions for current time period
        
        fio_keys=cfg['bf.keywords']
        
        print 'fio_keys', fio_keys
        
        
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
        
        
        
       # #c: longitude

        lon=cur_bf.fdesc.get_data('lon')
        lon=npy.array(lon)
        
        # #c: latitude
        
        lat=cur_bf.fdesc.get_data('lat')
        lat=npy.array(lat)
        
        # #c: regional map 
        
        reg_map=cur_bf.fdesc.get_data('map')
        reg_map=npy.array(reg_map)
        
        # #c: regional flux (orginal basis functions)

        reg_flux=cur_bf.fdesc.get_data('flux')
        
        reg_flux=npy.array(reg_flux)
        
        # #c: area for grid box
                                       
        area=cur_bf.fdesc.get_data('area')
        area=npy.array(area)
        
        # #c: read or re-set parent ID 
        
        bf_varname_lst=cfg['bf.varname_lst']
        
        if ('pid' in bf_varname_lst):
            
            # #c: if parent id has been read in 
            pid=cur_bf.fdesc.get_data('pid')
            pid=npy.array(pid)
            
        else:
            # ### set parent ID as one 
            
            dims=npy.shape(reg_map)
            
            if (npy.size(dims)>2):
                
                nlayer=dims[2]
                pid=npy.ones(nlayer)
                
            else:
                
                nlayer=int(npy.max(reg_map))
                if (nlayer==0):
                    nlayer=1
                    
                pid=npy.ones(nlayer)
        
        
                              
        
                                       
        # #S4: sampling original flux function to generate ensemble fluxes. 
        
        sample_f=cfg['sample.fsample']
        sample_keywords=cfg['sample.keywords']
        
        print 'map:', max(reg_map.flat), min(reg_map.flat)
        print 'flux:', max(reg_flux.flat), min(reg_flux.flat)
        
        ens_flux_dict=sample_f(lon, lat, reg_flux, \
                                   reg_map, \
                                   pid, \
                                   **sample_keywords)
        
        # #c: regional map for ensemble fluxes. 
        
        
        reg_map=ens_flux_dict['map'] 
        new_flux=ens_flux_dict['flux'] 
        
        print 'map2 :', max(reg_map.flat), min(reg_map.flat)
        print 'flux2 (kgC/m2/s):', max(reg_flux.flat), min(reg_flux.flat)
        
        # ttt=raw_input()
        
        # S5: write ensemble emission configuration file 
        
        # #c: layer (number) 
        
        nlayer=npy.size(reg_map[0,0,:])
        
        
        # #c: name for description file for saved ensemble fluxes. 
        # #c: this file will be used by EnKF to generate state vector etc 
        
        ens_path=cfg['ens.desc_path']
        ens_flnm=cfg['ens.desc_flnm']
        
        new_desc_flnm=ovio.construct_filename(ens_path, ens_flnm, \
                                                 XYYYYX=yyyy, XMMX=mm, XDDX=dd, \
                                                  XDOYX=doy)
        
        # #c: check whether a new description file for ensemble  outputs need to be  open 
        
        
        print new_desc_flnm
        
        if (new_desc_flnm<>desc_flnm):
            
            # if a new file is open 
            
            if (fl_desc<>None):
                fl_desc.close()
            
            desc_flnm=new_desc_flnm
            fl_desc=open(new_desc_flnm, 'w')
            line_0=r'# nstep:%3.3d; temporal_resolution:%3.3d d; npb_flux:%4.4d' % (nstep, time_step, nlayer)
            fl_desc.write(line_0+'\n')
            line_1='# year, doy,  nmem, emission_flnm'
            fl_desc.write(line_1+'\n')
            
        
        # #c: Short file name (without path) for sampled fluxes 
            
        
        
        new_ens_file=ovio.construct_filename("", cfg['ens.flnm'], \
                                                 XYYYYX=yyyy, \
                                                 XMMX=mm, XDDX=dd, \
                                                 XDOYX=doy)
        
        
        # #c: write description on flux to fl_desc, which includes 
        # #c: current year ;  current day of year ; layer number ; ensemble emission file name
        
        line_d=r'%4.4d, %3.3d, %4.4d, %-30s' % (cur_yyyy, cur_doy, nlayer, new_ens_file)
        fl_desc.write(line_d+'\n')
        
        
        # #S6: write ensemble flux to file 
        
        # #c: emission file name 
        
        new_ens_flnm=ovio.construct_filename(cfg['ens.path'], cfg['ens.flnm'], \
                                                 XYYYYX=yyyy, \
                                                 XMMX=mm, XDDX=dd, \
                                                 XDOYX=doy)
                
        # ##c: list of dimensions for flux ensemble 
        
        dim_name_lst=cfg['ens.dim_name_lst']
        dim_type_lst=cfg['ens.dim_type_lst']
        
        # ##c: use dimension variable constructor to construct var
        
        vdim_constructor=cfg['ens.fdim']
        
        dim_lst=vdim_constructor(dim_name_lst, dim_type_lst)
        
        # fill dimension data 
        # if the data is not in ens_flux_dict, the bf one will be used
        for dvar in dim_lst:
            
            dname=dvar.name
            
            if (dname in ens_flux_dict):
                data=ens_flux_dict[dname]
            else:
                data=cur_bf.fdesc.get_data(dname)
            
            data=npy.array(data)
                
            
            dvar.set_data(data)
        
        # ##c: list of variables to be write 
        
        vname_lst=cfg['ens.vname_lst']
        vdim_lst=cfg['ens.vdim_lst']
        vtype_lst=cfg['ens.vtype_lst']
        var_constructor=cfg['ens.fvar']

        var_lst=var_constructor(vname_lst, vdim_lst, vtype_lst)
        
        # # fill data
        
        for fvar in var_lst:
            
            vname=fvar.name
            # if vname is not in ens_flux_dict, it will be filled with 
            # name from dict
            
            if (vname in ens_flux_dict):
                data=ens_flux_dict[vname]
            else:
                data=cur_bf.fdesc.get_data(vname)
            
                
            data=npy.array(data)
            fvar.set_data(data)
        
        
                
        
        # #T9: write out the emission files
        
        ncfio.ncf_save_var(new_ens_flnm,  var_lst, dim_lst,  create_new=True)    
        
    if (fl_desc<>None):
        fl_desc.close()


#<<< TESTS >>>>

if __name__=='__main__':
    
    cfg_flnm='flux_ensemble_def.cfg'
    
    root_cfg=process_configuration_file(cfg_flnm)
    gen_ensemble_flux_by_reg(root_cfg)
    
    
