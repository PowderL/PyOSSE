""" the settings for generate an ensemble by sampling  
flux perturbations (i.e, the basis functions) at given time steps. 

It have 5 section: 


"""

import bf_m as bfm
import bf_file_m as bf_ncio
import numpy as npy
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio
import ESA.util.otool_ncfile_io as ncfio
import sample_reg_flux as sample_bf
import ESA.util.geo_constant as gc

class flux_ensemble_def_cl:
    def __init__(self,  root_menu):
        """
        Initialization
        """

        # PATH and file name 
        self.datapath=root_menu['flux_file.path']
        self.flnm=root_menu['flux_file.flnm']
        
 
        # #T2: Variables to be read
        
        
        self.varname_lst=root_menu['flux_file.varname_lst']
        self.varname_dict=root_menu['flux_file.varname_lst']

        # FUNCTIONS
        

        self.fopen=root_menu['flux_file.fopen']
        self.fread=root_menu['flux_file.fopen']
        self.fclose=root_menu['flux_file.fclose']
        self.fget=root_menu['flux_file.fget']
        self.keywords=root_menu['flux_file.keywords']
        
        # S3: flux  sampling procedure. 
        
        self.fsample=root_menu['sample.fsample']
        
        self.sel_reg_lst=root_menu['sample.sel_reg_lst']
        



    def sel_make_iovar(self):
        # constructure IO VAR list 
        em_lon_nm='lon'  # longitude name used in bf file  
        var_lon= ovio.io_var_cl(em_lon_nm,'f',  \
                                    [em_lon_nm], None)

        # #lat

        em_lat_nm='lat'   # latitude name 
        var_lat= ovio.io_var_cl(em_lat_nm,'f',  \
                                    [em_lat_nm], \
                                    None)
        
        # #layer

        em_layer_nm='layer'    # layer 
        var_layer= ovio.io_var_cl(em_layer_nm,'i',  \
                                      [em_layer_nm], \
                                      None)
        
        # #map

        em_map_nm='map'   # map 
        var_map= ovio.io_var_cl(em_map_nm,'f',  \
                                    [em_lon_nm, em_lat_nm, em_layer_nm], \
                                    None)
        
        # #flux

        em_flux_nm='flux'  # flux 
        var_flux= ovio.io_var_cl(em_flux_nm,'f',  \
                                     [em_lon_nm, em_lat_nm, em_layer_nm], \
                                     None)

        # #area
        em_area_nm='area'  
        var_area= ovio.io_var_cl(em_area_nm,'f',  \
                                     [em_lon_nm, em_lat_nm], \
                             None)

        # #pid

        em_pid_nm='pid'

        var_pid= ovio.io_var_cl(em_pid_nm,'i', [em_layer_nm], \
                            None)
        

        # #selected reg

        em_sel_nm='sel_reg_lst'

        var_sel= ovio.io_var_cl(em_sel_nm,'i',  [em_layer_nm], None)

        # #define dimension in the file 

        self.em_dim_lst=[var_lon,\
                             var_lat,\
                             var_layer]
        
        # #define variable list 

        self.em_var_lst=[var_map,\
                             var_flux, \
                             var_area,  \
                             var_pid, \
                             var_sel]
        
        # #path for emission ensemble file 
        
        em_path='./'

        # #output files for ensemble fluxes at different time periods  

        em_flnm_flux='CO2.EMISSION.XYYYYXDXDOYX.nc'

        # output configuration file
        
        em_flnm_cfg='cfg_em_reg_144.XYYYYX'

        
#<SECTION (SETTING)>

#  temporal resolution 
        
        em_time_resolution=32
        
        # number of time step 
        
        em_time_nstep=30  # two years 

    # True to move the day to the start of next month when necessary

        em_reset_at_month_st=True

# True to move the day to the start of next year  when necessary
        
        em_reset_at_year_st=True

# renormalize the each basis function to 1 GtC /per year
        
        em_flux_renorm=True
# ##c scaling factor. re-scale to GtC/yr
        
        
        time_in_second=24.*3600.
        mol_mass_factor=gc.mc/gc.An 
        flux_scale_factor=mol_mass_factor*time_in_second
# 1.0e-15: g to Pg; 1.0e4: m2 to cm2 (for area); 365: day of year
        flux_scale_factor=1.0e-15*1.0e4*365*flux_scale_factor  
        
        em_st_yyyy, em_st_mm, em_st_dd=2009,1,1


# S4: list for doys and years when the flux pulse start at each step 

#<SECTION (INITIALIZATION)>

em_doy_list=list()
em_yyyy_list=list()


# S3: the following codes for setup time steps 

#  local variables 

cur_dd=em_st_dd
cur_mm=em_st_mm
cur_yyyy=em_st_yyyy

cur_doy=tm.day_of_year(cur_yyyy, cur_mm, cur_dd)

em_doy_list.append(cur_doy)
em_yyyy_list.append(cur_yyyy)
ndays_year=tm.days_in_year(cur_yyyy)

# setup the 
for istep  in range(em_time_nstep+1):
    new_doy=cur_doy+em_time_resolution
    # check the time 
    

    if (new_doy>ndays_year):
        cur_yyyy=cur_yyyy+1

        if (em_reset_at_year_st):
            new_doy=1
        ndays_year=tm.days_in_year(cur_yyyy)

    new_dd, new_mm, new_yyyy=tm.doy_to_time_array(new_doy, cur_yyyy)

    if (new_mm<>cur_mm):
        cur_mm=new_mm
        if (em_reset_at_month_st):
            new_doy=tm.day_of_year(cur_yyyy, cur_mm, 1)
    cur_doy=new_doy
    
    em_doy_list.append(new_doy)
    em_yyyy_list.append(cur_yyyy)
    
    

                                   


