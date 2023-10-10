import ESA.util.time_module as tm
# default directories
run_path='./'
data_path=run_path+'/enkf_output/' # the model output

# S1: starting time for geos chem simulation

st_yyyy, st_mm, st_dd=2009,1,1
st_doy=tm.day_of_year(st_yyyy, st_mm, st_dd)

# S2: time resultion in day

temporal_resolution=32  #  the temporal resolution for inversions in days

# S3: lag window

inv_lag_window=5 #   

# S4: number of simulations for pulse-like BFs

nstep=12

# S5: configuration files for BF
# #c: the temporal resolution etc will be re-set by configuration file 

cfg_pb_file='./surface_flux/cfg_em_flux.2009'

# S6: force the system to re-initialize to be a fixed value when necessary 

new_restart=True
fixed_init_val=100.0e-6

# #S7: defaul restart file 

# #c: If new_restart=True, it will just be used as template. 
# #c: All the initial fields will be filled with fixed_init_val 

restart_file='./restart.20090101'  # 

# S8: time step for launching GEOS-Chem 

rerun_now=True
rerun_date='20090101'

# S9: selected run steps

select_runs=range(0, 2) # nstep)


# S10: setting for GEOS-Chem 

# #c: the maximal tracer number of GEOS-Chem 
maxtracer=80

# #c: the name for description file


enr_desc_file="test_ens_pos.dat"
