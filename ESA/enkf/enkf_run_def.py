import ems_cfg_file_m as ems_cfg
import run_cfg_file_m as run_cfg
import numpy as npy 


# setting for run a job
# default directories


# S1 starting time for tagged simulations

st_yyyy, st_mm, st_dd=2010,1,1
st_doy=tm.day_of_year(st_yyyy, st_mm, st_dd)

# S2  Temporal resultion (in days)

temporal_resolution=32  #  the temporal resolution for inversions in days

# S3  total steps of the simulations
### the whole period=temporal_resolution X nstep

nstep=12

# S4 lag window (in steps)

inv_lag_window=5 #  


#SECTION < EMISSIONS FILES > 

# ##configuration file access
ems_cfg_path='../surface/'
ems_cfg_file='cfg_em_reg_144.XYYYYX'

ems_cfg_fopen=ems_cfg.open_ems_cfg_file
ems_cfg_fclose=ems_cfg.close_ems_cfg_file
ems_cfg_fread=ems_cfg.read_ems_cfg_read
ems_cfg_fget=ems_cfg.get_ems_cfg_data

# ##T: column information 

ems_cfg_colnm_lst=ems_cfg.ems_cfg_colnm_lst
ems_cfg_coltype_lst=ems_cfg.ems_cfg_coltype_lst
ems_cfg_name_dict=ems_cfg.ems_cfg_colnm_dict
ems_cfg_io_keywords={}


#SECTION (run configuration file)

# ##configuration file access

run_cfg_path='../surface/'
run_cfg_flnm='ens_pos.XYYYYX'
run_cfg_fopen=rcfg_f.open_ens_cfg_file
run_cfg_fclose=rcfg_f.close_ens_cfg_file
run_cfg_open_write=rcfg_f.open_ens_cfg_write
run_cfg_write_head=rcfg_f.write_ens_colum_heads
run_cfg_write_data=rcfg_f.write_ens_colum_val

# ##T1: column names in the configuration 

run_cfg_colnm_lst=ecfg_f.ens_cfg_column_lst
run_cfg_coltype_lst=ecfg_f.ens_cfg_coltype_lst
run_cfg_name_dict=ens_colnm_dict
run_cfg_io_keywords={}






maxtracer=80


#SECTION (RUN)
run_path='./'   # executive codes
run_code='rungeos.sh' 

data_path='./'  # the model output

#SECTION (RESTART)

restart_flnm='restart.ENXEN_STX-ENXEN_ENDX.XYYYYXXMMXXDDX'
restart_fopen=rg.open_restart_file
restart_fread=rg.close_restart_file
restart_fupdate=rg.reset_restart_file 

#SECTION (CTM INPUTS)









