#! /geos/u23/epd7/bin/python2.7
"""
this program is used to pre-calculate 
the response functions for the pre-defined perturbations (basis functions)
This program will read in configurations on flux perturbation (BF) ensemble, and then
1): generate input.geos
2): generate restart file
3): launch geos_chem (rungeos.sh)


Notes:
=================================================================
1. Run schedule 
Ensemble Run is first made from the starting time T0 when the  first (file) set of mission 
perturbation BFs occur till T0+lag_window. 
Second run will be made from T0+1 (step) to T0+1+lag_window for the second (file) set of BFs. 
Third run  will cover T0+2 till T0+1+lag_window
 ... 
However,  when the size of each BF set is larger than maximal tracer number 
of GEOS-Chem CTM. The runs will be split into sub-runs. 

2. Emission perturbations as being used to force simulation
=========================================================================
First tracer will not be forced by any emissions. 
In co2_mod.f, tracer of number N corresponds to perturbations of BF ID: 
PBFLUX_ID=PBST+N-1. 
When no subrun is made. The total tracer number N=SIZE(BF)+1
When subruns are made. The total tracer number from all subruns is equal to  IZE(BF)+number of sub-runs

"""

# includes modules  
import input_geos_gen as igg      #  generate the input file for geos-chem model  
import sys                        # system commands 
import os                         #  operational system commands
import restart_gen as rg          # regenerate the restart file 
import ESA.util.time_module as tm          #  time conversion  
import time as systime 
import geos_chem_def as gcdf      # deflaut settings
import read_em_cfg as emcfg       # configurations for pre-define perturbation functions
from numpy import *




print '*'*80
print ''*30+'CO2 ENSEMBLE RUN DRIVER'+'*'*30
print '*'*80
print ' '*30


# S1: read in configuration for ensemble runs

# #c: restart now 
re_run_now=gcdf.rerun_now

# #c: restart date

re_run_date=gcdf.rerun_date

# #c: starting time for ensemble run (which may be different from re-run date)

yyyy=gcdf.st_yyyy  # the year  
mm=gcdf.st_mm     # the month 
dd=gcdf.st_dd     #  the day

# #c: time resolution  
# #c: exact temporal resolutions is controled by the emission basis function (BF) 

temp_res=gcdf.temporal_resolution # temporal resolution 
timestep=temp_res*24.0*3600.0 #  temporal resolutio in seconds
ntime=gcdf.nstep #  the last day of geos-chem simulation is given by ntime * temp_res 


# #: lag window

lag_window=gcdf.inv_lag_window

# #c: run steps

nstep=ntime

# #c: T/F for reset the initial distribitribution to the fixed value

new_restart=gcdf.new_restart


# #c: decription file on flux basis functions (BF) files

em_cfg=gcdf.cfg_pb_file

# S2: read in starting doy/yyyy/number of pb functions/stored files
em_doy_st, em_yyyy_st, em_npb, em_co2flnm=emcfg.read_em_conf(em_cfg)

print 'BF starting doys:', em_doy_st
print 'BF starting years:', em_yyyy_st
print  'Number of BF files:',  em_npb

# jjjj=raw_input()



# S3: output file for configurations to be used in inversion

ftt=open(gcdf.data_path+gcdf.enr_desc_file, "w")

# #c: desk clock for jobs starting

gmt=systime.gmtime()

# #c: the head lines in configuration file 


line='geos_chem run at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' % (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])
print line
ftt.write(line+'\n')

line=r'temp_res: %4.4d  nstep: %4.4d' % (temp_res, ntime)
ftt.write(line+'\n')
# #c: column heads
# # columns include /time step/member start/member end/year start/year end/doy start/doy end/BF file/extension of model output/method_for_select_member_for_inversion/

line='step,  mem_st,  mem_end,  year_st,  year_end,  doy_st,  doy_end,  co2flnm, output_name, index_method'
ftt.write(line+'\n')


# S4: do the ensemble forward simulations

# #:c 

maxtracer=gcdf.maxtracer


select_runs=range(0, nstep)
default_restart_file=gcdf.restart_file


# loop: loop over each step
      

for istep in range(0, nstep): 
    # #c: starting time 
    yyyy_st=em_yyyy_st[istep]
    doy_st=em_doy_st[istep]
    yyyy, mm, dd=tm.doy_to_time_array(doy_st, yyyy_st)

    # ##: tau_st in seconds
    tau_st=tm.get_tau(yyyy,mm, dd)

    # #c: BF set (file)
    npbf=em_npb[istep]
    pbflnm=em_co2flnm[istep]


    if (npbf>maxtracer):
        
        a_mst=[1, maxtracer+1] #  
        a_mend=[maxtracer+1,npbf+2]
    else:
        a_mst=[1]
        a_mend=[npbf+1]
        
    nsub_run=len(a_mst)
    
    
    
    print  '-'*10+'starting year, month, day:',  yyyy, mm, dd
    
    istep_end=min(istep+lag_window, len(em_doy_st)-1)
    
    # #T1: loop for sub_run to cover every segments of ensemble members
    
    # #c: positition of the first used perturbation ensemble 
    
    pbst=0
    
    
    
    
    
    for isub_run  in range(nsub_run): 
        mst=a_mst[isub_run]
        mend=a_mend[isub_run]
        pbst=mst-1
        
        line=r'%3.3d, %4.4d, %4.4d, %4.4d, %3.3d, %3.3d, %3.3d' % (istep, \
                                                                       mst, mend, \
                                                                       em_yyyy_st[istep], \
                                                                       em_yyyy_st[istep_end],\
                                                                       em_doy_st[istep], \
                                                                       em_doy_st[istep_end])
        
        print line
        print '======>Step 2: Generate input file<======'

    
        # #c: print time_start
        # #c: generate file extension for restart files
        
        enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (istep, mst, mend)
        # #c: finalize the information line for outputs/0 means the default indexing method

        line=line+ ','+pbflnm+','+enaf+',0\n'
        print line
        ftt.write(line+'\n')
        
        
        
        # #T2: loop run over each starting days at em_doy_st[istep] till at time of istep_end
        
        # #c: run is made from step: istep
        # #c: istep_shift denotes the run steps 

        istep_shift=0
        # loop: em_doy_st

        for end_date in em_doy_st[istep:istep_end]:
            
            yst=em_yyyy_st[istep+istep_shift]
            dst=em_doy_st[istep+istep_shift]
            yyyy, mm, dd=tm.doy_to_time_array(dst, yst)
            tau0=tm.get_tau(yyyy,mm, dd)
            tau0=tau0/3600.0
            
            tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
            tst=tst.replace('-', '')
            tst=tst.replace(':', '')
            time_start=tst[0:10]
            
            
            print '======>Step 3: Generate input.geos file<====='
            
            full_restart_name='restart.'+enaf+'.'+tst[0:8]+'00'
            
            # #c: bydefault, pbuse is equal to one
            pbuse_lst=ones(2)
            print full_restart_name
            # #T3: input.geos for CTM ensemble run
            
            igg.create_new_input_file(istep,\
                                          gcdf.run_path,\
                                          gcdf.data_path,\
                                          em_yyyy_st[istep+istep_shift:istep+istep_shift+2],\
                                          em_doy_st[istep+istep_shift:istep+istep_shift+2],\
                                          pbuse_lst,\
                                          pbst,\
                                          member_start=mst, \
                                          member_end=mend,\
                                          pbflnm=pbflnm,\
                                          time_start=time_start,\
                                          em_doy=em_doy_st[istep+istep_shift:istep+istep_shift+2],\
                                          em_yyyy=em_yyyy_st[istep+istep_shift:istep+istep_shift+2], \
                                          do_bk_run=4)
            
            os.system("mv input.geos.new input.geos")
            
            print '======>Step 4: Generate restart file<====='

            # #c: generate restart file
            
            
            ntracers=mend-mst+1
                                
            # #T4: restart file 

            # #c: file name 
            
            full_restart_name='restart.'+enaf+'.'+tst[0:8]+'00'
            
            # #c: the number of traces to be duplicated  in the restart file
            real_ntracers=ntracers
            
            if (istep_shift==0):
                # #c: at the first run step for each BF file, 3D distributions of all tagged tracers
                # #c: should be the same, either be fixed to a constant or be the same prior field    
                
                resflnm=default_restart_file

                if (new_restart):
                    do_cp=False
                else:
                    os_cp_cmd='cp '+resflnm+' '+'./'+full_restart_name
                    do_cp=True
            else:
                # #c: At the following time, restart files are ouputs from previous step. 

                resflnm='restart.'+enaf+'.'+tst[0:8]+'00'
                resflnm=gcdf.data_path+'/'+resflnm
                os_cp_cmd='cp '+resflnm+' '+'./'+full_restart_name
                print os_cp_cmd
                do_cp=True
            
            
            
            print '======>Step 5: launch GEOS-Chem <====='

            # #T4: Launch the CTM 
            
            if (istep in gcdf.select_runs):
                if (re_run_now or tst>=re_run_date):
                    # #c: generate or copy restart file
                    
                    if (do_cp):
                        os.system(os_cp_cmd)
                    else:
                        rsf=rg.geos_chem_restart_file_cl(resflnm)
                        
                        rsf.create_enr_restart_file(full_restart_name, \
                                                        real_ntracers, \
                                                        newtau=None, \
                                                        filling_value=gcdf.fixed_init_val, \
                                                        sel_idx=1,\
                                                        use_selected_tracer=False)
                        
        
                          
                
                    
                    # #c: run GEOS-Chem 
                    
                    os.system('sh ./rungeos.sh')
                
            istep_shift=istep_shift+1
            
        # loop: em_doy_st end
            
 

        # #c: set the position in BF file for next segement 
            
        pbst=mend-1

    # loop: range(nsub_run) end

# loop: range(0, nstep) end



 
    
 
    
# S7: close ensemble run file
ftt.close()


        


