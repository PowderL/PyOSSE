#! /geos/u23/epd7/bin/python2.7
"""
this program is used to pre-calculate the response functions for the pre-defined perturbations (basis functions)


"""
# includes modules  
import gc_input_gen as gig      #  generate the input file for geos-chem model  
import sys                        # system commands 
import os                         #  operational system commands
import restart_gen as rg          # regenerate the restart file 
import ESA.util.time_module as tm          #  time conversion  
import time as systime 
import geos_chem_def as gcdf      # deflaut settings for job running 

import read_em_cfg as emcfg       # configurations for pre-define perturbation functions
import numpy as npy 




print '*'*80
print ''*30+'CO2 ENSEMBLE RUN DRIVER'+'*'*30
print '*'*80
print ' '*30

# Setting from job definition file 

## T1:  read in configuration for ensemble run

re_run_now=gcdf.rerun_now
## T2:restart date

re_run_date=gcdf.rerun_date

## T3: starting year/mm/dd

yyyy=gcdf.st_yyyy  # the year  
mm=gcdf.st_mm     # the month 
dd=gcdf.st_dd     #  the day

## T4 temporal resolution 

temp_res=gcdf.temporal_resolution # temporal resolution 
# Notes: exact temporal resolutions is controlled by ensemble emission file
 
timestep=temp_res*24.0*3600.0 #  temporal resolutio in seconds

## T5:  steps;  notes: nstep=ntime
 
ntime=gcdf.nstep #  the last day of geos-chem simulation is given by ntime * temp_res 

nstep=ntime


## lag windows


lag_window=gcdf.inv_lag_window


## T5: switch for  initial file (restart file)

new_restart=gcdf.new_restart

## T6: ensemble emission configuration file (template)


ems_cfg_flnm=gcdf.cfg_ems_flnm

## read in starting 1) do 2) yyyyy 3) number of ensemble flux functions and emission file names

ems_cdf_d=gcdf.cfg_ems_open(gcdf.ems_cfg_path, 
                                yyyy, \
                                mm,\
                                dd)

                                

data=get_ems_cfg_data(fdesc, yyyy, mm, dd)


ems_doy_st_lst=data['doy']
ems_yyyy_st_lst=data['year']
ems_npb_lst=data['npb']
ems_flnm_lst=data['flnm']



print em_yyyy_st_lst
print em_npb_lst

fem=

# jjjj=raw_input()


##/c/desk clock for jobs starting

gmt=systime.gmtime()

##/c/output file for configurations to be used in inversion
ftt=open(gcdf.data_path+"ens_pos.dat", "w")

##/c/the head lines  

line='geos_chem run at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' % (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])
print line
ftt.write(line+'\n')
line=r'temp_res: %4.4d  nstep: %4.4d' % (temp_res, ntime)
ftt.write(line+'\n')
##/c/time step/member start/member end/year start/year end/doy start/doy end/BF file/extension of model output/method_for_select_member_for_inversion/
line='step,  mem_st,  mem_end,  year_st,  year_end,  doy_st,  doy_end,  co2flnm, output_name, index_method'
ftt.write(line+'\n')


#>/ps/do the ensemble simulations

maxtracer=80

select_runs=range(0, nstep)
default_restart_file=gcdf.restart_file

#>>/ps/loop0/over each step

for istep in range(0, nstep): # run through the emission time period 
    
    
    npbf=em_npb_lst[istep]
    pbflnm=em_flnm_lst[istep]

    yyyy_st=em_yyyy_st_lst[istep]
    doy_st=em_doy_st_lst[istep]
    yyyy, mm, dd=tm.doy_to_time_array(doy_st, yyyy_st)
    ##c/tau_st in seconds
    tau_st=tm.get_tau(yyyy,mm, dd)

    if (npbf>maxtracer):
        # separate to different batch of jobs
        
        #notes: only the second one in use and the position in CO2_mod. for pb_flux is calculated as
        # PBFLUX_ID=PBST+N-1. N starting from 2 to mend-mst+1
        
        last_mem=0
        a_mst=[1]
        a_mend=[maxtracer+1]
        nset=1
        last_mem=a_mend[nset]
        while (last_mem<npbf):
            a_mst=a_mst+[last_mem] #
            new_last_mem=min(npbf+nset+1, last_mem+maxtracer+1)
            
            a_mend=a_mem+[new_last_mem]
            last_mem=new_last_mem
            nset=nset+1
            
        
    else:
        a_mst=[1]
        a_mend=[npbf+1]
    
    nsub_run=len(a_mend)
    
    
   
    
    print  '-'*10+'starting year, month, day:',  yyyy, mm, dd
    
    istep_end=min(istep+lag_window, len(em_doy_st_lst)-1)
    
    #>>>/ps/loop0/loop1/sub_run to cover every segments of ensemble members

    ##/c/1/relative positition of the first used perturbation ensemble
    ##/c/2/(i.e, the second tracer) in the pb functions
    
    
    pbst=0
    
    
    for isub_run  in range(nsub_run): 
        
        mst=a_mst[isub_run]
        mend=a_mend[isub_run]
        
        pbst=mst-1
        
        line=r'%3.3d, %4.4d, %4.4d, %4.4d, %3.3d, %3.3d, %3.3d' % (istep, \
                                                                   mst, mend, \
                                                                   em_yyyy_st_lst[istep], em_yyyy_st_lst[istep_end],
                                                                   em_doy_st_lst[istep], \
                                                                   em_doy_st[istep_end])
        print line
        print '======>Step 2: Generate input file<======'
        
    
        if (gcdf.make_restart):
            # # c/generate file extension for restart files
        
            enaf=gcdf.make_extension(yyyy, nn, dd, istep, mst, mend)
            # #c finalize the information line for outputs/0 means the default indexing method

            
            
        line=line+ ','+pbflnm+','+enaf+',0\n'
        print line
        ftt.write(line+'\n')
        
        
        
        #  run over each days between istep to istep_end
        
        istep_shift=0
        for end_date in em_doy_st[istep:istep_end]:
            
            yst=em_yyyy_st_lst[istep+istep_shift]
            dst=em_doy_st_lst[istep+istep_shift]
            yyyy, mm, dd=tm.doy_to_time_array(dst, yst)
            tau0=tm.get_tau(yyyy,mm, dd)
            tau0=tau0/3600.0
            
            tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
            tst=tst.replace('-', '')
            tst=tst.replace(':', '')
            time_start=r'%4.4d%3.3d' % (yst, dst)
            print '======>Step 3: Generate restart file<====='
            
            full_restart_name='restart.'+enaf+'.'+tst[0:8]+'00'
            
            ##c/bydefault, pbuse is equal to one
            pbuse=ones(2)
            print full_restart_name
            ##c/generate configuration file for ensemble run
            
            igg.create_new_input_file(gcdf.run_path,\
                                          gcdf.data_path,\
                                          istep,\
                                          em_yyyy_st_lst[istep+istep_shift],\
                                          em_doy_st_lst[istep+istep_shift],\
                                          em_yyyy_st_lst[istep+istep_shift],\
                                          em_doy_st_lst[istep+istep_shift],\
                                          pbuse,\
                                          pbst,\
                                          member_start, \
                                          member_end,\
                                          pbflnm,\
                                          em_yyyy_st_lst[istep+istep_shift],\
                                          em_doy_st_lst[istep+istep_shift],\
                                          em_yyyy_st_lst[istep+istep_shift],\
                                          em_doy_st_lst[istep+istep_shift])
            
            os.system("mv input.geos.new input.geos")
            
            print '======>Step 3: Generate restart file<====='

            #>>>>/ps/generate restart file
            

            ntracers=mend-mst+1
                                
            ##/c/restart file names
            
            full_restart_name='restart.'+enaf+'.'+tst[0:8]+'00'

            ##/c/the number of traces to be duplicated  in the restart file
            real_ntracers=ntracers
            
            if (istep_shift==0):
                ##/c/1/at the beginning of the  time line for each PBF, mixing ratio of all tracers
                ##/c/2/should be fixed to a constant or the same 3D distribution   
                
                resflnm=default_restart_file
                if (new_restart):
                    do_cp=False
                else:
                    os_cp_cmd='cp '+resflnm+' '+gcdf.run_path+'/'+full_restart_name
                    do_cp=True
            else:
                # copy  restart files from data path to rerun path
                
                resflnm='restart.'+enaf+'.'+tst[0:8]+'00'
                resflnm=gcdf.data_path+'/'+resflnm
                os_cp_cmd='cp '+resflnm+' '+gcdf.run_path+'/'+full_restart_name
                print os_cp_cmd
                do_cp=True
            
            
            
                
            
            # launch the CTM 
            
            if (istep in gcdf.select_runs):
                if (re_run_now or tst>=re_run_date):
                    if (do_cp):
                        os.system(os_cp_cmd)
                    else:
                        rsf=rg.geos_chem_restart_file(resflnm)
                        rsf.mod_restart_file(1, full_restart_name, real_ntracers, \
                                                 fixed_value=gcdf.fixed_init_val, do_regrid=False,  \
                                             keep_first=False)
                    
                            
                
                
                    if (gcdf.run_ctm):
                        os.system('sh ./rungeos.sh')
            
            istep_shift=istep_shift+1
            
        #>>>>/ps/loop0/loop1/loop2/end
        ##/c/1/for each segment, the first one is used as references, no
        ##/c/2/corresponding in the
        
        pbst=mend-1
    #>>>/ps/loop0/loop1/end

#>>/ps/loop0/end


 
    
 
    
##c/ close ensemble run file
ftt.close()


        


