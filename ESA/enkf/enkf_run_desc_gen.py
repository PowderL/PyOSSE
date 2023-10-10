"""
this program is to generate configurations and drive tagged simulations 

"""

# includes modules  
import sys                           # system commands 
import os                            #  operational system commands
import ESA.util.time_module as tm    #  time conversion  
import time as systime   
import enkf_def as emdf      # default settings or forward (tagged) simulations 
import read_em_cfg as emcfg          # configurations for pre-define perturbation functions
import numpy import npy 

# S1 starting time 

yyyy=emdf.st_yyyy  # the year  
mm=emdf.st_mm     # the month 
dd=emdf.st_dd     #  the day

# S2 temporal resolution (the exact temporal resolutions is controled by the emission Basis function)

temp_res=emdf.temporal_resolution # temporal resolution 
ntime=emdf.nstep #  


# S3  lag window

lag_window=emdf.inv_lag_window



# S5 perturbation configuration files (BF)

em_cfg=emdf.cfg_pb_file

## read in list for tracer, starting doy yyyy, number of BF functions, and filenames


ems_traname_lst, ems_doy_st_lst, ems_yyyy_st_lst,  ems_npb_lst, ems_flnm_lst=emcfg.read_em_conf(em_cfg)

print em_doy_st
print em_yyyy_st
print em_npb



##/c/desk clock for jobs starting

gmt=systime.gmtime()

## output file for configurations to be used in inversion

ftt=open(emdf.data_path+'ens_cfg.'+str(st_yyyy)+'.dat', "w")

## headers  

line='geos_chem run at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' % (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])
print line
ftt.write(line+'\n')
line=r'temp_res: %4.4d, nstep: %4.4d' % (temp_res, ntime)
ftt.write(line+'\n')
## column names: step; member start; member end; year start ; year end; doy start; doy end; BF file; extensions for model output;  tracer; method_for_select_member_for_inversion; 


line='step,  mem_st,  mem_end,  year_st,  year_end,  doy_st,  doy_end,  ems_flnm, ctm_flnm, tracer; index_method'

ftt.write(line+'\n')


#>/ps/do the ensemble simulations

maxtracer=80

select_runs=range(0, nstep)
default_restart_file=emdf.restart_file

#>>/ps/loop0/over each step

for istep in range(0, nstep): # run through the emission time period 
    npbf=ems_npb[istep]
    pbflnm=ems_flnm[istep]
    
    yyyy_st=ems_yyyy_st[istep]
    doy_st=ems_doy_st[istep]
    
    if (npbf>maxtracer):
        #notes: only the second one in use and the position in CO2_mod. for pb_flux is calculated as
        # PBFLUX_ID=PBST+N-1. N starting from 2 to mend-mst+1
        
        
        a_mst=[1, maxtracer+1] #  
        a_mend=[maxtracer+1,npbf+2]
    else:
        a_mst=[1]
        a_mend=[npbf+1]
        
    nsub_run=len(a_mst)
    
    
    
    print  '-'*10+'starting year, month, day:',  yyyy, mm, dd
    
    istep_end=min(istep+lag_window, len(em_doy_st)-1)
    
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
                                                                   em_yyyy_st[istep], em_yyyy_st[istep_end],
                                                                   em_doy_st[istep], \
                                                                   em_doy_st[istep_end])
        print line
        print '======>Step 2: Generate input file<======'

    
        ##c/print time_start
        ##c/generate file extension for restart files
        
        enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (istep, mst, mend)
        ##c finalize the information line for outputs/0 means the default indexing method

        line=line+ ','+pbflnm+','+enaf+',0\n'
        print line
        ftt.write(line+'\n')
        
        
        
        #>>>>/ps/loop0/loop1/loop2/run over each days till at time of istep_end
        
        istep_shift=0
        for end_date in em_doy_st[istep:istep_end]:
            
            yst=em_yyyy_st[istep+istep_shift]
            dst=em_doy_st[istep+istep_shift]
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
            
            igg.create_new_input_file(istep,\
                                      em_yyyy_st[istep+istep_shift:istep+istep_shift+2],\
                                      em_doy_st[istep+istep_shift:istep+istep_shift+2],\
                                      pbuse,\
                                      pbst,\
                                      member_start=mst, \
                                      member_end=mend,\
                                      co2flnm=pbflnm,\
                                      time_start=time_start,\
                                      em_doy=em_doy_st[istep+istep_shift:istep+istep_shift+2],\
                                      em_yyyy=em_yyyy_st[istep+istep_shift:istep+istep_shift+2])
            
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
                    os_cp_cmd='cp '+resflnm+' '+'./'+full_restart_name
                    do_cp=True
            else:
                ##/c/1/at the following time, restart files are ouputs from previous 
                resflnm='restart.'+enaf+'.'+tst[0:8]+'00'
                resflnm=emdf.data_path+'/'+resflnm
                os_cp_cmd='cp '+resflnm+' '+'./'+full_restart_name
                print os_cp_cmd
                do_cp=True
            
            
            
                
            
            #>>>>>/ps/loop0/loop1/loop2/lauch the CTM 
            
            if (istep in emdf.select_runs):
                if (re_run_now or tst>=re_run_date):
                    if (do_cp):
                        os.system(os_cp_cmd)
                    else:
                        rsf=rg.geos_chem_restart_file(resflnm)
                        rsf.mod_restart_file(1, full_restart_name, real_ntracers, \
                                             fixed_value=emdf.fixed_init_val, do_regrid=False,  \
                                             keep_first=False)
                    
                            
                
                
            
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


        


