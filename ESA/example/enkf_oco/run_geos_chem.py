"""
  class used to launch and control GEOS-Chem tagged 'single' tracer forecast calculation 

  Classes:
  ======================================
  1. geos_chem_cl:class used to manage GEOS-Chem simulation 
  
"""
import sys
import os
import numpy as npy
import ESA.util.time_module as tm
import ESA.util.bpch2_rw_smp as brw
import time as systime
import ESA.util.otool_var_io as ovar
    
class geos_chem_cl:
    
    """
    Members:
    ====================================================
    
    1. runpath:<str>: run directory 
    2. datapath:<str>: output data directory 
    3. tra_st: <int>:  starting tracer 
    4. tra_end:<int>:  ending tracer
    5. finput_gen:<func>: function used to generate input.geos
    6. ext:<str>: extension for output file names 
    7. run_type:<int>: type of run 
    8. runscript:<str>: to make restart file 
    10. doy_st:<int>: starting doy
    11. yyyy_st:<int>: starting year
    13. doy_end:<int>: end doy
    14. yyyy_end:<int>: end year
        
    
    Functions:
    -------------------------------------------------
    1. __init__: initialization
    2. set_time: set starting and ending time for forecasting
    3. rerun_geos_chem: run GEOS-Chem by one step
    4. prepare_next_restart_file: make changes to first tracer according to analysis increments
    
    
    """

    def __init__(self, \
                     step,doy_st, yyyy_st, \
                     doy_end, yyyy_end,\
                     runpath, \
                     datapath, \
                     tra_st,\
                     tra_end,\
                     finput_gen,\
                     ext,\
                     run_type,\
                     runscript\
                     ):
    
        """
        Initialization
        
        Inputs:
        ==========================================================

        1. runpath:<str>: run directory 
        2. datapath:<str>: output data directory 
        3. tra_st: <int>:  starting tracer 
        4. tra_end:<int>:  ending tracer
        5. finput_gen:<func>: function used to generate input.geos
        6. ext:<str>: extension for output file names 
        7. run_type:<int>: type of run 
        8. runscript:<str>: to make restart file 
        9. istep:<int>: step number 
        10. doy_st:<int>: starting doy
        11. yyyy_st:<int>: starting year
        13. doy_end:<int>: end doy
        14. yyyy_end:<int>: end year
        

    
        """

        self.runpath=runpath
        self.datapath=datapath
        self.tra_st=tra_st
        self.tra_end=tra_end
        self.finput_gen=finput_gen
        self.ext=ext
        self.run_type=run_type 
        self.runscript=runscript
        # date time 
        
        self.doy_st=doy_st
        self.doy_end=doy_end
        self.step=step

        self.yyyy_st=yyyy_st
        self.yyyy_end=yyyy_end

        
        
    def set_time(self, step,doy_st, yyyy_st, \
                     doy_end, yyyy_end):
        """ set time 
        1. run_type:<int>: type of run 
        2. runscript:<str>: to make restart file 
        3. doy_st:<int>: starting doy
        4. yyyy_st:<int>: starting year
        5. doy_end:<int>: end doy
        6. yyyy_end:<int>: end year

        """
    
        self.doy_st=doy_st
        self.doy_end=doy_end
        self.step=step
        
        self.yyyy_st=yyyy_st
        self.yyyy_end=yyyy_end
        
    
    def rerun_geos_chem(self, step,doy_st, yyyy_st, \
                            doy_end, yyyy_end):
        
        """
        run GEOS-Chem for  given time period 

        Input:
        ========================================
        
        1. step:<int>: step number 
        2. doy_st:<int>: starting doy
        3. yyyy_st:<int>: starting year
        4. doy_end:<int>: end doy
        5. yyyy_end:<int>: end year
        
        
        """
        
        
        self.set_time(step,doy_st, yyyy_st, \
                     doy_end, yyyy_end)
        
        

    
        print '*'*50
        print ''*20+'CO2 ENSEMBLE RUN DRIVER'+'*'*20
        print '*'*50
        print ' '*50
        
        print '======>Step 1: Generate co2 emission data<======'
    
        
        # S1 --------------- simulation period -------------------
        
        
        enr_doy_st=[doy_st, doy_end]
        enr_yyyy_st=[yyyy_st,yyyy_end]
        
        
        yst=yyyy_st
        dst=doy_st
        
        tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
        tst=tst.replace('-', '')
        tst=tst.replace(':', '')
        time_start=tst[0:8]
        
        # S2: fill values (no pb will be used)

        pbuse_lst=[0, 0]
        pbst=1

        

        # S3: generate input.geos

        print self.run_type
        
        self.finput_gen(self.step,\
                            self.runpath,\
                            self.datapath,\
                            enr_yyyy_st,\
                            enr_doy_st,\
                            pbuse_lst,\
                            pbst,\
                            member_start=self.tra_st, \
                            member_end=self.tra_end,\
                            newext=self.ext,\
                            time_start=time_start,\
                            em_doy=enr_doy_st,\
                            em_yyyy=enr_yyyy_st, \
                            do_bk_run=self.run_type)

        os.system('mv input.geos.new input.geos')
        

        # S4: launch GEOS-Chem 
        
        os.system('sh '+self.runscript)
    
    def prepare_next_restart_file(self, newdata, resflnm=None, action=0):
        """
        
        Make change to first tracer of the restart file

        Inputs:
        -------------------------------------------------
        1. newdata:<array, (nlon, nlat, nlvl)>:  changes to first tracer of restart file 
        2. resflnm:<str>: the restart file name 
        3. action:<int>: 0: newdata will be added to concentrations of the first one 
        
        
        """
        # S1: read in original restart file for next time step 
        
        
        if (resflnm==None):
            # prepare restart file for next time step 
            
            tst=tm.doy_to_utc(self.doy_end, sec=0, yyyy=self.yyyy_end)
            tst=tst.replace('-', '')
            tst=tst.replace(':', '')
            full_restart_name='restart.'+self.ext+'.'+tst[0:8]+'00'
        
        else:
            full_restart_name=resflnm
        
        
        newfile=self.runpath+"/"+full_restart_name
        full_restart_name=self.datapath+"/"+full_restart_name
        
        bpch2=brw.bpch2_file_rw(full_restart_name, 'r', do_read=1)
        bplist, found=bpch2.get_data(None, None, None, None)
        bpdata=bplist[0]
        new_data0=bpdata.data
        
        print '--------prepare restart file------' 
        print 'full_restart_name:', full_restart_name
        print ' before change'
        print 'max:', npy.max(new_data0.flat)
        print 'min:', npy.min(new_data0.flat)

        # S2: modify the data
        
        
        if (action==0):
            # #c: add in  
            new_data0=new_data0
        
        elif (action==1):
            # #c: replace
            
            new_data0=new_data0+newdata

        else:
            new_data0=newdata

        print ' after change'
        print 'npy.max:', npy.max(new_data0.flat)
        print 'min:', npy.min(new_data0.flat)

        # S3: save data to restart file 
            
        print 'new restart file:', newfile
        
        bpch2_w=brw.bpch2_file_rw(newfile, 'w')
        ifunit=95
        bpch2_w.open_bpch2_w(ifunit,'test')
        bpdata.data=new_data0
        
        bpdata.ntracer=1
        bpdata.write_to_bpch2_file(ifunit)
        
        itracer=2
        for bpdata in bplist[1:]:
            # second tracer
            bpdata.ntracer=itracer
            bpdata.write_to_bpch2_file(ifunit)
            itracer=itracer+1
       
        
        
        bpch2_w.close_bpch2_file() 
        
        
        

    
    
    def transfer_restart_file(self, enr_st_lst, \
                                  enr_end_lst, \
                                  enr_yst_lst, \
                                  enr_step_lst,\
                                  use_set_lst,\
                                  xinc, \
                                  path,\
                                  flnm\
                                  ):
        
        """
        calculate increments to CO2 concentrations from analysis increment 
        
        
        Inputs:
        =========================================================
        1. enr_st_lst:<list, t:int>: starting ensemble member number
        2. enr_end_lst:<list, t:int>: ending ensemble member number
        3. enr_yst_lst:<list, t:int>: starting year
        4. enr_step_lst:<list, t:int>: starting step
        5. use_set_lst:<list, t:int>: the set data will be used 
        6. xinc:<array, (nx)>: the number of coefficient 
        7. path:<str>: location of the data 
        8. flnm:<str>: names for ensemble restart file 
        Returns:
        ========================================================
        1. new_data:<array, (nlon, nlat, nlvl)>: increment of the CO2 concentrations
        
        
        """
        
        yyyy, doy=self.yyyy_end, self.doy_end
        nset=len(enr_st_lst)
        sdate=r'%4.4dD%3.3d' % (yyyy, doy)
        yyyy, mm,dd=tm.doy_to_time_array(doy, yyyy)
        tst=tm.doy_to_utc(doy, sec=0, yyyy=yyyy)
        tst=tst.replace('-', '')
        tst=tst.replace(':', '')
        # the total number of ensemble 
    
        # #c: ic is the index inside xinc, as CTM outputs includes redundant tracers such as  
        # #c: reference etc, which are not part of the perturbation ensemble. 
        
        ic=0
        
        # S1: loop over each perturbation data set can calculate the results 
        
        for eid in use_set_lst:
            
            # S2: read in restart file from GEOS-Chem ensemble run 
            
            est=enr_st_lst[eid]
            eend=enr_end_lst[eid]
            yyyy=enr_yst_lst[eid]
            enr_step=enr_step_lst[eid]
            print '====prepare_next_restart_file===='
            print 'use ensemble set (step, enr_st, enr_end, year):', enr_step, est, eend, yyyy 
            
            full_flnm=ovar.construct_enr_filename(path, flnm,\
                                                      yyyy, dd, mm,\
                                                      yyyy, \
                                                      doy, \
                                                      enr_step, \
                                                      est, \
                                                      eend)
            
            
            print 'read data from ensemble restart file: ', full_flnm
            
            
            bpch2=brw.bpch2_file_rw(full_flnm, 'r', do_read=1)
            bplist, found=bpch2.get_data(None, None, None, None)
            
            ne=len(bplist)
            
            if (ic==0):
                new_data=npy.zeros(npy.shape(bplist[0].data), float)
            
            # #c: the first is reference data without perturbation 
            
            base_data=npy.array(bplist[0].data)
            
            print 'reference data max, min:', enr_step, npy.max(base_data.flat), npy.min(base_data.flat)    
            
            icut=0
            
            if (est>1):
                # #c: last one is not useful for following set of each step 
                icut=1
                
            # #c: adding perturbation to
            
            # S3:  accumulated contributions from each member included in the restart file 
            for ie in range(1, ne-icut):
                data=bplist[ie].data
                data=data-base_data
                new_data=new_data+data*xinc[ic]
                ic=ic+1
            
            # loop: ie in range(1, ne-icut):  end 
                
        # loop: eid end 
        print 'last ic', ic
        print 'increment to restart file (max, min)', npy.max(new_data.flat), npy.min(new_data.flat)
        print '===prepare_next_restart_file (END)==='
        
        return new_data



        


        
