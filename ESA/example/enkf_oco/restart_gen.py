import  ESA.util.bpch2_rw_smp as brw
from pylab import *
import ESA.util.gen_plots as gpl
import numpy as npy 
import ESA.util.time_module as tm
class geos_chem_restart_file_cl:
    """
    Class for managing restart files for GEOS-Chem run 
   
    Members:
    ==================================================
    
    1. ntracers:<int>: number of tracers
    2. org_flnm:<str>: file name for original restart file 
    3. bpch2:<bpch2_file_rw>: class for bpch2 file  
    

    Functions
    ==================================================
    1. __init__: initialization 
    2. create_enr_restart_file: create restart file for ensemble runs 
    
    """
    
    def __init__(self, org_flnm=None):
    
        self.ntracers=0
        
        if (org_flnm<>None):
            bpch2=brw.bpch2_file_rw(org_flnm, 'r', do_read=1)
            print bpch2.stat
            self.org_flnm=org_flnm
            self.bpch2=bpch2
        else:
            self.org_flnm=None
            self.bpch2=None
        
    
    def create_enr_restart_file(self, newfile, \
                                    ntracers, \
                                    newtau=None, \
                                    filling_value=1.0e-6, \
                                    sel_idx=1,\
                                    use_selected_tracer=False):
        
        """ make a new_restart file for ensemble run 

        Inputs:
        ------------------------------------
        
        1. newfile:<str>: name of the new restart file 
        2. ntracers:<integer>: number of tracers 
        3. newtau:<float>: time 
        4. filling_value:< float>: values used to fill the restart data
        5. sel_idx:<int>: the ID for tracer in the original restart file, 
        which field will be used as template. 
        6. use_selected_tracer:<T/F>: if true, all the initial 
        fields will be intialized as the selected one. 
        
        
        Returns:
        -------------------------------------
        1. ntracers:<int>: the number of tracers included in the restart file 
        
        """

        # S1: read in the selected tracer from original file 

        data_list, found=self.bpch2.get_data(None, sel_idx, None, None)
        if (found[0]==0):
            print 'No record found in original file'
            return None
        
        bpdata=data_list[0]
        
        # S2: open file to write 
        
        bpch2_w=brw.bpch2_file_rw(newfile, 'w')
        ifunit=95
        
        bpch2_w.open_bpch2_w(ifunit,'test')
        # #c: change time to the required one 
        if (newtau<>None):
            tau1, tau0=bpdata.get_attr(['tau0', 'tau1'])
            tau1=newtau+tau1-tau0
            print 'restart_gen> new tau0:', tau0
            print 'restart_gen> new tau1:', tau1
            
            tau0=newtau
            bpdata.set_attr('tau0', tau0)
            bpdata.set_attr('tau', tau1)
        
        sl=r'write restart data for %d tracers into file ' % ntracers
        print 
        print '='*20+'     restart_gen    '+'='*20
        print sl+newfile.strip()
        
        # S3: write data for ntracer 

        if (not use_selected_tracer):
            # #c: set values to the filling values
            bpdata.data=bpdata.data-bpdata.data+filling_value
            
        # #c: write the fields for ntracers

        for i in range(0, ntracers):
            bpdata.ntracer=i+1
            traname=bpdata.get_attr(['name'])
            print i+1, traname[0].strip(), bpdata.category.strip(), bpdata.ntracer, \
                bpdata.unit.strip(), shape(bpdata.data), \
                max(bpdata.data.flat), min(bpdata.data.flat)
            
            bpdata.write_to_bpch2_file(ifunit)
        # S4: close file
        
        bpch2_w.close_bpch2_file()
        print '-'*30+'Done'+'-'*30
        print '='*80
        
        return ntracers
    
        
        # write data back 
        

# <<< TEST >>> 

if (__name__=="__main__"):
    
    rsf=geos_chem_restart_file_cl('restart.ST001.EN0001-EN0002.2009010100')
    ntracers=28
    yyyy=2009
    mm=1
    dd=1
    doy=tm.day_of_year(yyyy=yyyy, mm=mm, dd=dd)
    tau0=tm.doy_to_tai85(doy, sec=0, yyyy=yyyy)
    tau0=tau0/3600.0
    
    newfile='restart.20090101'
    
    rsf.create_enr_restart_file(newfile, \
                                    ntracers,\
                                    newtau=tau0, \
                                    filling_value=1.0e-8, \
                                    sel_idx=1,\
                                    use_selected_tracer=False)
    
    

    
    
    
    
