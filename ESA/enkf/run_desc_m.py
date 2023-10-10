"""
    Class for decoding ensemble run description 
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15

    Classes:
    =================================================
    1. enr_desc_cl:  Information about ensemble runs   

"""
import numpy as npy
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.otool_descfile_io as ofo
import ESA.util.time_module as tm
import ESA.util.otool_var_io as ovio
import run_desc_file_m as rdescio


# ====== < CLASSES >=======

class enr_desc_cl:
    
    """
    
    Information about ensemble runs    
    
    Members:
    ==========================================
    # current time 
    
    1.  yyyy:<int>: year
    2.  mm: <int>: month
    3.  dd: <int>: day 
    
    4.  flnm_ems_lst:<list, t:list>: list of emission file names
    5.  flnm_mod_lst:<list, t:list>: list of tracer file names

    # current coverage 

    6.  step_lst:<list, t:int>: steps
    7.  st_year_lst:<list, t:int>: starting year 
    8.  end_year_lst:<list, t:int>: ending year

    9.  st_doy_lst:<list, t:int>: starting day of year 
    10.  end_doy_lst:<list, t:int>: ending day of year
    
    11.  st_no_lst:<list, t:int>:  first ensemble member  
    12. end_no_lst:<list, t:int>: last ensemble member 
    13. nset: set of the number 
    14. attr_dict:<dict>: attributes
    
    15. enr_table_lst:<list>: list of ensemble run  description tables
    
    # file access
    16. fdesc: file descript for ensemble run description 
    17. fopen: open ensemble run description file 
    18. fread: read ensemble run description  file 
    19. fget:  get description table lis for current time 

    Functions:
    ===========================================================
    1. __init__:initialization
    2. load_enr_table: read description table 
    3. get_enr_table_list: check (or add in) description table for given time  (such as year)
    
    4. extract_coverage:  get ensemble runs covering a given time 
            
    5. copy: make a copy of itself 
    
    """
    
    def __init__(self, flnm,\
                     datapath, \
                     yyyy, mm, dd,\
                     fopen=rdescio.open_enr_desc_file,\
                     fread=rdescio.read_enr_desc_file,\
                     fclose=rdescio.close_enr_desc_file,\
                     fget=rdescio.get_enr_desc_table,\
                     fio_keywords={},\
                     **keywords): 
        """
        
        1.  yyyy:<int>: year
        2.  mm: <int>: month
        3.  dd: <int>: day 
        4.  enr_table_lst:<list, t:recarray>: configure table for ensemble runs 
        
        ---required columns for each members
        --->year_st:<int>: starting year
        --->doy_st:<int>: starting doy 
        --->year_end:<int>: end year
        --->doy_end:<int>: end doy
        --->mem_st:<int>: the first ensemble no 
        --->mem_end:<int>: the last ensemble no 
        
        """
        
        # S1: 
        self.flnm=flnm
        self.datapath=datapath
        
        # S2: members for current coverage
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd
        self.step_lst=None
        self.st_year_lst=None
        
        self.flnm_ems_lst=None
        self.flnm_mod_lst=None
        
        self.end_year_lst=None
        self.st_no_lst=None
        self.end_no_lst=None
        self.idx_mod_lst=None
        
        self.nset=0
        
        # S3: attribute
        self.attr_dict={}
        self.attr_dict.update(keywords)
        
        
        # S4: IO functions
        self.fopen=fopen
        self.fclose=fclose
        self.fread=fread
        self.fget=fget
        
        self.enr_table_lst=None
        self.fdesc=None
        self.fio_keywords=fio_keywords
        
        self.load_enr_table(**fio_keywords)
        
        if (self.fdesc<>None):
            self.get_enr_table(yyyy, mm, dd,**fio_keywords)
        
        print len(self.enr_table_lst)
        
        if (self.enr_table_lst<>None):
            self.extract_coverage(yyyy, mm, dd)
            
            
            
    def load_enr_table(self,**fio_keywords):
        
        """
        Read ensemble description files into fdesc
        
        Inputs:
        =================================
        1. fio_keywords:<dict>: keywords for fio
        
        """
        
        self.fdesc=self.fopen(self.flnm, \
                                  self.datapath, \
                                  self.yyyy, \
                                  self.mm, \
                                  self.dd,\
                                  **fio_keywords)
        
        self.fdesc=self.fread(self.fdesc, **fio_keywords)
        
        
    
    def get_enr_table(self, yyyy, mm, dd, **fio_keywords):
        """
        check ensemble run description file 
        
        Inputs:
        -----------------------------------------------
        1. yyyy, mm, dd:<int>: time 

        Returns:
        -----------------------------------------------
        1. enr_table:<recarray>: description table for ensemble runs 
        which covers given date. 
        
        
        """
        # S1: check whether the table list covering the latest date (i.e., year)
        
        
        enr_table=self.fget(self.fdesc,  yyyy, mm, dd, **fio_keywords)
        # S2: get the table list 
        print 'len again', len(self.fdesc.data_lst)
        self.enr_table_lst=self.fdesc.data_lst
        
        
        return enr_table
    
    
    def extract_current_time_step(self, yyyy, mm, dd):
        """ extract informations for runs covering given time 
        
        Inputs:
        ---------------------------------------------
        1. yyyy, mm, dd:<int>: time 
    
        Returns:
        ---------------------------------------
        1. st_yyyy, st_doy, end_yyyy, end_doy:<list, t>:starting year, doy and ending year, doy
        
        
        """
        enr_table=self.get_enr_table(yyyy, mm, dd)
        st_yyyy=enr_table[ 'year_st']
        st_doy=enr_table['doy_st']
        end_yyyy=enr_table[ 'year_end']
        end_doy=enr_table['doy_end']
        

        st_yyyy=st_yyyy.astype('int')
        st_doy=st_doy.astype('int')

        end_yyyy=end_yyyy.astype('int')
        end_doy=end_doy.astype('int')
        
        st_date=st_yyyy*1000+st_doy
        end_date=end_yyyy*1000+end_doy
        
        st_date=st_date.tolist()
        end_date=end_date.tolist()
        print 'st', st_date
        print 'end', end_date
        
        st_date=set(st_date)
        st_date=list(st_date)
        st_date=npy.array(st_date)
        st_date=npy.sort(st_date)
        
        
        
        # print 'end_date', end_date
        
        end_date=set(end_date)
        end_date=list(end_date)
        end_date=npy.array(end_date)
        end_date=npy.sort(end_date)

        
        
        st_yyyy=st_date/1000
        st_yyyy=st_yyyy.astype(int)
        st_doy=st_date-1000*st_yyyy
        
        end_yyyy=end_date/1000
        end_yyyy=end_yyyy.astype(int)
        end_doy=end_date-1000*end_yyyy
        

        return st_yyyy, st_doy, end_yyyy, end_doy
    
        
        
        
        
    def extract_coverage(self, yyyy, mm, dd, enr_table_lst=None):
        
        """ extract informations for runs covering given time 
        
        Inputs:
        ---------------------------------------------
        1. yyyy, mm, dd:<int>: time 
        2. enr_table_lst:<list, t:recarray>: configure table for ensemble runs 
        
        ---required columns for each members
        --->year_st:<int>: starting year
        --->doy_st:<int>: starting doy 
        --->year_end:<int>: end year
        --->doy_end:<int>: end doy
        --->mem_st:<int>: the first ensemble no 
        --->mem_end:<int>: the last ensemble no 
        
        Returns:
        ----------------------------------------------
        1. self.nsets:<int>: number of the 
        
        Notes:
        
        """
        # S1: load enr_table_lst
        if (enr_table_lst==None):
            # #c: forced to load new table  if necessary
            enr_table=self.get_enr_table(yyyy, mm, dd)
            enr_table_lst=self.enr_table_lst
        
        
        
        # S2: initialize list for current coverage
        
        self.step_lst=[]
        self.st_year_lst=[]
        self.end_year_lst=[]
        self.st_no_lst=[]
        self.end_no_lst=[]
        self.flnm_ems_lst=[]
        self.flnm_mod_lst=[]
        self.idx_mod_lst=[]
        
        self.st_doy_lst=[]
        self.end_doy_lst=[]

        self.mm=mm
        self.yyyy=yyyy
        self.dd=dd
        
        
        self.nset=0
        
        # S3: loop over table_lst
        
        for enr_table in enr_table_lst:
            
            st_yyyy=enr_table[ 'year_st']
            st_doy=enr_table['doy_st']
            end_yyyy=enr_table[ 'year_end']
            end_doy=enr_table['doy_end']
            st_mem=enr_table[ 'mem_st']
            end_mem=enr_table['mem_end']
            
            flnm_ems=enr_table['flnm_ems']
            flnm_mod=enr_table['flnm_mod']
            
            idx_mod=enr_table['index_method']
                        
            
            nentry=npy.size(st_yyyy)
            doy=tm.day_of_year(yyyy, mm, dd)
            # covert time to yyyyDOY
            
            cur_time=yyyy*1000+doy
            
            st_time=st_yyyy*1000+st_doy
            end_time=end_yyyy*1000+end_doy
            
            step=enr_table['step']
            usd_idx=npy.where((st_time<=cur_time) & (end_time>cur_time))
            usd_idx=npy.squeeze(usd_idx)
            self.nset=self.nset+npy.size(usd_idx)
            
            # S4: append records covering current day 
            
            for i in usd_idx:
                self.step_lst.append(int(step[i]))
                self.st_year_lst.append(int(st_yyyy[i]))
                self.end_year_lst.append(int(end_yyyy[i]))

                self.st_doy_lst.append(int(st_doy[i]))
                self.end_doy_lst.append(int(end_doy[i]))

                self.st_no_lst.append(int(st_mem[i]))
                self.end_no_lst.append(int(end_mem[i]))
                self.flnm_ems_lst.append(flnm_ems[i])
                self.flnm_mod_lst.append(flnm_mod[i])
                
                self.idx_mod_lst.append(int(idx_mod[i]))
                
                
            # loop i end
        
        # loop table end
        
        
        return self.nset


#====< Tests >======

if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
    dd=2
    mm=3
    yyyy=2009
    flnm='enr_cfg_XYYYYX.dat'

    
    
    run_desc=enr_desc_cl(flnm,\
                             datapath, \
                             yyyy, mm, dd,\
                             fopen=rdescio.open_enr_desc_file,\
                             fread=rdescio.read_enr_desc_file,\
                             fclose=rdescio.close_enr_desc_file,\
                             fget=rdescio.get_enr_desc_table)
    
    mm=6
    enr_table=run_desc.get_enr_table(yyyy, mm, dd)
    print enr_table
    
    
    nset=run_desc.extract_coverage(yyyy, mm, dd)
    print 'coverage' 
    
    print run_desc.st_year_lst
    print run_desc.st_doy_lst
    print run_desc.end_doy_lst
    
    
    st_yyyy, st_doy, end_yyyy, end_doy=run_desc.extract_current_time_step(yyyy, mm, dd)
    print st_yyyy
    print st_doy
    print end_yyyy
    print end_doy
    
    
    print nset
    
