""" class for reading and writing GEOS-Chem binary punch files
    history:
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.22
    History: v0.95, 2013.03.02
    

    
    Classes:
    =============================================
    1. diag_info_cl: class for diaginfo file
    2. tracer_info_cl: class for tracer info file 
    3. gcfile_desc_cl: class for accessing GEOS-Chem bpch2 files
    

    Functions:
    ================================================
    
    1.  is_bpch_matched: check whether it is the required tracer 
    2.  read_bpch_to_data_list: read 3D files from 
    bpch file to list oif ctm_field
    
    3.  read_bpch_to_profile_list: read 3D files from bpch file 
    and sample them along latitude and logitudes.  
    
    4.  open_bpch2_write: open a new bpch2 file 
    
    5.  close_bpch2_file: close one bpch2 file 
    
    6.  write_record_to_bpch2: write a data set into bpch2 file 
    
    7.  close_bpch2_file: close one bpch2 file 
    
"""

import ESA.util.bpch2_rw_py as brw

import ESA.util.time_module as tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.util.pres_m as pres_m

import ESA.util.message_m as msm 
import ESA.atmosphere.ctm_field_m as cfld
import ESA.atmosphere.ctm_profile_m as cprof
import ESA.atmosphere.gc_grid_3d as gg3d
import numpy as npy

import otool_gcfile_io as ogco
 

datapath='/home/lfeng/local_disk/otool_data/enkf_output/XYYYYX/'
flnm='ts_satellite.STXENRSTEPX.ENXEMSTX-ENXEMENDX.XYYYYXXMMXXDDX.bpch'
ftracerinfo='tracerinfo.STXENRSTEPX.ENXEMSTX-ENXEMENDX.dat'
fdiaginfo='diaginfo.STXENRSTEPX.ENXEMSTX-ENXEMENDX.dat'

yyyy, mm, dd=2009, 2, 12


fdesc=ogco.gcfile_desc_cl(flnm, ftracerinfo, \
                              fdiaginfo,\
                              yyyy, mm, dd,\
                              is_enr_file=True,\
                              categorys=None,\
                              tracers=None,\
                              taus=None,\
                              tranames=None)

fdesc.set_file_path(datapath)
fdesc.set_filename_format(flnm)


fdesc.set_tracerinfo_file_path(datapath)
fdesc.set_tracerinfo_filename_format(ftracerinfo)

fdesc.set_diaginfo_file_path(datapath)
fdesc.set_diaginfo_filename_format(fdiaginfo)

enr_yst_lst=[2009,2009]
enr_dst_lst=[1,1]
enr_em_st_lst=[1,81]
enr_em_end_lst=[81,146]
enr_step_lst=[0,0]


nset=fdesc.read_enr_file(yyyy, mm, dd, \
                             enr_yst_lst,\
                             enr_dst_lst,\
                             enr_step_lst,\
                             enr_em_st_lst,\
                             enr_em_end_lst)



print 'nset', nset 

traname_lst, category_lst, tid_lst, sid_lst, co2_lst=fdesc.get_data(traname='CO2')




traname_lst, category_lst, tid_lst, sid_lst, sp_lst=fdesc.get_data(traname='PSURF')

print traname_lst[0:2]
print tid_lst[0:2]

print len(sp_lst)

sp=sp_lst[0]
print npy.shape(sp)





    
