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
import ESA.util.otool_var_io as ovio
import ESA.util.pres_m   as pm

import numpy as npy


#=======<<< FUNCTIONS >>>======



def is_bpch_matched(tcat, ttracer, tname, ttau0, ttau1,\
                        categorys=None,tracers=None, \
                        taus=None, tranames=None):

    """ 
    
    Check whether a bpch set meeting the selection criteria
    
    
    Inputs:
    -------------------------------------------------
    1. tcat:<str>: category 
    2. ttracer:<integer>: tracer id 
    3. tname:<str>: name of the tracer 
    4. ttau0:<float>: start tau
    5. ttau1:<float>: end tau 

    6. categorys:<list, t:str>: required  category 
    7. tracers:<list, t:integer>: required tracer numbers
    8. taus: <list, t:float>:  required tau 
    9. tranames:<list, t:float>: required names

    Returns:
    ==================================================
    1. matched: <True>: if the criteria are met. 
    
    Notes:
    =================================================
    1. tranames are the short names  
    2. ttracer are the ID stored in the bpch2 file, 
    ---instead of the ones after being shifted by category offset 
    
    3. fill values are treated as switching off the criteria. 
    
    
    """

    # duplicate inputs if given into array 
    
    
    array_cat=None
    array_tracer=None
    array_tau=None
    array_traname=None
    
    matched=True
    
    # S1 convert values to be  array 
    
    ## category 
    if (categorys<>None):
        array_cat=npy.array(categorys)
        nval_cat=npy.size(array_cat)
    else:
        nval_cat=0
        
    ## tracer ID 
        
    if (tracers<>None):
        
        array_tracer=npy.array(tracers)
        nval_tracer=npy.size(array_tracer)
    else:
        nval_tracer=0
    
    ## taus
        
    if (taus<>None):
        array_tau=npy.array(taus)
        nval_tau=npy.size(array_tau)
    else:
        nval_tau=0
    ## trace name (not the full name)
    
    if (tranames<>None):
        array_traname=npy.array(tranames)
        nval_traname=npy.size(array_traname)
        
    else:
        nval_traname=0
    
    nval=max(nval_cat, nval_tracer,  nval_tau,  nval_traname)
        
    if (nval==0):
        # no selection is required 
        return matched

    
    score_cat=npy.ones(nval)
    score_tracer=npy.ones(nval)
    score_tau=npy.ones(nval)
    score_traname=npy.ones(nval)
    
    if (array_cat<>None):
        # if category comparison is chosen
        chosen=(array_cat==tcat.strip()) | (array_cat==oob.fill_val_str)
        # print len(tcat)
        # print len(tcat.strip())
        # chosen=(array_cat==tcat) | (array_cat==oob.fill_val_str)
        #
        
        score_cat=npy.where(chosen,1,0)
        

    # 
    
    if (array_tracer<>None):
        # if tracer id comparison is chosen
        
        chosen=(array_tracer==ttracer) | (array_tracer==oob.fill_val_int)
        score_tracer=npy.where(chosen, 1,0)
    
    if (array_tau<>None):
        # if tau selection is given 
        score_tau=npy.where(((array_tau>=ttau0) & (array_tau<ttau1)), 1,0)
        chosen_idx=npy.where(array_tau==oob.fill_val)
        score_tau[chosen_idx]=1
    
				
    if (array_traname<>None):
        # if trace names are compared
        
        chosen=(array_traname==tname.strip()) | (array_traname==oob.fill_val_str)
        
        score_traname=npy.where(chosen, 1,0)
        
    
    final_score=score_cat*score_tracer*score_tau*score_traname
    # if all the criteria have been ticked 
    
    final_score=npy.sum(final_score)
    if (final_score==0):
        matched=False
    
    return matched





class tracer_info_cl:
    
    """ class for infomation of the tracers 
    included in the GEOS-Chem simulation
    
    Members:
    -----------------------------------------------------
    1. name:<list, t:str>: short name of tracers  
    2. fullname:<list, t:str>: full name of tracer  
    3. tracer:<list/array, t:integer>: id of the tracer  
    4. scale:<list, t:float>: scaling factor 
    5. unit: <unit, t:str>: unit of the tracers 
    6. c: <list, t:integer>: category id of the tracer. 
    7. molew:<list,t:float>: mass of tracer
    
    
    Functions:
    ----------------------------------------------
    
    1. __init__: initialization
    2. read_info: decode lines into different fields 
    3. get_tracer_info:check information of the tracer 
    4. load_tracer_info(self, flnm): Read tracer infor from file
        
        
    """
    

    def __init__(self, flnm):

        self.name=list()  # tracer name 
        self.fullname=list()   # tracer full name  
        self.tracer=list()     # tracer id 
        self.scale=list()     # scaling factor 
        self.unit=list()     # unit 
        self.c=list()       # category id 
        self.molew=list()  # molecular weight 
        
        file_found=False

        if (flnm.strip()<>""):
                
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                    
                msg='Tracer info file '+flnm+' not found'
                msm.show_err_msg(msg)
                
                file_found=False
        
        if (file_found):
            lines=fl.readlines()
            fl.close()
            self.read_info(lines)
        
    
    def read_info(self, lines):
        
        """decode lines in tracer info file into different fields
        Inputs:
        ----------------------------------------------
        1. lines:<list, str>: lines in tracerinfo file 
        --- a typical line include
        name, fullname, molecular mass, scaling factor, tracer_id, unit, category id 

        """

        
        for line in lines:
            if (line[0]<>'#'):  
                # if not a comment line 
                sname, sfullname, smw, sc, stra, sscal, sunit= \
                    line[0:8], line[9:39], line[39:49], \
                    line[49:53], line[53:62], line[62:72], line[73:] 
                
                self.name.append(sname.strip())
                self.fullname.append(sfullname.strip())
                self.molew.append(float(smw))
                self.c.append(float(sc))
                self.scale.append(float(sc))
                self.tracer.append(int(stra))
                self.unit.append(sunit.strip())
        
        # converted tracer ID list to an array 
        self.tracer=npy.array(self.tracer)
        
    

    def get_tracer_info(self, tracer_in):
        
        """ check information of the tracer 
        Inputs:
        ----------------------------------------
        1. tracer_in:<integer>: (real) tracer id 
        
        Returns:
        ------------------------------------
        1. name:<str>: name of the tracer
        2. fullname:<str>: full name of the tracer 
        3. molew: <float>: molecular mass
        4. scale: <float>: scaling factor 
        5. unit: <str>: unit 

        Notes:
        ----------------------------
        1. If tracer_in is not found in self.tracer, returns will be filled with "???" or with oob.fill_val
        
        
        
        """
        fill_str="???"
        
        if (len(self.tracer)==0):
            ## if no traceinfo file found 
            
            return fill_str, fill_str, \
                oob.fill_val, oob.fill_val, fill_str
        
        
        idx=npy.where(tracer_in==self.tracer)
        
        if (npy.size(idx)>=1):
            ##  if found 
            tid=idx[0]
            name=self.name[tid]
            fullname=self.fullname[tid]
            scale=self.scale[tid]
            unit=self.unit[tid]
            c=self.c[tid]
            molew=self.molew[tid]
            return name, fullname, molew, scale, unit 
        else:
            ## if not found 
            
            return fill_str, fill_str, oob.fill_val, \
                oob.fill_val, fill_str
        
        


        

    def load_tracer_info(self, flnm):
        """
        Read tracer infor from file
        
        Inputs:
        ---------------------------------------
        1. flnm:<str>: file name for tracer info 
        
        """
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                msg='No tracer info file: '+flnm.strip()
                msm.show_err_msg(msg)
                file_found=False
        
        if (file_found):
            lines=fl.readlines()
            fl.close()
            self.read_info(lines)

    



class diag_info_cl:
    """
    class for infomation of the tracers 
    included in the GEOS-Chem simulation
    
    Members:
    -----------------------------------------------------
    1. name:<list, t:str>: short name of tracers  
    2. fullname:<list, t:str>: full name of tracer  
    3. tracer:<list/array, t:integer>: id of the tracer  
    4. scale:<list, t:float>: scaling factor 
    5. unit: <unit, t:str>: unit of the tracers 
    6. c: <list, t:integer>: category id of the tracer. 
    7. molew:<list,t:float>: mass of tracer
    
    
    Functions:
    ----------------------------------------------
    1. __init__: initialization
    2. read_info: decode lines into different fields.
    3. get_offset:  get the offset for given category 
    4. load_diag_info: read diag info from file
    

    """
    def __init__(self, flnm):
        """
        initialization

        Inputs:
        ------------------------------------
        1. flnm :<str>:  name of the diaginfo file
        
        """
        self.category=list()
        self.comment=list()
        self.offset=list()
        
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                msg='Diag file '+flnm+' not found'
                file_found=False
        
        if (file_found):
            lines=fl.readlines()
            fl.close()
            self.read_info(lines)




    def read_info(self,lines):

        """ read information from lines in diaginfo file  
        
        Inputs:
        -------------------------------
        1. lines: <list, t:str>: list of lines 
        --typical line includes:  offset, category and comment 
        
        """
        
        for line in lines:
            # line=line.strip()
            if (line[0]<>'#'):
                ## if not a comment line 
                
                soffset, scategory, scomment= line[0:8], line[9:49], line[49:] 
                
                # print  sname, sfullname, smw, sc, stra, sscal, sunit
                self.offset.append(int(soffset))
                self.category.append(scategory.strip())
                self.comment.append(scomment.strip())
    

    def get_offset(self, category):
        
        """
        get the offset for given category 
        Inputs:
        -----------------------
        1. category:<str>: tracer category 
        

        Returns:
        1. offset:<integer>: offset of the category 
        
        """


        nid=len(self.category)
        found=0
        scate=category.strip()
        for id in range(nid):
            if (scate in self.category[id]):
                found=1
                offset=self.offset[id]
                
        if (found):
            return offset
        else:
            return oob.fill_val_int
    
    def load_diag_info(self, flnm):
        
        """ read diag info from file
        
        Inputs:
        -------------------------------
        1. flnm:<str>: file name 

        
        """
        
        file_found=False
        
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                msg='no diag info file: ',  flnm.strip()
                msm.show_err_msg(msg)
                file_found=False
        if (file_found):
            lines=fl.readlines()
            fl.close()
            self.read_info(lines)



            
def read_bpch_to_data_list(flnm, ftracerinfo="", \
                               fdiaginfo="",\
                               categorys=None,tracers=None, \
                               taus=None, tranames=None\
                               ):
    
    
    """ read in the bpch data which meet the selection criteria 
    
    Inputs:
    -------------------------------------------------------------
    1. flnm:<str>: name of the bpch2 file 
    2. ftracerinfo:<str>: name of the tracerinfo file 
    3. fdiaginfo:<str>:name of the diag info file 
    4. categorys:<list, t:str>: required  category 
    5. tracers:<list, t:integer>: required tracer numbers
    6. taus: <list, t:float>:  required tau 
    7. tranames:<list, t:float>: required names


    Returns:
    -----------------------------------------------
    1. data_lst:<list, t:ctm_field_cl>: list of class object of ctm_field_cl.


    Notes:
    -----------------------------------------
    1. tracers is the ID stored in the bpch2 file,  instead of the 'real'
    --- tracer ID  after being shifted by category offsets read from diaginfo file 
    2. Detailed information such as traname, trafullname, tramolew, trascale and traunit 
    ---are read from tracerinfo file according to the real tracer ID 
    
    3. fill value means  the switching off the criteria. 
    
    
    """

    
    data_lst=[]


    funit=199
    
    # S1  open bpch2 file 
    
    
    print '--otool_gc_file: read_bpch_to_data_list'
    print '---bpch2 file (1) :', flnm.strip()
    
    
    fti,title,stat = brw.open_bpch2_for_read(funit,flnm)

    if (stat<>0):
        msg=r'error in read_bpch_to_data_list: %4.4d',  stat
        msm.show_err_msg(msg)
        return data_lst

    # S2 open tracer and diag info 
    
    print '---diaginfo file (2) :', fdiaginfo.strip()
    diaginfo=diag_info_cl(fdiaginfo)
    
    print '---tracerinfo file (3)  :', fdiaginfo.strip()
    tracerinfo=tracer_info_cl(ftracerinfo)
    
    
    # S3 loop over the data  
    
    vtracer_id,vhalfpolar,vcentre180,vni,vnj,vnl,\
        vifirst,vjfirst,vlfirst,vlonres,vlatres,\
        vtau0,vtau1,vmodelname,vcategory,vunit,vreserved,\
        vdata_array,stat = brw.read_bpch2_record(funit)
    
    while (stat==0):
            
        ## T1  get offset for the category 
        
        
        vcategory=vcategory.strip()
        
        offset=diaginfo.get_offset(vcategory)
            
        if (offset==oob.fill_val_int):
            ###c if no offset is found 
            offset=0
                
        
        real_id=vtracer_id+offset
        
        ## T2 information on the tracer
        
        traname, trafullname, tramolew, trascale, traunit=\
            tracerinfo.get_tracer_info(real_id)
        
        
        traname=traname.strip()
        trafullname=trafullname.strip()
        traunit=traunit.strip()
        
       
        
        
        vdata=npy.zeros([vni, vnj, vnl], float)                                                
        
        ## T3 check whether it meets the selection 
        vcategory=vcategory.strip()
        
        matched=is_bpch_matched(vcategory, vtracer_id,  \
                                    traname, vtau0, vtau1, \
                                    categorys=categorys,tracers=tracers,\
                                    taus=taus, tranames=tranames)
        
        
        if (matched):
            ## T4 setup grid 
            
            
            dgrid= gg3d.gc_grid_cl(vni, vnj, vnl,\
                                       ifirst=vifirst, \
                                       jfirst=vjfirst, \
                                       lfirst=vlfirst,\
                                       halfpolar=vhalfpolar, \
                                       centre180=vcentre180, \
                                       lonres=vlonres, \
                                       latres=vlatres)
            
            
            
            ## T5 resize data set to from ifirst, jfirst, and lfirst
            
            vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[(vifirst-1):(vifirst+vni-1), \
                                                       (vjfirst-1):(vjfirst+vnj-1), \
                                                       (vlfirst-1):(vlfirst+vnl-1)]
            
            ## T6 construct the class to  contain the data 
                        
            bpdata=cfld.ctm_field_cl(traname,\
                                         vdata, \
                                         unit=vunit.strip(), \
                                         ctm_grd=dgrid,\
                                         is_xyz_data=True,\
                                         extra_axis=None,\
                                         tracer_id=vtracer_id,\
                                         category=vcategory.strip(),\
                                         traname=traname,\
                                         tau0=vtau0,\
                                         tau1=vtau1, \
                                         modelname=vmodelname.strip(),\
                                         reserved=vreserved, \
                                         offset=offset, \
                                         fullname=trafullname, \
                                         molew=tramolew, \
                                         scale=trascale, \
                                         traunit=traunit, \
                                         full_id=real_id)
            
                              
            ##  added the data to the list 
            
            data_lst.append(bpdata)
            
        ###C if matched 
            
        ## read next record
        
        vtracer_id,vhalfpolar,vcentre180,vni,vnj,vnl,\
            vifirst,vjfirst,vlfirst,vlonres,vlatres,\
            vtau0,vtau1,vmodelname,vcategory,vunit,vreserved,\
            vdata_array,stat = brw.read_bpch2_record(funit)
    
    # S4 close file 
    
    stat=brw.close_bpch2_file(funit)
    print '--otool_gc_file: read_bpch_to_data_list exit'
    
    return data_lst
    




def read_bpch_to_profile_list(olon, \
                                  olat,\
                                  flnm,     \
                                  ftracerinfo="", \
                                  fdiaginfo="",\
                                  categorys=None,\
                                  tracers=None, \
                                  taus=None, \
                                  tranames=None\
                                  ):
    
    
    """ read in the bpch data which meet the selection criteria 
    Inputs:
    -------------------------------------------------------------
    1. olon:<array>: longitude
    2. olat:<array>: latitude 
    
    3. flnm:<str>: name of the bpch2 file 
    4. ftracerinfo:<str>: name of the tracerinfo file 
    5. fdiaginfo:<str>:name of the diag info file 
    6. categorys:<list, t:str>: required  category 
    7. tracers:<list, t:integer>: required tracer numbers (Note 1)

    8. taus: <list, t:float>:  required tau 
    9. tranames:<list, t:str>: required names
    

    Returns:
    -----------------------------------------------
    1. data_lst:<list, t:ctm_profile_cl>: list of class object of ctm_field_cl.
    
   
    Notes:
    -----------------------------------------
    1. tracers is the ID stored in the bpch2 file,  instead of the 'real'
    --- tracer ID  after being shifted by category offsets read from diaginfo file 
    2. Detailed information such as traname, trafullname, tramolew, trascale and traunit 
    ---are read from tracerinfo file according to the real tracer ID 
    
    3. fill value means  the switching off the criteria. 
    

    """

    
    data_lst=[]


    funit=199
    
    # S1  open bpch2 file 
    
    print '--otool_gc_file: read_bpch_to_profile_list'
    print '---bpch2 file (1) :', flnm.strip()

    fti,title,stat = brw.open_bpch2_for_read(funit,flnm)

    if (stat<>0):
        msg=r'error in read_bpch_to_profile_list: %4.4d',  stat
        msm.show_err_msg(msg)
        return data_lst

    # S2 open tracer and diag info 
    
    print '---diaginfo file (2) :', fdiaginfo.strip()
    diaginfo=diag_info_cl(fdiaginfo)
    
    print '---tracerinfo  file (3) :', ftracerinfo.strip()
    tracerinfo=tracer_info_cl(ftracerinfo)
    
    
    # S3 loop over the data  
    
    vtracer_id,vhalfpolar,vcentre180,vni,vnj,vnl,\
        vifirst,vjfirst,vlfirst,vlonres,vlatres,\
        vtau0,vtau1,vmodelname,vcategory,vunit,vreserved,\
        vdata_array,stat = brw.read_bpch2_record(funit)
    
    while (stat==0):
            
        ## T1  offset for the category 
        
        vcategory=vcategory.strip()
        
        offset=diaginfo.get_offset(vcategory)
            
        if (offset==oob.fill_val_int):
            # if offset not found 
            
            offset=0
                
        
        real_id=vtracer_id+offset
        
        ## T2 information on the tracer
        
        traname, trafullname, tramolew, trascale, traunit=\
            tracerinfo.get_tracer_info(real_id)
        
        
        traname=traname.strip()
        trafullname=trafullname.strip()
        traunit=traunit.strip()
        
       
        
        
        vdata=npy.zeros([vni, vnj, vnl], float)                                                
        
        ## T3 check whether it meets the selection criteria 
        
        
        matched=is_bpch_matched(vcategory, vtracer_id,  \
                                    traname, vtau0, vtau1, \
                                    categorys=categorys,\
                                    tracers=tracers,\
                                    taus=taus, tranames=tranames)
        
        if (matched):
            ## T4 setup grid 
            dgrid= gg3d.gc_grid_cl(vni, vnj, vnl,\
                                       ifirst=vifirst, \
                                       jfirst=vjfirst, \
                                       lfirst=vlfirst,\
                                       halfpolar=vhalfpolar, \
                                       centre180=vcentre180, \
                                       lonres=vlonres, \
                                       latres=vlatres)
            
            ## T5 resize data set to start from ifirst, jfirst, lfirst
            
            
            hintpl=dgrid.get_hintpl(olon, olat)
            vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[(vifirst-1):(vifirst+vni-1), \
                                                       (vjfirst-1):(vjfirst+vnj-1), \
                                                       (vlfirst-1):(vlfirst+vnl-1)]
           
            profiles=hintpl.get_profile(vdata)
            
            
            
            
            ## T6 construct the class to  contain the data 
                        
            bpprof=cprof.ctm_profile_cl(traname,\
                                            profiles, \
                                            unit=vunit.strip(), \
                                            lon=olon,\
                                            lat=olat,\
                                            tau=vtau0,\
                                            mask_val=oob.fill_val,\
                                            is_xz_data=True,\
                                            extra_axis=None,\
                                            tracer_id=vtracer_id,\
                                            category=vcategory.strip(),\
                                            traname=traname,\
                                            tau0=vtau0,\
                                            tau1=vtau1, \
                                            modelname=vmodelname.strip(),\
                                            reserved=vreserved, \
                                            offset=offset, \
                                            fullname=trafullname, \
                                            molew=tramolew, \
                                            scale=trascale, \
                                            traunit=traunit, \
                                            full_id=real_id)
            
            
                              
            ##  added the data to the list 
            
            data_lst.append(bpprof)
            
            ## return the data 
            
        vtracer_id,vhalfpolar,vcentre180,vni,vnj,vnl,\
            vifirst,vjfirst,vlfirst,vlonres,vlatres,\
            vtau0,vtau1,vmodelname,vcategory,vunit,vreserved,\
            vdata_array,stat = brw.read_bpch2_record(funit)
    
    stat=brw.close_bpch2_file(funit)
    print '--otool_gc_file: read_bpch_to_profile_list: exit'

    return data_lst






def open_bpch2_write(funit, flnm, title=""):

    """ open one bpch2 file to write 
    Inputs:
    -------------------------
    1. funit:<int>: file unit handle 
    2. flnm:<str>: file name 
    3. title:<str>: title of the data set to be written. 
    

    Returns:
    ---------------------
    1. stat:<int>: status of the operation


    """

    stat=brw.open_bpch2_for_write(funit,flnm, title)
    return stat

def close_bpch2_file(funit):

    """ open one bpch2 file to write 
    Inputs:
    -------------------------
    1. funit:<int>: file unit handle 
    
    Returns:
    ---------------------
    1. stat:<int>: status of the operation

    """

    if (funit<>None):
        stat=brw.close_bpch2_file(funit)
    return stat

def write_record_to_bpch2(funit, ctm_data,  **keywords):
    
    """ write data into one open bpch2 file 
    
    Inputs:
    -----------------------------------------
    
    1. funit:	 integer- the unit number of the bpch2 file 
    2. ctm_data:<ctm_field_cl/array>: gridded data with the ctm_grd 
    3. keywords: Additional inputs to replace what stored in data 
    --- reserved keywords 
    --->'ctm_grd':<gc_grid_cl>: grid 
    --->'unit':<str>: unit        
    --->'modelname': <str>: model name 
    --->'tau0':<float>: starting tau for the data set 
    --->'tau1':<float>: end tau for the data set
    ---'>tracer_id:':<int>: tracer number 
    --->'category': <str>: tracer category 
    --->'reserved': <str>: reserved fields.
    
    
    Returns:
    ---------------------
    1. stat:<int>: status of the operation

    """
    

    # S1 fill the information 
    
    if ('ctm_grd' in keywords):
        grid=keywords['grid']
    else:
        grid=ctm_data.get_grid()
    
    
    ## T1 grid set up 

    ix=grid.ix
    jx=grid.jx	
    lx=grid.lx
    
    lonres=grid.lonres
    latres=grid.latres
    ifirst=grid.ifirst
    jfirst=grid.jfirst
    lfirst=grid.lfirst
    
    halfpolar=grid.halfpolar
    centre180=grid.centre180
    
    ## T2 tracer identity  
    
    if ('tracer_id' in keywords):
        ntracer=keywords['tracer_id']
    else:
        ntracer=ctm_data.get_attr('tracer_id')
    
    
        
    if ('category' in keywords):
        category=keywords['category']
    else:
        category=ctm_data.get_attr('category')
    
    if ('unit' in keywords):
        unit=keywords['unit']
    else:
        unit=ctm_data.unit
    
    
    ## T3 time 
        
    if ('tau0' in keywords):
        tau0=keywords['tau0']
    else:
        try: 
            tau0=ctm_data.get_attr('tau0')
        except:
            tau0=oob.fill_val

    
    
    if ('tau1' in keywords):
        tau1=keywords['tau1']
    else:
        try:
            tau1=ctm_data.get_attr('tau1')
        except:
            tau1=oob.fill_val
    
            

    ## T4 modelname  
    
    
    if ('modelname' in keywords):
        modelname=keywords['modelname']
    else:

        # modelname may not been necessary 
        try:
            modelname=ctm_data.get_attr('modelname')
        except:
            modelname=""
    
    
    ## reserve 
    
    if ('reserved' in keywords):
        
        reserved=keywords['reserved']
    else:

        # modelname may not been necessary 
        try:
            reserved=ctm_data.get_attr('reserved')
        except:
            reserved=""
    
    
    reserved=reserved.strip()
    if (len(reserved)==0):
        # ##T: reserved can not be empty 
        
        reserved="000000"
    

    # S2  write data to file         
    
    data=ctm_data.get_data()
    
    stat = brw.write_bpch2_data(funit,modelname,\
                                    category,reserved, \
                                    lonres,latres,\
                                    halfpolar,centre180,\
                                    ntracer,unit,tau0,tau1,\
                                    ifirst,jfirst,lfirst,data)
    
    
    return stat



#==========<<< CLASSES >>>==============



class gcfile_desc_cl:
    """
    Members:
    ====================================
    1. flnm:<str>: file name
    2, ftracerinfo:<str>: name for tracerinfo file
    3. fdiaginfo:<str>: name for diaginfo file
    4. categorys:<list, t:str>: list of categorys to be read
    5. tracers:<list,t:int>: list of tracers to be read 
    6. taus:<list, t:float>: list of taus to be read 
    7. tranames:<list, t:str>:list of tracer name
    8. fio_keywords: <dict>: extra inputs for file IO
    9. attr_dict:<dict>: attributes
    
    10. mask_val:<float>: fillings for bad or missing values
    11. data_lst:<list>: list of list of model data or profiles read from bpch2 file 
    12. ref_lst:<list, t:int>:list of the references.
    13. flnm_lst:<list, t:str>: list of the file names.
    14. tracerinfo_flnm_lst:<list, t:str>: list of the tracerinfo file names.
    15. diaginfo_flnm_lst:<list, t:str>: list of the diaginfo file names.
    
    16. cur_flnm: <str>: current filename 
    17. cur_tracerinfo_flnm: <str>: current tracer info file 
    18. cur_diaginfo_flnm:<str>: current diag info file 
    19. cur_ref:<current ref>: current reference number 
    
    

    
    Functions:
    =================================================
    1.  __init__: initialization 
    2. append_data: add the data set to lists of the class
    3. del_data: Remove data from the list 
        
    4. set_file_path: set path for data file
    5. set_filename_format: set template for filename

    6. set_tracerinfo_file_path: set path for tracerinfo  file
    7. set_tracerinfo_filename_format: set template for tracerinfo filename
    
    8. set_diaginfo_file_path:set path for diaginfo  file
    9. set_diaginfo_filename_format:set template for diaginfo filename
    10. construct_filename: construction file name using keywords 

    11. construct_diaginfo_filename: construction name for diag info file using keywords
    12. set_attr: assign one attribute
    13. get_attr: get one attribute
    14. read_file: Read file from GEOS-Chem runs
    15. read_file_to_profile_list:Read outputs from GEOS-Chem  runs to list of profiles (ctm_profile_cl)

    # ensemble runs 
    
    16. construct_enr_filename_dict:  Form dictionary (extensions) to construct name for ensemble runs
    17. read_enr_file: Read file from ensemble runs to a list of fields (gc_field_3d)
    18. read_enr_file_to_profile_list:Read outputs from ensemble run to list of profiles (ctm_profile_cl)
    
    19. get_data: find the data meeting criteria on traname, category and tracer id etc 
    
    # pressure grid

    20. compute_mod_pres: compute from read surface pressure 
       or directly fetch read model pressure
    
    
    
    
    """
    
    
    def __init__(self, flnm, ftracerinfo, \
                     fdiaginfo,\
                     yyyy, mm, dd,\
                     is_enr_file=False,\
                     categorys=None,\
                     tracers=None,\
                     taus=None,\
                     tranames=None,\
                     mask_val=oob.fill_val,\
                     fio_keywords={}, **keywords):
        
        """
        
        initialization 

        Inputs:
        --------------------------------------------------------------
        1. flnm:<str>: file name
        2, ftracerinfo:<str>: name for tracerinfo file
        3. fdiaginfo:<str>: name for diaginfo file
        4. categorys:<list, t:str>: list of categorys to be read
        5. tracers:<list,t:int>: list of tracers to be read 
        6. taus:<list, t:float>: list of taus to be read 
        7. tranames:<list, t:str>:list of tracer name
        8. mask_val:<float>: values for bad or missing data
        9. fio_keywors: <dict>: extra inputs for file IO
        10. keywords:<dict>: attributes
        
        """
        self.flnm=flnm
        self.yyyy=yyyy
        self.mm=mm
        self.dd=dd

        self.ftracerinfo=ftracerinfo
        self.fdiaginfo=fdiaginfo

        # storage 
        
        self.categorys=categorys
        self.tracers=tracers
        self.taus=taus
        self.tranames=tranames
        self.is_enr_file=is_enr_file
        
        
        self.fio_keywords=fio_keywords
        self.attr_dict={}
        
        if ('ot_type' in self.attr_dict):
            pass
        else:
            self.attr_dict.update({'ot_type':oob.ot_fdesc})

        self.attr_dict.update(keywords)
        self.data_lst=[]
        
        
        self.mask_val=mask_val

        self.cur_flnm=flnm
        self.cur_ref=0
        self.flnm_lst=[]
        
        self.tracerinfo_flnm_lst=[]
        self.diaginfo_flnm_lst=[]
        
        self.ref_lst=[]
        self.nset=0
    

    

    def append_data(self, flnm, tracerinfo_flnm, \
                        diaginfo_flnm, data_lst, ref):
        
        """
        
        add the data set to lists of the class
        
        Inputs:
        -------------------------------------------
        1. flnm: <str>:   file name 
        2. tracerinfo_flnm: <str>:   name for tracerinfo file 
        3. diaginfo_flnm: <str>:  name for diaginfo file  
        
        4. data: <array>:  raw data read from the file 
        5. ref:<str/numeric>: id (such as time etc) for the data set 
        
        """
        if (ref in self.ref_lst):
            for idx in self.ref_lst:
                if (ref==self.ref_lst[idx]):
                    self.flnm_lst[idx]=flnm
                    self.tracerinfo_flnm_lst[idx]=tracerinfo_flnm
                    self.diaginfo_flnm_lst[idx]=diaginfo_flnm
                    
                    self.data_lst[idx]=data
                    self.ref_lst[idx]=ref
                    break
        
        else:
            self.flnm_lst.append(flnm)
            self.tracerinfo_flnm_lst.append(tracerinfo_flnm)
            self.diaginfo_flnm.append(diaginfo_flnm)
            
        
            self.data_lst.append(data_lst)
            self.ref_lst.append(ref)
            self.nset=self.nset+1

    def del_data(self, idx, by_ref=False):
        
        """
        Remove data from the list 
        
        Inputs:
        -----------------------------------------------------
        
        1.idx :<int/obj>: index of the data set to be deleted 
        2. by_ref:<T/F>: if ture , index will be treated as ref 
        ---values for looking up ref_lst
        
        Returns:
        ---------------------------------------------------
        self.nset:<int>: the length of the current data lst 
        
        """
        
        if (by_ref):
            
            idx=self.find_index(idx)
            if (idx==self.mask_val):
                msg='Not found'
                msm.show_err_msg(msg)
                return self.nset
        
        
        del self.flnm_lst[idx]
        del self.tracerinfo_flnm_lst[idx]
        del self.diaginfo_flnm_lst[idx]
        
        del self.data_lst[idx]
        del self.ref_lst[idx]
        self.nset=self.nset-1
        
        return self.nset

    def set_file_path(self, path):

        """ set path for data file 
        Inputs:
        -----------------------------------
        1. path:<str>: path  for file 
        
        """
        self.set_attr('filepath', path)
    
    def set_filename_format(self, flnm_format):
        """set template for filename 
        Inputs:
        -----------------------------------
        1. flnm_format:<str>: format for file name 
        
        
        """

        self.set_attr('flnmformat', flnm_format)



    def set_tracerinfo_file_path(self, path):

        """ set path for tracerinfo  file 
        Inputs:
        -----------------------------------
        1. path:<str>: path  for tracerinfo file 
        """
        
        self.set_attr('tracerinfo_filepath', path)
    
    def set_tracerinfo_filename_format(self, flnm_format):
   
        """set template for tracerinfo filename 
        Inputs:
        -----------------------------------
        1. flnm_format:<str>: format for tracerinfo file name 
        
        
        """
        
        self.set_attr('tracerinfo_flnmformat', flnm_format)
    



    def set_diaginfo_file_path(self, path):

        """ set path for diaginfo  file 
        Inputs:
        -----------------------------------
        1. path:<str>: path  for tracerinfo file 
        """
        
        self.set_attr('diaginfo_filepath', path)
    
    def set_diaginfo_filename_format(self, flnm_format):
   
        """set template for diaginfo filename 
        Inputs:
        -----------------------------------
        1. flnm_format:<str>: format for diaginfo file name 
        
        
        """
        
        self.set_attr('diaginfo_flnmformat', flnm_format)
    


    def construct_filename(self, **keywords):
    
        """ construction file name using keywords 

        Inputs:
        ----------------------------
        
        1. keywords:<dict>: dictionary for file format and path
        ---expected Keynames:
        --->1. flnmformat:<str>: file format
        --->2. filepath:<str>: file path
        --->3. words used to replace dummies in file format.
        

        Returns
        ---------------------------------------
        1. sflnm:<str>: file name 
        
        
        """
        
        sflnm=self.get_attr('flnmformat')
        
        if ('filepath' in self.attr_dict):
            filepath=self.get_attr('filepath')
        else:
            filepath=""
        
        
        # construct the name 
        
        sflnm=ovio.construct_filename(filepath, sflnm, **keywords)
        
        return sflnm
    
    def construct_tracerinfo_filename(self, **keywords):
    
        
        """ construction name for tracer info file  use the keyword 

        Inputs:
        ----------------------------
        1. keywords:<dict>: dictionary for file format and path
        ---expected Keynames:
        --->1. tracerinfo_flnmformat:<str>: file format
        --->2. tracerfo_filepath:<str>: file path
        --->3. words used to replace dummies in file format.
        
        
        Returns
        ---------------------------------------
        1. sflnm:<str>: file name 
        

        """
        
        sflnm=self.get_attr('tracerinfo_flnmformat')
        if ('tracerinfo_filepath' in self.attr_dict):
            filepath=self.get_attr('tracerinfo_filepath')
        else:
            filepath=""
        
        # construct the name 
        
        sflnm=ovio.construct_filename(filepath, sflnm, **keywords)
        
        return sflnm
    
    def construct_diaginfo_filename(self, **keywords):
        
        """ construction name for diag info file using keywords
 
        Inputs:
        ----------------------------
        1. keywords:<dict>: dictionary for file format and path
        ---expected Keynames:
        --->1. diaginfo_flnmformat:<str>: file format
        --->2. diaginfo_filepath:<str>: file path
        --->3. words used to replace dummies in file format.
        

        Returns
        ---------------------------------------
        1. sflnm:<str>: file name 
        
        
        """
        
        sflnm=self.get_attr('diaginfo_flnmformat')
        if ('diaginfo_filepath' in self.attr_dict):
            filepath=self.get_attr('diaginfo_filepath')
        else:
            filepath=""
        
        # construct the name 
        
        sflnm=ovio.construct_filename(filepath, sflnm, **keywords)
        
        return sflnm

    
                    
    def set_attr(self, attr_name, value):
        """  assign one attribute
        Inputs:
        -------------------------------------------------------
        1. attr_name: <str>: attribute name
        2. value: <obj>:attribute value
        """
        
        self.attr_dict.update({attr_name:value})
        
        
    def get_attr(self, attr_name):
        """  get one attribute
        
        Inputs:
        ------------------------------------------------
            1. name:<str>: attribute name
        Returns:
        -----------------------------------------------
            1. val: <obj>: attribute value
         """
        
        if (attr_name in self.attr_dict):
            val=self.attr_dict[attr_name]
        else:
            val=oob.fill_val
        
        return val

    def construct_enr_filename_dict(self, \
                                        yyyy, dd, mm,\
                                        enr_yyyy, \
                                        enr_doy, \
                                        enr_step, \
                                        enr_em_st, \
                                        enr_em_end):
        
                                        
        """
        Form dictionary (extensions) to construct name for ensemble runs

        Input:
        ============================================
        1. enr_yyyy:<int>: year for ensemble run start 
        2. enr_doy:<int>:  doy for ensemble run start
        3. enr_step:<int>: step for ensemble run 
        4. enr_em_st:<int>: firt one in the ensemble 
        5. enr_em_end:<int>: last one in the ensemble 
        6. yyyy, dd, mm:<int>: output time 
        
        Notes:
            1. the default extension will be like: 
                STXENRSTEPX.ENXEMSTX-XEMEND.XYYYYXXMMXXDDX
                for example:
                 ST001.EN0001-EN0048.2009.01.01
                 stands for 
                 step=1
                 ensemble member of 0001-0049
                 yyyyy=2009, mm=1, dd=1
                 
                 
        """

        # S1: convert to string 
        
        em_yyyy=r'%4.4d' % enr_yyyy
        em_step=r'%3.3d'% enr_step 
        em_doy=r'%3.3d'% enr_doy
        
        em_st=r'%4.4d' % enr_em_st
        em_end=r'%4.4d' % enr_em_end
        
        yyyy=r'%4.4d' % yyyy
        mm=r'%2.2d' % mm
        dd=r'%2.2d' % dd
        
        
        # S2: build dictionary 
        
        exts_dict={'XEMSTX':em_st, \
                       'XEMENDX':em_end, \
                       'XENRYYYYX':em_yyyy, \
                       'XENRDOYX':em_doy, \
                       'XENRSTEPX':em_step, \
                       'XYYYYX':yyyy, \
                       'XMMX':mm, \
                       'XDDX':dd}
        
        
        return exts_dict
    
    
        
    def read_file(self, yyyy, mm, dd, \
                      ext="",\
                      **keywords):
        
        """
        Read outputs from GEOS-Chem runs to a list of 3D fields 
        
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. ent:<str>: string for file extension 
        3. keywords:<dict>: extra inputs:
        
        Returns:
        ==========================================================
        1. bp_lst:<list, t:gc_field_3d>: list of  data sets. 
        
        
        """
        data_lst=[]
        flnm_lst=[]
        tracerinfo_flnm_lst=[]
        diaginfo_flnm_lst=[]
        ref_lst=[]
        nset=1
        

        # S1: get file name 
        
        doy=tm.day_of_year(yyyy, mm, dd)
        
        exts_dict={'XYYYYX':yyyy, 'XMMX':mm, 'XDDX':dd, 'XDOYX':doy, 'XEXTX':ext}
        
        full_name=self.construct_filename(**exts_dict)
        fdiag=self.construct_diaginfo_filename(**exts_dict)
        ftracer=self.construct_tracerinfo_filename(**exts_dict)
        
        # S2: set storage list 

        ref_lst.append(iset)
        flnm_lst.append(full_name)
        tracerinfo_flnm_lst.append(ftracer)
        diaginfo_flnm_lst.append(fdiag)
        
        # S3: read in data
        bp_lst=read_bpch_to_data_list(full_name, ftracerinfo=ftracer, \
                                          fdiaginfo=fdiag,\
                                          categorys=self.categorys,tracers=self.tracers, \
                                          taus=self.taus, tranames=self.tranames)
        
            
            
        data_lst.append(bp_lst)
        
        # S4: update class member 
        
        self.data_lst=data_lst
        self.flnm_lst=ref_lst
        self.tracerinfo_flnm_lst=tracerinfo_flnm_lst
        self.diaginfo_flnm_lst=diaginfo_flnm_lst
        self.ref_lst=ref_lst
        self.nset=len(ref_lst)
        
        return bp_lst
    

    
    def read_file_to_profile_list(self,yyyy, mm, dd, \
                                      olon, olat,\
                                      ext="",\
                                      **keywords):
        
        """
        Read outputs from forecast runs to a list of ctm_profile_cl
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. olon:<array>: longitude
        3. olat:<array>:latitude 
        4. ext:<str>: file extension 
        
        5. keywords:<dict>: extra inputs:
        
        Returns:
        ==========================================================
        1. bp_lst:<list, t:prof_lst>: list of found profiles. 
        
        
        """
        
        # S1: get file name 
        
        doy=tm.day_of_year(yyyy, mm, dd)
        
        exts_dict={'XYYYYX':yyyy, 'XMMX':mm, 'XDDX':dd, 'XDOYX':doy, 'XEXTX':ext}
        
        full_name=self.construct_filename(**exts_dict)
        fdiag=self.construct_diaginfo_filename(**exts_dict)
        ftracer=self.construct_tracerinfo_filename(**exts_dict)
        
        # S2: storage for outputs 
        
        iset=0

        ref_lst=[]
        flnm_lst=[]
        tracerinfo_flnm_lst=[]
        diaginfo_flnm_lst=[]
        data_lst=[]
        
        ref_lst.append(iset)
        flnm_lst.append(full_name)
        tracerinfo_flnm_lst.append(ftracer)
        diaginfo_flnm_lst.append(fdiag)
        
        # S3: read in profiles 
        
        bp_lst=read_bpch_to_profile_list(olon, \
                                             olat,\
                                             full_name, \
                                             ftracerinfo=ftracer, \
                                             fdiaginfo=fdiag,\
                                             categorys=self.categorys,tracers=self.tracers, \
                                             taus=self.taus, tranames=self.tranames)
        
        
            
        data_lst.append(bp_lst)
        
        # ##L: end loop on set 
        
        self.data_lst=data_lst
        self.flnm_lst=ref_lst
        self.tracerinfo_flnm_lst=tracerinfo_flnm_lst
        self.diaginfo_flnm_lst=diaginfo_flnm_lst
        self.ref_lst=ref_lst
        self.nset=len(ref_lst)
        
        return bp_lst
    

    
    
    

        
    def read_enr_file(self, yyyy, mm, dd, \
                          enr_yst_lst,\
                          enr_dst_lst,\
                          enr_step_lst,\
                          enr_em_st_lst,\
                          enr_em_end_lst,\
                          **keywords):
        
        """
        Read outputs from ensemble runs to list of ctm_data_cl
        
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. enr_yyyy_lst:<int>: year for ensemble run start.
        3. enr_doy_lst:<int>:  doy for ensemble run start.
        4. enr_step_lst:<int>: step for ensemble run.
        5. enr_em_st_st:<int>: firt one in the ensemble. 
        6. enr_em_end_lst:<int>: last one in the ensemble. 
        7. keywords:<dict>: extra inputs:
        
        Returns:
        ==========================================================
        1. self.nset:<int>: number of data set. 
 
        """
        data_lst=[]
        flnm_lst=[]
        tracerinfo_flnm_lst=[]
        diaginfo_flnm_lst=[]
        ref_lst=[]
        nset=len(enr_em_st_lst)
        
        for iset  in range(nset):
            
            enr_em_st=enr_em_st_lst[iset]
            enr_em_end=enr_em_end_lst[iset]
            enr_yyyy=enr_yst_lst[iset]
            enr_step=enr_step_lst[iset]
            enr_doy=enr_dst_lst[iset]
            exts_dict=self.construct_enr_filename_dict(\
                yyyy, dd, mm,\
                    enr_yyyy, \
                    enr_doy, \
                    enr_step, \
                    enr_em_st, \
                    enr_em_end)
            
            
            full_name=self.construct_filename(**exts_dict)
            
            fdiag=self.construct_diaginfo_filename(**exts_dict)
            
            ftracer=self.construct_tracerinfo_filename(**exts_dict)
        
            ref_lst.append(iset)
            flnm_lst.append(full_name)
            tracerinfo_flnm_lst.append(ftracer)
            diaginfo_flnm_lst.append(fdiag)
        
                
            bp_lst=read_bpch_to_data_list(full_name, ftracerinfo=ftracer, \
                                              fdiaginfo=fdiag,\
                                              categorys=self.categorys,tracers=self.tracers, \
                                              taus=self.taus, tranames=self.tranames)
            
            
            
            data_lst.append(bp_lst)
            
        # ##L: end loop on set 
        
        self.data_lst=data_lst
        self.flnm_lst=ref_lst
        self.tracerinfo_flnm_lst=tracerinfo_flnm_lst
        self.diaginfo_flnm_lst=diaginfo_flnm_lst
        self.ref_lst=ref_lst
        self.nset=len(ref_lst)
        
        
        
        return self.nset

    
    
    def read_enr_file_profile_list(self,yyyy, mm, dd, \
                                       olon, olat,\
                                       enr_yst_lst,\
                                       enr_dst_lst,\
                                       enr_step_lst,\
                                       enr_em_st_lst,\
                                       enr_em_end_lst,\
                                       **keywords):
        
        """
        Read outputs from ensemble runs to list of ctm_profile_cl
        
        Input:
        ============================================================
        1. yyyy, mm, dd:<int>: year month day 
        2. olon:<array>: longitude
        3. olat:<array>:latitude 
        4. enr_yyyy_lst:<int>: year for ensemble run start.
        5. enr_doy_lst:<int>:  doy for ensemble run start.
        6. enr_step_lst:<int>: step for ensemble run.
        7. enr_em_st_st:<int>: firt one in the ensemble. 
        8. enr_em_end_lst:<int>: last one in the ensemble. 
        9. keywords:<dict>: extra inputs:
        
        Returns:
        ==========================================================
        1. self.nset:<int>: number of data set. 
 
        """
        data_lst=[]
        flnm_lst=[]
        tracerinfo_flnm_lst=[]
        diaginfo_flnm_lst=[]
        ref_lst=[]
        nset=len(enr_em_st_lst)
        
        for iset  in range(nset):
            
            enr_em_st=enr_em_st_lst[iset]
            enr_em_end=enr_em_end_lst[iset]
            enr_yyyy=enr_yst_lst[iset]
            enr_step=enr_step_lst[iset]
            enr_doy=enr_dst_lst[iset]
            exts_dict=self.construct_enr_filename_dict(\
                yyyy, dd, mm,\
                    enr_yyyy, \
                    enr_doy, \
                    enr_step, \
                    enr_em_st, \
                    enr_em_end)
            
            
            full_name=self.construct_filename(**exts_dict)
            
            
            fdiag=self.construct_diaginfo_filename(**exts_dict)
            
            ftracer=self.construct_tracerinfo_filename(**exts_dict)
            
            ref_lst.append(iset)
            flnm_lst.append(full_name)
            tracerinfo_flnm_lst.append(ftracer)
            diaginfo_flnm_lst.append(fdiag)
            
            bp_lst=read_bpch_to_profile_list(olon, \
                                                 olat,\
                                                 full_name, \
                                                 ftracerinfo=ftracer, \
                                                 fdiaginfo=fdiag,\
                                                 categorys=self.categorys,tracers=self.tracers, \
                                                 taus=self.taus, tranames=self.tranames)
            
            
            
            data_lst.append(bp_lst)
        
        # ##L: end loop on set 
            
        self.data_lst=data_lst
        self.flnm_lst=ref_lst
        self.tracerinfo_flnm_lst=tracerinfo_flnm_lst
        self.diaginfo_flnm_lst=diaginfo_flnm_lst
        self.ref_lst=ref_lst
        self.nset=len(ref_lst)
        
        
        
        return self.nset

    

    def get_data(self,traname=None, \
                     category=None, \
                     tracer=None):
        
        """
        select the data from data_lst
        
        
        Inputs:
        =======================================
        1. traname:<list, t:str>: tracer name
        2. category:<list, t:str>: category 
        3. tracerID:<list, t:int>: tracer ID 
        

        Returns:
        ==========================================
        1. traname_lst:<list, t:str>: trace name 
        2. category_lst:<list, t:str>: category 
        3. tid_lst:<list, t:int>: tracer ID 
        4. sid_lst:<list, t:int>:  Index of the data in the list
        5. found_data_lst:<list, t:class bpdata or profile>: 

        """

        found_data_lst=[]
        idx=0
        traname_lst=[]
        category_lst=[]
        tid_lst=[]
        sid_lst=[]
        
        for iref in self.ref_lst:
            bp_lst=self.data_lst[idx]
            for bpdata in bp_lst:
                if (bpdata.is_matched(traname=traname, \
                                          tracer_id=tracer,\
                                          category=category)):
                    
                    
                    found_data_lst.append(bpdata.data)
                    traname=bpdata.get_attr('traname')
                    traID=bpdata.get_attr('tracer_id')
                    category=bpdata.get_attr('category')
                    
                    traname_lst.append(traname)
                    category_lst.append(category)
                    tid_lst.append(traID)
                    sid_lst.append(idx)
            
                    
                    
            idx=idx+1
        
        return traname_lst, category_lst, \
            tid_lst, sid_lst, found_data_lst
    
            
    
    def compute_mod_pres(self,traname='PSURF', \
                             use_sp=True,\
                             geos_ver=5, \
                             use_reduced=1, \
                             **keywords):
        
        """
        
        compute or fetch pressure at model level centres 
      
        Inputs:
        -------------------------------------
        1. traname:<str>: traname for pressure or sp 
        2. use_sp:<T/F>: True, surface pressure will be used to calculate pressures
        3. geos_ver:<int>: GEOS-Chem version
        4. use_reduced:<int>: 1 means reduced model vertical grid is used  

        keywords:
        =========================================
        ---Reserved keys:
        ---category:<str>: category of the surface pressure or pressure
        ---tracer:<int>: tracer ID of surface pressure or model pressure
        ---tau:<float>: time of surface pressure or model pressure
        
        
        
        Returns:
        ----------------------------------
        1. mod_res: <array,(nx, [ny],  [lx])>: pressure at model level or surface pressure (see Notes 1)
        
        Notes:
        -----------------------------------------------------------
        1. when use_sp is set to False, the return will be the surface pressure
        
        
        """
        
        
        category=None
        if ('category' in keywords):
            category=keywords['category']
        
        tracer=None
        if ('tracer' in keywords):
            tracer=keywords['tracer']
            
        
        tau=None
        
        if ('tau' in keywords):
            tracer=keywords['tau']
            
        
        traname_lst, category_lst, \
            tid_lst, sid_lst, found_data_lst=self.get_data(traname=traname, category=category, tracer=tracer)
        
        if (len(found_data_lst)>0):
            gp=found_data_lst[0]
        else:
            msg='No data found for:'+traname
            msm.show_err_msg(msg)
            return None
        
        if (use_sp):
        
            if (geos_ver==5):
                ap, bp=pm.get_geos5_ap_bp(use_reduced)
            
            elif (geos_ver==4):
                ap, bp=pm.get_geos4_ap_bp(use_reduced)
        
    
            
            mod_pres=pm.get_mod_pres(gp, ap, bp)
            mod_pres=npy.squeeze(mod_pres)
                
            return  mod_pres
        
        else:
            return gp
        
        
    
    


        
#=========<<< TESTS >>>=================



if (__name__=='__main__'):
    
    datapath='/home/lfeng/local_disk/otool_data/enkf_output/2009/'
    flnm='ts_satellite.ST000.EN0001-EN0081.20090101.bpch'
    flnm=datapath+flnm
    print flnm
    
    
    fdiaginfo='diaginfo.ST000.EN0001-EN0081.dat'
    fdiaginfo=datapath+fdiaginfo

    ftracerinfo='tracerinfo.ST000.EN0001-EN0081.dat'
    ftracerinfo=datapath+ftracerinfo
    
    print 'read to data'
    
    data_lst=read_bpch_to_data_list(flnm, ftracerinfo, \
                                       fdiaginfo,\
                                       categorys=None,tracers=None, \
                                       taus=None, tranames=['CO2']\
                                       )
    
    olon=npy.arange(-60, 60, 20)
    olat=npy.arange(-60, 60, 20)
    
    print 'read to profile'
    
    prof_lst=read_bpch_to_profile_list(olon, olat, \
                                           flnm, ftracerinfo, \
                                           fdiaginfo,\
                                           categorys=None,tracers=None, \
                                           taus=None, tranames=['CO2']\
                                           )
    

    print 'Number of data:', len(data_lst)
    data=data_lst[0]
    
    print data.dims
    print data.ndim
    print data.name

    tau0=data.get_attr('tau0')
    print tau0
    
    tau1=data.get_attr('tau1')
    print tau1

    ggrid=data.ctm_grd
    lons=ggrid.get_lon()
    print lons

    lats=ggrid.get_lat()
    print lats
    

    print ggrid.lonres

    print ggrid.ix, ggrid.jx, ggrid.lx 
    
    print ggrid.latres

        
    
    lon_edge=ggrid.get_lon_edge()
    lat_edge=ggrid.get_lat_edge()

    print npy.size(lon_edge)
    
    print lon_edge
    
    
    print npy.size(lat_edge)
    print lat_edge


    print 'profiles--------------'
    
    for data in prof_lst:
        print 'name:', data.name
        print 'shape:', data.dims
        print 'dim:', data.ndim
        tracer_id=data.get_attr('tracer_id')
        print 'tid:', tracer_id
        
        traname=data.get_attr('fullname')
        print 'tname:', traname
        
        category=data.get_attr('category')
        print 'category:', category
        
        tau0=data.get_attr('tau0')
        yyyy, mm, dd, hh, mi, sec=tm.tau_to_time_array(tau0)
        print 'Starting time:', yyyy, mm, dd, hh, mi, sec
        tau0=data.get_attr('tau1')
        yyyy, mm, dd, hh, mi, sec=tm.tau_to_time_array(tau0)
        
        print  'End time:', yyyy, mm, dd, hh, mi, sec
    
  
  
    
    prof_lst=read_bpch_to_profile_list(olon, olat, \
                                           flnm, ftracerinfo, \
                                           fdiaginfo,\
                                           categorys=[oob.fill_val_str],tracers=[1], \
                                           taus=None, tranames=[oob.fill_val_str]\
                                           )
    
    
    
  

            
    
    print 'Number of data:', len(prof_lst)
    for data in prof_lst:
        print 'name:', data.name
        print 'shape:', data.dims
        print 'dim:', data.ndim
        tracer_id=data.get_attr('tracer_id')
        print 'tid:', tracer_id
        
        traname=data.get_attr('fullname')
        print 'tname:', traname
        

        category=data.get_attr('category')
        print 'category:', category
        
    print '-----------last ------------------'
    
    prof_lst=read_bpch_to_profile_list(olon, olat, \
                                           flnm, ftracerinfo, \
                                           fdiaginfo,\
                                           categorys=['', ''],\
                                           tracers=[oob.fill_val_int,oob.fill_val_int], \
                                           taus=[oob.fill_val, oob.fill_val], \
                                           tranames=['CO2', 'PS']\
                                           )
    
    
    
    
    
            
    
    print 'Number of data:', len(prof_lst)
    for data in prof_lst:
        if (data.name=='PS'):
            ## surface pressure
            
            print data.dims
            ps=data.data
            ps=npy.squeeze(ps)
            
            
            ap, bp=pres_m.get_geos4_ap_bp(reduced_grid=1)
            mod_pres=pres_m.get_mod_pres(ps, ap, bp)
            print mod_pres[0,:]
            print npy.shape(mod_pres)
            
        print ' '
            
        
        print 'name:', data.name
        print 'shape:', data.dims
        print 'dim:', data.ndim
        tracer_id=data.get_attr('tracer_id')
        print 'tid:', tracer_id
        
        traname=data.get_attr('fullname')
        print 'tname:', traname
        

        category=data.get_attr('category')
        print 'category:', category
        
        
