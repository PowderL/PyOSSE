""" Class for variables used for (netcdf) file IO

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.22
    History: v0.95, 2013.02.28 
    
    
    Classes:
    =============================================================
    1. io_var_cl: class for variables used to do file I/O. 
    
    
    Functions
    ============================================================
    1. construct_filename: constructure file name 
    2. create_var_read_list: create a list of IO var to be read from 
       (netcdf) file
    3. create_var_dim_list: create a list of IO var used to define dimension
    4. create_var_write_list: create a list of var to be saved into (netcdf) file
        
"""



 
import numpy as npy
import otool_obj as oob

import message_m as msm

###1 ============= io_var_cl
 
class io_var_cl:
    """
    Class for file io 
    
    Members:
    -----------------------------------------------
    1.  name: <str>:Name of the variable 
    2.  type: <str>: type of the variable
    3.  dimname:<list, t:str>: list if variable  types
    4.  var:   <array>: data 
    5.  varattr: <dict>: attributes
    6.  colno: <int>: column number of the data 
    7.  attr_rd_lst:<list>: the attribute to be read
    8.  dims:<tuple>: shape of the var
    9.  ndim:<integer>: number of the dimension 
   10.  mask_val: <float>: filling value for missing or bad data
   11.  scale:  <float>: scaling value 
   12.  shift:  <float>: shift value 
   13.  grd:    <gp_grd_cl>: grid 
   
   

   Functions:
   --------------------------------------------------
   1. set_attr_rd_list: set attributes to be read
   2. set_dimname: set dimenion names for the variable 
   3. set_name: set variable name 
   4. set_type: set variable type 
   5. set_colno: set column number
   6. get_attr: get attribute 
   7. set_attr: set attribute
   
   8. set_mask: set mask values
   9. set_offset: set offset values
   10. post_process: re-scale and add shift to the data
   11. set_data: set the data set 
   12. __getitem__: get var at given index 
   13. shape: return shape of the var
   14. copy: make a copy of the var
   
   15. set_attr: set attributes
   16. get_attr: get_attributes
   
   
   """
    
    def __init__(self, varname, vartype, vardim, var, \
                     varattr={},  \
                     colno=oob.fill_val_int,\
                     attr_rd_lst=[], \
                     mask_val=oob.fill_val, \
                     scale=1.0,\
                     offset=0.0,\
                     grd=None ):
        
        """ initialization  
        
        Inputs:
        ------------------------------------------
        1.  varname: <str>:Name of the variable 
        2.  vartype: <str>: type of the variable
        3.  vardim:<list, t:str>: list if variable  types
        4.  var:   <array>: data 
        5.  varattr: <dict>: attributes
        6.  colno: <int>: column number of the data 
        7.  attr_rd_lst:<list>: the attribute to be read
        8.  mask_val:<float>: filling for bad or missing values
        9.  scale:<float>: scaling factor 
        10. offset:<float>: offset 
        11. grd:<gp_grid_cl>: grid of the data 

        """

        self.name=varname
        self.type=vartype
        self.dimname=vardim
        self.var=None
        self.dims=None
        self.ndim=None
        self.scale=scale
        self.offset=offset
        self.grd=grd

        
        
        self.set_data(var)
        if (varattr<>None):
            self.varattr=dict(varattr)
        else:
            self.varattr={}
        
        self.colno=colno
        
        self.attr_rd_lst=attr_rd_lst
        
        if ('_fillvalue' in self.varattr):
            pass
        else:
            self.varattr.update({'_fillvalue':mask_val})
        
        if ('units' in self.varattr):
            pass
        else:
            if (varname=='level'):
                self.varattr.update({'units':'hPa'})
            elif (varname=='longitude'):
                self.varattr.update({'units':'degrees_east'})
            elif (varname=='latitude'):
                self.varattr.update({'units':'degrees_north'})


        if ('ot_type' in self.varattr):
            pass
        else:
            self.varattr.update({'ot_type':oob.ot_iovar})
        
            
        
    def set_grid(self, grd):
        """
        setup the grid for the var 
        
        """
        self.grd=grd.copy()
        

        
    def set_attr_rd_list(self,  attr_rd_lst):
        """ set list for attributes to be read
        
        Inputs:
        ---------------------------------------
        1.  attr_rd_lst:<list>: the attribute to be read
        """
        self.attr_rd_lst=attr_rd_lst


    def set_dimname(self,  dimname):
        """ set dimension Names
        
        Inputs:
        ---------------------------------------
        1.  dimname:<list, t:str>: dimension names
        """
        
        self.dimname=dimname
        
        
    def set_name(self, name):
        """ set dimension Names
        
        Inputs:
        ---------------------------------------
        1.  name:<str>: name of the var 
        """
        
        self.name=name

        
    def set_type(self, stype):
        """ set dimension Names
        
        Inputs:
        ---------------------------------------
        1.  stype:<>: var type 
        """
        
        self.type=stype
    
   
    def set_mask(self, mask_val):
        """ set mask value 
        
        Inputs:
        ---------------------------------------
        1.  mask_val:<float>: mask_val 
        """
        self.mask_val=mask_val

    
    def set_scale(self, scale):
        """ set scaling factor 
        
        Inputs:
        ---------------------------------------
        1.  scale:<float>: scaling factor 
        """
        self.scale=scale

    def set_offset(self, offset):
        """ set offset values 
        
        
        Inputs:
        ---------------------------------------
        1. offset:<float>: offset  
        """
        self.offset=offset


    def post_process(self, data, new_maskval):
        
        """ re-scale and add shift to the data 
        
        Inputs:
        ---------------------------------------
        1. data:<array>: data  
        2. new_mask_val:<float>: new masks for data or missing values
        
        Returns:
        --------------------------------------------------
        1. data: <array>: data after been processed 
        """
        
        if (oob.get_ot_type(data)==oob.ot_array):
            usd_idx=npy.where(data<>self.mask_val)
            data[usd_idx]=self.scale*data[usd_idx]+self.offset
            usd_idx=npy.where(data==self.mask_val)
            data[usd_idx]=new_maskval
            
        elif (oob.get_ot_type(data)==oob.ot_list):
            ndata=len(data)
            for idx in range(ndata):
                if (data[idx]==self.mask_val):
                    data[idx]=new_maskval
                else:
                    data[idx]=self.scale*data[idx]+self.offset
        else:
            if (data==self.mask_val):
                data=new_maskval
            else:
                data=self.scale*data+self.offset
        
        return data



        
    def set_colno(self, colno):
        """ set column number 
        
        Inputs:
        ---------------------------------------
        1.  colno:<integer>: column of the variable  in the file 
        """
        self.colno=colno

        
        
    def get_attr(self, attr_name):
        
        """ find the attribute
        
        Inputs:
        ----------------------------------
        1. attr_name:<str>: name of the attribute
        
        Returns:
        ------------------------------------------
        1. attr: <obj>: if the attribute exists. 
        """
        
        if (attr_name in self.varattr):
            return self.varattr[attr_name]
        else:
            try:
                attr=getattr(self.var, attr_name)
                return attr
            except:
                msg=attr_name
                msm.show_err_msg(msg, msm.msm_no_attr)
                return None
        
   
    
    def set_attr(self, attr_name, attr_value):
        
        """ set attribute
        Inputs:
        ----------------------------
        1. attr_name:<str>: name of the attribute
        2. attr_value:<obj>: the value of the attribute 
        """
        
        self.varattr.update({attr_name:attr_value})
    
    def del_attr(self, attr_name):
    
        """ delete  the attribute
        
        Inputs:
        ----------------------------------
        1. attr_name:<str>: name of the attribute
        
        """
        
        if (attr_name in self.varattr):
            del self.varattr[attr_name]
        else:
            msg=attr_name
            msm.show_err_msg(msg, msm.msm_no_attr)
    
    

    
    def set_data(self, data, copy=True):
        """ set the data set  
        
        Inputs:
        ----------------------------
        1. data:<array>: data
        
        """
        
        if (data<>None):
            if (copy):
                self.var=npy.array(data)
                # print 'set_data copy', npy.shape(data)
            else:
                self.var=data
            
            
            self.dims=npy.shape(self.var)
            self.ndim=npy.size(self.dims)
        else:
            self.dims=[]
            self.ndim=0
            self.var=data
        #       print 'var shape', npy.shape(self.var)
            
            
    def __getitem__(self, index):
        """ get value at given index 

        Inputs:
        -----------------------------------------------
        1. index:<integer/tuple>: index

        Returns:
        ---------------------------------------------------
        1. var[index]:<float/array>: values at given index
        
        """
        
        return self.var[index]
    

    def shape(self):
        """
        get shape of the data set 
        Returns:
        ---------------------------------------------------
        1. dims:<tuple/list>: shape of the var
        """
        
        return self.dims
    
    def copy(self):
        """ make a copy 

        Returns:
        1. new_var:<io_var_cl>: a copy of the var

        """
        
        new_var=io_var_cl(self.name, self.type, self.dimname, self.var, \
                              varattr=self.varattr,  \
                              colno=self.colno,\
                              attr_rd_lst=self.attr_rd_lst, \
                              mask_val=oob.fill_val)
        return new_var



    


def construct_filename(filepath, filenameformat, **keywords):
    """ construction file name use the keywords 
    Inputs:
    ------------------------------------------------
    1. filepath:<str>:template for file path 
    2. filenameformat:<str>:template for file name

    3. keywords:<dict>: the parts should be replaced 


    Returns:
    -------------------------------------------------
    sflnm:<str>: full filename 
    

    Notes:
    --------------------
    Keywords should in format of XKEYWORDX
    
    

    
    """
    
    # construct the name 
    sflnm=filenameformat
    spath=filepath
    
    for keyname in keywords:
        keyval=keywords[keyname]
        if (oob.get_ot_type(keyval)<>oob.ot_string):
            # not a string 
            if (keyname=='XYYYYX'):
                keyval=r'%4.4d' % keyval
            elif (keyname=='XMMX'):
                keyval=r'%2.2d' % keyval
            elif (keyname=='XDDX'):
                keyval=r'%2.2d' % keyval
            elif (keyname=='XHHX'):
                keyval=r'%2.2d' % keyval
            elif (keyname=='XDOYX'):
                keyval=r'%3.3d' % keyval
                
            else:
                keyval=str(keyval)
                
        
            
        sflnm=sflnm.replace(keyname, keyval)
        spath=spath.replace(keyname, keyval)
    spath=spath.strip()
    
    if (len(spath)>0):
        sflnm=spath+'/'+sflnm
    
    return sflnm
    
                    
def construct_enr_filename(path, flnm,\
                               yyyy, dd, mm,\
                               enr_yyyy, \
                               enr_doy, \
                               enr_step, \
                               enr_em_st, \
                               enr_em_end):
        
                                        
    """
    construct file name for ensemble runs
    
    Input:
    ============================================
    1. path:<str>: directory
    2. flnm:<str>: name 
    3. enr_yyyy:<int>: year for ensemble run start 
    4. enr_doy:<int>:  doy for ensemble run start
    5. enr_step:<int>: step for ensemble run 
    7. enr_em_st:<int>: firt one in the ensemble 
    8. enr_em_end:<int>: last one in the ensemble 
    9. yyyy, dd, mm:<int>: output time 
    
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
        
        
    full_flnm=construct_filename(path, flnm, **exts_dict)
    return full_flnm

    


def create_var_read_list(vname_lst, attr_read_lst=None):

    """
    
    Make a list for IO var to be read from 
    
    
    Inputs:
    ======================================================
    1. vname_lst: <list, t:str>: list of variables to read in 
    2. attr_read_lst: <list, t:list>: attributes of vars to be read in 
    
    
    Returns:
    ================================================
    1, var_lst:<list, io_var_cl>: list of var 
    
    Notes:
    """
    # S1: Set attribute read list
    
    empty_lst=[]
    nvar=len(vname_lst)
    
    if (attr_read_lst==None):
        var_attr_lst=nvar*[empty_lst]
    else:
        var_attr_lst=attr_read_lst
        

    var_lst=[]

    # S2: creat var for each varname 
    
     
    for ivar in range(nvar):
        # #c: name 
        vname=vname_lst[ivar]
        
        # #c: attribute to be read
        
        vattr=var_attr_lst[ivar]
        
        # #c: create var 
        
        new_var=io_var_cl(vname, 'f', [], None, attr_rd_lst=vattr)
        
        # #c: adding to the list
        
        var_lst.append(new_var)

    
    # loop ivar
    
    
    return var_lst


def create_var_dim_list(vname_lst, vtype_lst=None, vdata_lst=None):
    
    """
    construct io a list for dimension definition 
    
    Inputs:
    ------------------------------------------
    1. vname_lst:<list, t:str>: name of dimensions used in the file 
    2. vtype_lst:<list, t:str>: data type of dimension 
    3. vdata_lst:<list, t:array>: data of dimension 
    
    Outputs:
    ------------------------------------------
    1. var_lst:<list, io_var_cl>: list of dimension  var 
    
    """
    var_lst=[]
    nvar=len(vname_lst)
    # S1: loop over 

    for ivar in range(nvar):
        # #c: Name
        vname=vname_lst[ivar]
        # #c: Type
        
        vtype='f'
        if (vtype_lst<>None):
            vtype=vtype_lst[ivar]

        # #c: data
        vdata=None
        
        if (vdata_lst<>None):
            vdata=vdata_lst[ivar]
    
        # S2 : create var
            
        var= io_var_cl(vname,vtype,  \
                           [vname], \
                           vdata)
        # S3: adding to the new 
        var_lst=var_lst+[var]

    # loop ivar
    
    return var_lst




def create_var_write_list(vname_lst, vdim_lst, vtype_lst=None, vdata_lst=None):
    
    
    """
    construct io avr list for variables to be stored n ensemble fluxes 
    
    Inputs:
    ------------------------------------------
    1. vname_lst:<list, t:str>: list of variable names
    2. vdim_lst:<list, t:list>: list of dimension names for each variable  (see Note 1)
    3. vtype_lst:<list, t:str>: list of variable type
    4. vdata_lst:<list, t:array>: list of data
    
    1. Outputs:
    ------------------------------------------
    var_lst:<list, io_var_cl>: list of IO var 
    
    Notes: 
    ==========================================================
    
    1. if one element in vdim_lst is a string, this element will be converted to a list 
    afterm being split using separator ';' 
    
    
    """
    # S1: set list 
    var_lst=[]
    nvar=len(vname_lst)
    
    for ivar in range(nvar):
        
        # S1: var information 
        
        # #c: Name 
        vname=vname_lst[ivar]
        
        # #c: dimension 
        
        vdim=vdim_lst[ivar]
        
        if (oob.get_ot_type(vdim)==oob.ot_string):
            # single string
            terms=vdim.split(';')
            dim_s=[]
            
            for item in terms:
                dim_s.append(item)
            
            vdim=dim_s
        print 'vdim', vdim
        
        
            
        
        
        # #c: type

        if (vtype_lst<>None):
            vtype=vtype_lst[ivar]
        else:
            # default is float 
            
            vtype='f'
        
        
        # #c: data
        if (vdata_lst<>None):
            vdata=vdata_lst[ivar]
        else:
            vdata=None


        

        # #S2: create val
        
        var= io_var_cl(vname,vtype,  \
                           vdim, \
                           vdata)
        
        # add
        var_lst=var_lst+[var]
    
    # loop ivar
    
    return var_lst
    
    

if (__name__=='__main__'):
    
    print '>>>>1 test creat var'

    
    
    vlon= io_var_cl('longitude','f',  ['longitude'], None)
    vlon.set_attr_rd_list(['units'])
    vlat= io_var_cl('latitude','f',  ['latititude'], None)
    vlat.set_attr_rd_list(['units'])
    
    vtime=io_var_cl('time', 'f', ['time'], None)
    vtime.set_attr_rd_list(['units'])

    print '>>> test create var list'
    
    
    # t1: create read list 
    vname_lst=['lon', 'lat']
    vtype_lst=['f', 'f']
    vattr_lst=['units', 'units']
    
    var_lst=create_var_read_list(vname_lst, vattr_lst)
    
    # t2: create dimension
    lon=npy.arange(-180, 180, 5.0)
    lat=npy.arange(-90, 90, 4.0)
    vdata_lst=[lon, lat]
    nlon=npy.size(lon)
    nlat=npy.size(lat)
    
    dim_lst=create_var_dim_list(vname_lst, vtype_lst, vdata_lst)
    
    
    vname_lst=['pres', 'area']
    vtype_lst=['f', 'f']
    vdim_lst=['lon;lat', 'lon;lat']
    pres=npy.zeros([nlon,nlat], float)
    area=npy.zeros([nlon,nlat], float)
    vdata=[pres, area]
    
    
    var_out_lst=create_var_write_list(vname_lst, vdim_lst, vtype_lst, vdata_lst)
    
    
    
    
    
    
  
    
