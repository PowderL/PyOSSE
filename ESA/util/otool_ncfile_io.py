"""
module for file IO. It contains codes for reading and writing netcdf files;

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.05.18
    History: v0.95, 2013.02.15
    

Classes:
=============================================================
1. ncfile_desc_cl: class for ncfile access



Functions:
==============================================================
1. ncf_read: read variables into arrays  from netcdf file.
2. ncf_write: create a netcdf file, and write variables into it.
3. ncf_write_byvar: create a netcdf file, and write variables given 
---as a list of class io_var_cl. 
4. ncf_dims_def_byvar: create a netcdf file, and define dimensions 
---using a list of io_var_cl.
5. ncf_dims_add_byvar: add dimensions using a list of io_var_cl 
---to a netcdf file. 
6. ncf_dims_def:  create a netcdf file, and define dimensions.
7. ncf_data_write: write data into a netcdf file.
8. ncf_save_var: save var into a netcdf file.
9. ncf_data_write: write data into a netcdf file. 
10. ncf_save_var: save var into a netcdf file.


"""

import netCDF4 as snf 
import numpy as npy
import otool_obj as oob
import otool_var_io as ovio
import message_m as msm


#===< FUNCTIONS >====      
def ncf_get_attr(ncflnm, varname, attrname):
    """ read variables into arrays from a netcdf file 
    Inputs:
    ------------------------------------------------
    1. ncflnm:<str>: file name.
    2. varname:<str>: var name to be read.
    3. attrname:<str>: attribute name 
    
    Returns:
    -----------------------------------------
    1. attr_val:<obj>: attribute of variable 
    
    Notes:
    1. if varname=='', global attribute will be seek 
    """
    
    ncfile=snf.Dataset(ncflnm)
    if (varname<>''):
        # #c: variable attribute 
        var=ncfile.variables[varname]
        attr=getattr(var, attrname)
    else:
        # #c: global  attribute
        attr=getattr(ncfile, attrname)
    
    ncfile.close()
        
    return attr
    


def ncf_read(ncflnm, varnames):
    """ read variables into arrays from a netcdf file 
    Inputs:
    ------------------------------------------------
    1. ncflnm:<str>: file name.
    2. varnames:<list, t:str>: var name to be read.

    Returns:
    -----------------------------------------------
    1. outvars:<list, t:array>: list of the variables read from netcdf file. 

    """

    # S1 open file 
    
    ncfile=snf.Dataset(ncflnm)
    
    # S2 read data
    
    tc=list()
    if (type(tc)==type(varnames)):
        # name in a list 
        outvars=list()
        nvar=len(varnames)
        for ivar in range(nvar):
            vname=varnames[ivar]
            var=ncfile.variables[vname]
            var=npy.array(var)
            var=npy.squeeze(var)
            outvars.append(var)
    else:
        # a single name 
        outvars=ncfile.variables[varnames]
        outvars=npy.array(outvars)
        outvars=npy.squeeze(outvars)
    ncfile.close()
    
    return outvars

def ncf_write(ncflnm, dimNames, dimTypes, dimVars, \
                  varNames, varTypes, varDatas):
    
    """ write data into a netcdf file 
    Inputs:
    --------------------------------------------------------------
    1. ncflnm:<str>: file name 
    2. dimNames:<list, t:str>: dimension name 
    3. dimTypes:<list, t:str>:dimension type 
    4. dimvars:<list, t:array>:dimension data
    
    5. varNames:<list, t:str>: var name
    6. varTypes:<list, t:str>:var type
    7. varDatas:<list, t:array>: data
    """
    
    # S1 create dimension 
    ncf_dims_def(ncflnm,dimNames,dimTypes, dimVars)
    

    # S2 create var and write data 

    tc=list()
    if (type(tc)==type(varNames)):
        
        # list 

        nvar=len(varNames)
        for ivar in range(nvar):
            varName=varNames[ivar]
            varType=varTypes[ivar]
            varData=varDatas[ivar]
            ncf_data_write(ncflnm,dimNames,varName, varType, varData)
    else:
        # single 
        
        varName=varNames
        varType=varTypes
        varData=varDatas
        ncf_data_write(ncflnm,dimNames,varName, varType, varData)

def ncf_write_byvar(ncflnm, dimNames, dimTypes, dimVars, varinfos):

    """ write data contained by io_var_cl into a netcdf file 
    
    Inputs:
    --------------------------------------------------------------
    1. ncflnm:<str>: file name 
    2. dimNames:<list, t:str>: dimension name 
    3. dimTypes:<list, t:str>:dimension type 
    4. dimvars:<list, t:array>:dimension data
    
    5. varinfos:<list, t:geos_info>: list of vars
    """
    
    # S1 create dimension
    
    ncf_dims_def(ncflnm,dimNames,dimTypes, dimVars)
    
    # S2 create and write var 

    tc=list()
    
    if (type(tc)==type(varinfos)):
        # list 
        nvar=len(varinfos)
        for varinfo in varinfos:
            varName=varinfo.name
            varType=varinfo.type
            varData=varinfo.var
            vardimnames=varinfo.dimname
            varattr=varinfo.varattr
            ncf_data_write(ncflnm,vardimnames,\
                               varName, varType, varData,varattr)
                
    else:
        # single 
        
        varinfo=varinfos
        varName=varinfo.name
        varType=varinfo.type
        varData=varinfo.var
        vardimnames=varinfo.dimname
        varattr=varinfo.varattr
        
        ncf_data_write(ncflnm,vardimnames,\
                           varName, varType, varData,varattr)
    


                               


        

def ncf_dims_def_byvar(filename,var_lst):
    
    """
    Create dimensions from list of var 
    Inputs:
    ------------------------------------
    1. filename:<str>: name of netCDF file 
    2. var_list:<list, t:io_var_cl>: list of geos_var
    """
    
    
    # S1 open file to write 


    file_nc=snf.Dataset(filename, 'w')


    nsizes=list()
    
    if (type(nsizes)==type(var_lst)):
        ndim=len(var_lst)
    else:
        ndim=1
        var_lst=[var_lst]

    
    
    nsizes=list()
    vcontain=list()

    # S2 define dimension, and write variables 
    
    for id in range(ndim):
        gvar=var_lst[id]
        
        dvar=gvar.var
        
        dname=gvar.name
        
        # print 'dname, dvar:', dname, npy.shape(dvar)
        
        dsize=npy.size(dvar)
        dim_1d=(dname,)
        ## T1 create dimension 
        file_nc.createDimension(dname,dsize)
        st=gvar.type
        ## T2 create variable 
        
        var1d=file_nc.createVariable(dname, st, dim_1d)
        nsizes.append(dsize)
        ## T3 assign data 
        for i in range(dsize):
            var1d[i]=float(dvar[i])
            # print var1d[i]
        ## T4 assign attributes
        
        for attr_name in gvar.varattr:
            attr_val=gvar.varattr[attr_name]
            setattr(var1d, attr_name,attr_val)
        
        vcontain.append(var1d)
        file_nc.sync()

    file_nc.close()




def ncf_dims_add_byvar(filename,var_lst):
    
    """
    Add dimensions from list of var 
    Inputs:
    ------------------------------------
    1. filename:<str>: name of netCDF file 
    2. var_list:<list, t:io_var_cl>: list of geos_var
    """
    
    
    # S1 open file to append

    
    file_nc=snf.Dataset(filename, 'a')
    
    nsizes=list()
    
    if (type(nsizes)==type(var_lst)):
        ndim=len(var_lst)
    else:
        ndim=1
        var_lst=[var_lst]
        
    vcontain=list()

    # S2 define dimension, and write variables 
    
    for id in range(ndim):
        gvar=var_lst[id]
        
        dvar=gvar.var
        dname=gvar.name
        dsize=npy.size(dvar)
        dim_1d=(dname,)
        ## T1 create dimension 
        file_nc.createDimension(dname,dsize)
        st=gvar.type
        ## T2 create variable 
        
        var1d=file_nc.createVariable(dname, st, dim_1d)
        nsizes.append(dsize)
        ## T3 assign data 
        for i in range(dsize):
            var1d[i]=float(dvar[i])
            # print var1d[i]
        ## T4 assign attributes
        
        for attr_name in gvar.varattr:
            attr_val=gvar.varattr[attr_name]
            setattr(var1d, attr_name,attr_val)
        
        vcontain.append(var1d)
        file_nc.sync()

    file_nc.close()



def ncf_dims_def(filename,dimNames,dimTypes, dimVars):
    """  Define dimensions into a netcdf file 
    
    
    Inputs:
    --------------------------------------------------------------
    1. ncflnm:<str>: file name 
    2. dimNames:<list, t:str>: dimension name 
    3. dimTypes:<list, t:str>:dimension type 
    4. dimvars:<list, t:array>:dimension data
    
 
    """
    # S1 open file to write 
    tx=[]
    
    file_nc=snf.Dataset(filename, 'w')
    
    # S2 create dimension 
    if (type(tx)==type(dimNames)):
        ndim=len(dimNames)
    else:
        # single 
        ndim=1
        dimNames=[dimNames]
        dimTypes=[dimTypes]
        dimVars=[dimVars]
    
    nsizes=list()
    vcontain=list()
    
    
    for id in range(ndim):

        dvar=dimVars[id]
        dname=dimNames[id]
        dsize=size(dvar)
        dim_1d=(dname,)
        ## T1  create dimension 
        file_nc.createDimension(dname,dsize)
        
        
        ## T2 create variable for the dimension
        
        st=dimTypes[id]
        var1d=file_nc.createVariable(dname, st, dim_1d)

        nsizes.append(dsize)
        
        for i in range(dsize):
            var1d[i]=float(dvar[i])
        ## T3 add some attributes
            
        
        setattr(var1d, '_fillValue', -999.0)
        
        if (dname=='level'):
            setattr(var1d, 'units', 'hPa')
        elif (dname=='longitude'):
            setattr(var1d, 'units',  "degrees_east")
        elif (dname=='latitude'):
            setattr(var1d, 'units',  "degrees_north")
            
        vcontain.append(var1d)
        file_nc.sync()
    file_nc.close()





def ncf_data_write(filename,dimNames,varName, \
                       varType, varData, varattr=None):
    """ Write variable into netcdf file
    
    Inputs:
    --------------------------------------
    1. filename:<str>: netcdf file name
    2. dimNames:<list, t:str>: name of the var dimension
    3. varName:<str>: name of the var
    4. varData:<array,>: data 
    5. varType:<str>: type of the data 
    6. varattr:<dict>: attributes 
    """
    
    # S1 open file to append
    file_nc=snf.Dataset(filename, 'a')
    
    # S2 create variable 
    
    ndim=len(dimNames)
    
    print 'ncf_data_write--var name:', varName
    print 'ncf_data_write--var shape:', npy.shape(varData)
    
    
    if (ndim>1):
        dim_nd=tuple(dimNames)
    else:
        dim_nd=(dimNames[0],)
    
    
    
    cdf_var=file_nc.createVariable(varName, varType, dim_nd)
    
    # S3 assign data to variable 
    
    var=varData
    nsizes=npy.shape(var)
    idx=slice(None)
    cdf_var[idx]=var[idx]
    
    if (varattr<>None):
        for attrname in varattr:
            
            attrval=varattr[attrname]
            
            setattr(cdf_var, attrname, attrval)
            
    del var
    
    file_nc.close()



def ncf_save_var(ncflnm, var_lst, \
                     dim_lst=None, create_new=False, glb_attr=None):

    """
    save var info a file 
    Inputs:
    -------------------------------------------------------
    1. ncflnm:<str>: file name 
    2. var_lst:<list, t:io_var_cl>: datas
    3. dim_lst:<list, t:io_var_cl>: dimensions
    4. create_new:<T/F>: True if a new file will be created
    5. glb_attr:<dict>: global attributes
    
    
    """

    # S1 create file and dimension if necessary 
    
    if (dim_lst<>None):
        if (create_new):
            ncf_dims_def_byvar(ncflnm, dim_lst)
        else:
            ncf_dims_add_byvar(ncflnm, dim_lst)
    
    tc=list()

    # S2 create var and save data
    

    if (type(tc)==type(var_lst)):
        # list 
        nvar=len(var_lst)
        
        for varinfo in var_lst:
            
            # #c: loop over each variable 
            
            varName=varinfo.name
            varType=varinfo.type
            varData=varinfo.var
            vardimnames=varinfo.dimname
            varattr=varinfo.varattr
            
            # #c: write data into netcdf file 
            
            
            ncf_data_write(ncflnm,vardimnames,\
                               varName, \
                               varType, varData,varattr=varattr)
    else:
        # single 
        
        
        varinfo=var_lst
        
        
        varName=varinfo.name
        varType=varinfo.type
        varData=varinfo.var
        vardimnames=varinfo.dimname
        varattr=varinfo.varattr
        
        ncf_data_write(ncflnm,vardimnames,\
                           varName, \
                           varType, varData,varattr)
    
    # S3: assign global attributes
    if (glb_attr<>None):
        ncf_set_glb_attr(ncflnm, glb_attr)
        
    
def ncf_set_glb_attr(ncflnm, glb_attr):
    
    """
    
    Inputs:
    ------------------------------------------------
    1. ncflnm:<str>: file name.
    2. glb_attr:<dict>: global attributes.
    """
    
    ncfile=snf.Dataset(ncflnm, 'a')

    for attr_name in glb_attr:
        attr_val=glb_attr[attr_name]
        setattr(ncfile, attr_name, attr_val)

    ncfile.close()
    

def ncf_read_glb_attr(ncflnm, attr_name):
    
    """
    
    Inputs:
    ------------------------------------------------
    1. ncflnm:<str>: file name.
    2. attr_name:<str>: global attribute name .

    Returns:
    -----------------------------------
    1. attr_val:<obj>: global attributes
    """
    
    ncfile=snf.Dataset(ncflnm)

    attr_val=getattr(ncfile, attr_name)
    ncfile.close()
    
    return attr_val



        

def ncf_read_var(ncflnm, var_lst, to_array=True):
    
    """ read variables into var list from a netcdf file 
    
    Inputs:
    ------------------------------------------------
    1. ncflnm:<str>: file name.
    2. var_lst:<list, t:str>: var name to be read.

    Returns:
    -----------------------------------------------
    1. nread:<obj>: 
    ---if to_array=True, nread will be the number of read
    ---if to_array=False, nread will be the filehandle 
    
 
    
    
    """

    # S1 open file 
    
    ncfile=snf.Dataset(ncflnm)
    nread=0
    # S2 read data
    
    tc=list()
    if (type(tc)==type(var_lst)):
        # list 
        
        nvar=len(var_lst)
        for ivar in range(nvar):
            gvar=var_lst[ivar]
            vname=gvar.name
            print 'netcdf read var:', ivar, vname
            
            fvar=ncfile.variables[vname]
            dimNames=fvar.dimensions
            gvar.set_dimname(dimNames)
            
            
            for attr_name in gvar.attr_rd_lst:
                try:
                    attr_val=getattr(fvar, attr_name)
                    gvar.set_attr(attr_name, attr_val)
                    
                except:
                    
                    msm.show_err_msg('No attribute '+attr_name+' in'+vname) 
                    
                
            if (to_array):
                gvar.set_data(fvar)
            else:
                gvar.set_data(fvar, False)
            
            nread=nread+1
            
    else:
        
        # a single name 
        gvar=var_lst
        vname=gvar.name
        fvar=ncfile.variables[vname]
        dimNames=fvar.dimensions
        gvar.set_dimname(dimNames)
        
        for attr_name in gvar.attr_rd_lst:
            attr_val=getattr(fvar, attr_name)
            gvar.set_attr(attr_name, attr_val)
            
        gvar.set_data(fvar)
        nread=1
        
    if (to_array):
        ncfile.close()
        return nread
    
    else:
        return ncfile
    



#===< CLASSES >====

class ncfile_desc_cl:
    
    """ class for netcdf file access
    
    Members:
    -----------------------------------------------
    1. mask_val:<float>: filling for bad or missing data
    2. var_lst:<list, t:str>: list of var to be read 
    3. varname_dict:<dict>: translation between  name in netcdf file and name in the code. 
    4. varname_lst:<list>: list of var name
    
    
    5. read_to_array:<T/F>: if ture, data will be read from file as array (full loaded), otherwise
    --- data will be just stored as pointer
    
    6. nread:<obj>: number of vars, or the file handle 
    7. attr_dict:<dict>: attribute dictionarys
    8. grd:<ctm_grd_cl>: grid 
    10. ref:<int/str>: value to identify the object
    11. flnm:<str>: filename format
    12. path:<str>: path format
    13. full_flnm:<str>: the actual file name 
    
    
    Functions
    --------------------------------------------
    1.  __init__: initialization
    2. read_var_from_file: read var into the class 
    3. append_var: append var to the var list of the class
    
    4. get_lon: get longitude if existing 
    5. get_lat: get latitude if existing 
    6. get_time: get time if existing 
    7. get_data: get var data 
    8. get_var: get var

    9. set_file_path:  set path for data file

    10. set_filename_format: set template for filename 

    11. construct_filename: make file name 
    12. get_attr: get attribute
    13. set_attr: set attribute
    14. close_file: close file 
    15. set_grid: set grid 
    
    

    Attributes (reserved words)
    -----------------------------------------------------------
    1. flnm: <str>: full file name 
    2. filepath:  <str>: path of the file name 
    3. flnmformat:  <str>: template for filename
    
    """
    
    def __init__(self,\
                     flnm="",\
                     ref=0,\
                     path="",\
                     var_lst=[],\
                     varname_dict={},\
                     mask_val=oob.fill_val, \
                     read_to_array=True,\
                     fio_keywords={},\
                     **keywords):
        

        """ read txt table files into var and return a record table
        
        Inputs:
        -----------------------------------------
        1. flnm:<str>: filename 
        2. datapath:<str>: data path for the file 
        
        3. var_lst:<list, t:io_var_cl_>:list of var to be read
        4. varnames_dict:<dict>:translation between default var name and name 
        in netcdf file
        3. mask_val:<float>: filling for bad or missing data
        5. read_to_array:<T/F>: if ture, the data will be fully read and 
        --- converted into numpy array 
        6. fio_keywords:<dict>: extra inputs for file I/O

        7. keywords:<dict>: attributes

        
        """
        self.flnm=flnm
        self.path=path
        self.full_filename=flnm
        
        self.mask_val=mask_val
        self.var_lst=var_lst
        self.varname_lst=[]
        if (self.var_lst<>None):
            for var in self.var_lst:
                self.varname_lst.append(var.name)
        
        
        self.varname_dict=varname_dict
        
        self.read_to_array=read_to_array
        self.nread=0
        self.fio_keywords=fio_keywords
        self.grd=None
        
        self.attr_dict={}
        
        if ('ot_type' in self.attr_dict):
            pass
        else:
            self.attr_dict.update({'ot_type':oob.ot_ncfdesc})
        
        # update attributes
            
        self.ref=mask_val
        for keyname in keywords:
            keyval=keywords[keyname]
            self.set_attr(keyname, keyval)

        
    
    def read_var_from_file(self, flnm=""):
        
        """ read txt table files into var and return a record table
        
        Inputs:
        -----------------------------------------
        1. flnm:<str>: file name 
        
        2. var_lst:<list, t:io_var_cl_>:list of var to be read
        
        Returns:
        -----------------------------------------
        1. nread:<int>:number of variables been read. 

        Notes:
        1. 
        """
        
        if (flnm==""):
            flnm=self.full_filename
        else:
            self.set_attr('filename', flnm)
        
        self.nread=ncf_read_var(flnm, self.var_lst, self.read_to_array)
        
        return self.nread
    
    def append_var(self, var):
        
        """append var to the var_lst of the class

        Inputs:
        ----------------------------------------
        1. var:<io_var_cl>: variable 
                
        """
        self.var_lst.append(var)
        self.varname_lst.append(var.name)
        
    
    def set_grid(self, grd):
        """
        setup the grid for the var 
        
        """
        self.grd=grd
        
    def get_grid(self):
        
        """
        get the grid 
        Returns:
        ---------------------------------------------
        1. self.grd:<ctm_grid_cl>: grid 
        
        
        """
        return self.grd

        
        
    def get_lon(self, varname='lon'):
        
        """get longitude 
        Inputs:
        ----------------------------------------
        1. varname:<str>: name of variable longitude 
        
        Returns:
        -----------------------------
        1. data:<array/obj>: longitude
        
        """
        
        
        data=self.get_data(varname)
        return data
  
        
        
    def get_lat(self, varname='lat'):
    
        """get latitude
        
        Inputs:
        ----------------------------------------
        1. varname:<str>: name of variable latitude 
        
        Returns:
        -----------------------------
        1. data:<array/obj>: latitude
        
        """
        
        
        data=self.get_data(varname)
        return data
  
        
    
    def get_time(self, varname='time'):
        """get time 
        
        Inputs:
        ----------------------------------------
        1. varname:<str>: name of variable time 
        
        Returns:
        -----------------------------
        1. data:<array/obj>: time
        
        """
        
        data=self.get_data(varname)
        return data
    
                       
    def get_data(self, varname):
        
        """get value
        
        Inputs:
        ----------------------------------------
        1. varname:<str>: name of variable
        
        Returns:
        -----------------------------
        1. data:<array/obj>: data of the variable
        
        """
        
        data=None
        # print varname
        
        if (varname in self.varname_dict):
            # translate from varname to name used 
            # to read and store data
            varname=self.varname_dict[varname]
            
        for gvar in self.var_lst:
            if (gvar.name==varname):
                data=gvar.var
                break
            
        return data
    
    
    def set_file_path(self, path):
        """ set path for data file 

        Inputs:
        ----------------------------------------
        1. path:<str>: file path
        
        """
        
        self.set_attr('filepath', path)
        self.path=path
    
    def set_filename(self, flnm):
        """ set path for data file 

        Inputs:
        ----------------------------------------
        1. flnm:<str>: file name
        
        """
        
        self.set_attr('filename', flnm)
        
        
    
    def set_filename_format(self, flnm_format):
        """set template for filename 

        Inputs:
        ---------------------------
        1. flnm_format:<str>: file name 
        
        """
        
        self.set_attr('flnmformat', flnm_format)
        self.flnm=flnm_format
        
    
    def construct_filename(self, **keywords):
        
        """ construction file name use the keywords """
        
        sflnm=self.get_attr('flnmformat')
        if ('filepath' in self.attr_dict):
            filepath=self.get_attr('filepath')
        else:
            filepath=""
        
        # construct the name 
        
        
        sflnm=ovio.construct_filename(filepath, sflnm, **keywords)
        
        return sflnm
    
                    
    def set_attr(self, attr_name, value):
        
        """  assign one attribute to axis
        
        Inputs:
        -------------------------------------------------------
        1. attr_name: <str>: attribute name
        2. value: <obj>:attribute value
        
        """
        
        self.attr_dict.update({attr_name:value})
        if (attr_name=='filename'):
            self.full_filename=value
        
        
    def __getitem__(self, index):
        """ overide []
        Inputs:
        -------------------------------
        1. index:<str/int>: var name or the index in the var_lst 
        
        Returns:
        ----------------------------
        1. gvar:<io_var_lst>: variable for given name or index
        

        """

        if (oob.get_ot_type(index)==oob.ot_string):
            # translate from varname to name used to read in data
            
            if (index in self.varname_dict):
                index=self.varname_dict[index]
            
            for gvar in self.var_lst:
                if (gvar.name==index):
                    return gvar
        else:
            return self.var_lst[index]
    
    def get_var(self, varname):
        """ check the var 

        Inputs:
        -------------------------------
        1. varname:<str>: var name 
        
        Returns:
        ----------------------------
        1. gvar:<io_var_lst>: variable for given name or index
        

        """
        return self.__getitem__(varname)
    
    

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
    
    def close_file(self):
        """
        close file if necessary
        """
        
        if (not self.read_to_array):
            if (oob.get_ot_type(self.nread)<>oob.ot_int):
                self.nread.close()


#===< TESTS >====    

if (__name__=='__main__'):
    
    print '>>>>1 test netcdf read' 
    datapath='/home/lfeng/local_disk/otool_data/ecmwf/'
    flnm=datapath+'cld_ecmwf.200301.nc'
    # setup io_var_cl
    
    
    vlon= ovio.io_var_cl('longitude','f',  ['longitude'], None)
    vlon.set_attr_rd_list(['units'])
    vlat= ovio.io_var_cl('latitude','f',  ['latititude'], None)
    vlat.set_attr_rd_list(['units'])
    
    vtime=ovio.io_var_cl('time', 'f', ['time'], None)
    vtime.set_attr_rd_list(['units'])
   
    
    
    # read infomation into these vars
    
    print '--- read from:', flnm 
    vlist=ncf_read_var(flnm, [vlon, vlat, vtime])

    units=vtime.get_attr('units')
       
    print '-----lon:', vlon[0:10]
    print '-----lat:',vlat[0:10]
    print '-----dimension name', vlon.dimname
    
    print '------shape of time', vtime.shape()
    

    print '>>>>2 copy var' 

    vlon_cp=vlon.copy()
    vlon_cp.set_name('lon')
    
    vlat_cp=vlat.copy()
    vlat_cp.set_name('lat')
    
    print '>>>>>3 test write'
    
    ncflnm='test_out.nc'
    
    
    ncf_save_var(ncflnm, [vlon_cp], \
                     [vlon, vlat, vtime], \
                     create_new=True)
    ncf_save_var(ncflnm, [vlat_cp])
    
    
    
    
    
    
    
    
    
        
    
