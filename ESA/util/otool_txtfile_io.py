""" Class for text file IO. It contains codes for reading ascii table file 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.11.28
    History: v0.95, 2013.01.22
    

    Classes:
    =============================================================
    1. file_desc_cl: class for ascii file access. 
    

"""

import numpy as npy
import otool_obj as oob
import otool_var_io as ovio
import gp_grid_m as grid_m

import message_m as msm
import line_process_m as lnpm 
#===== < CLASSES >======

class file_desc_cl:
    
    """ class for description of ascii file  
    Members:
    =====================================================
    1. header_start:<int>: header starting row 
    2. data_start:<int>: data starting row 
    3. comment_begin:<str>: beginning character of comment lines
    4. mask_val:<float>: filling for bad or missing data
    5. colnm_lst:<list, t:str>: list of column names in the file
    6. dtable:<recarray>: recarray read from files 
    7. data:<array>: array read from files 
    8. attr_dict:<dict>: attribute dictionarys
    9. grd:<gp_grid_cl>: grid of the data 
    10.gdata:<array>: gridded data 
    
        
    
    Functions
    ==================================
    1.  __init__: initialization
    2. is_comment_line: check whether a line is for comments  
    3. read_table_from_file: read columns of a text file  into a table (recarray)
    4. read_table_from_file_byvar: read columns named by a var_lst 
    ---from a text file  into a table (recarray)
    5. read_table_from_lines_byvar: decode lines into a table 
    ---(recarray) defined by a var_lst
    
    6. read_table_from_lines: decode lines into a table (recarray) (wrapper for functions in line_process_m)
    
    7. read_matrix_from_lines: decode lines into matrix (array). (wrapper for functions in line_process_m)
     
    
    8. find_colno: get column number 
    9. get_colno_list: find the column numbers for a list of names
    10. get_colno_list_byvar: find the column numbers for a list of var 
    11. get_colnm_list_byvar: find the column names for a list of var 
    12. get_coltype_list_byvar: find the types for a list of var 
    
    13. get_lon: get longitude if existing 
    14. get_lat: get latitude if existing 
    15. get_time: get time if existing 
    16. get_column_data: get data for a given column 
    17. convert_matrix_to_grd: convert matrix to grid data
    
    18. read_header_from_lines:read header from given lines
    19. locate_section: find the beginning (or end) of section 
    


    # file 
    20. set_filename:  set the current filename and the current reference ID. 
    21. set_file_path: set the file path 
    22.  set_filename_format: set template for filename
    23. construct_filename: use template and keywords
    
    
    # attributes
    24. set_attr: set attributes
    25. get_attr: get attributes
    
    Attributes (reserved words)
    -----------------------------------------------------------
    1. flnm: <str>: full file name 
    2. filepath:  <str>: path of the file name 
    3. flnmformat:  <str>: template for filename
    
    """
    
    def __init__(self, header_start, \
                     data_start, \
                     comment_begin='#', \
                     delim=',', \
                     mask_val=oob.fill_val, \
                     colnm_lst=[], **keywords):
        """
        Inputs:
        ----------------------------------------------------------
        1. header_start:<integer>: header starting line Number 
        2. data_start:<integer>: header starting line Number
        3. comment_begin:<str>:  comment starting letter 
        4. delim:<str>: separator between items 
        5. mask_val: <float>: fillings for bad or missing values
        6. colnm_lst:<list, t:int>: list of columns to be read
        7. keywords:<dict>: attributes
        
        """
        
        self.header_start=header_start
        self.data_start=data_start
        self.comment_begin=comment_begin
        self.mask_val=mask_val
        self.colnm_lst=colnm_lst
        self.dtable=None
        self.data=None
        self.grd=None 
        self.gdata=None
        
        
        self.delim=delim
        
        self.attr_dict={}
        
        
        if ('ot_type' in self.attr_dict):
            pass
        else:
            self.attr_dict.update({'ot_type':oob.ot_fdesc})

        # update attributes
            
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
        
    
    def is_comment_line(self, line):
   
        """ check a line is just for comments
        
        Inputs:
        ----------------------------------
        1. line: <str>: line 
        
        Returns:
        -------------------------------------
        1. is_comment:<T/F>: True if the line start with comment_begin
        

        """
        
        line=line.strip()
        comm_len=len(self.comment_begin)
        line_head=line[0:comm_len]
        
        if (line_head==self.comment_begin):
            return True
        else:
            return False


    def read_table_from_file_byvar(self, flnm, var_lst):
       
        """ read txt table files into var and return a record table
        
        Inputs:
        -----------------------------------------
        1. flnm:<str>: file name 
        2. var_lst:<list, t:io_var_cl_>:list of var to be read
        

        Returns:
        ----------------------------
        1. dtable:<recarray>: data 
        
        """

        colno_lst=list()
        colnm_lst=list()
        nvar=len(var_lst)
        # collect the column name and column no needs to be read
        
        for ivar in range(nvar):
            gvar=var_lst[ivar]
            colnm=gvar.name
            colnm_lst.append(colnm)
            # if column number is not given, it will be searched from colnm_lisr 
            
            colno=gvar.colno
            
            if (colno==int(self.mask_val)):
                colno=self.find_colno(colnm)
                colno_lst.append(colno)
        
        
        
        self.dtable=npy.genfromtxt(flnm, comments=self.comment_begin, \
                                       skip_header=self.data_start, \
                                       usecols=colno_lst, names=colnm_lst)
        
        
        return self.dtable
    
    
    
    
    def locate_section(self, lines, marker):
    
        """find the first line starting with the given marker 
        Inputs:
        ----------------------
        1. lines:<list, t:str>: list of lines
        2. marker:<str>: marker for the line 

        Returns:
        -----------------------------
        1. pos:<int>: the location  of the line starting with marker 


        """
        ln_marker=len(marker)
        pos=0
        for line in lines:
            line=line.strip()
            nlen=len(line)
            if (nlen>ln_marker):
                starter=line[0:ln_marker]
                if (starter==marker):
                    return pos
            pos=pos+1
        return int(self.mask_val)
    
    def read_table_from_file(self, flnm, colnm_lst,  colno_lst):
        """ read txt table files into var and return a record table
        
        Inputs:
        -----------------------------------------
        1. flnm:<str>: file name 
        2. var_lst:<list, t:io_var_cl_>:list of var to be read
        
        
        Returns:
        ----------------------------
        1. dtable:<recarray>: data 
        
        """

        
        self.dtable=npy.genfromtxt(flnm, comments=self.comment_begin, \
                                       skip_header=self.data_start, \
                                       usecols=colno_lst, names=colnm_lst)
        
        
        return self.dtable

    def read_header_from_lines(self, lines, \
                                   head_st=None, \
                                   head_end=None, marker=None):
        
        """
        Read header from lines 
        Inputs:
        -----------------------------------
        1. lines:<list, t:str>: list of lines to be read
        2. head_st:<int/str>: start of the head section 
        3. head_end:<int/str>: end of the head section 
        4. marker:<str>: marker for header lines 
        
        Returns:
        -----------------------------------------
        1. head_line:<list, t:str>: head lines (see Note 1)

        Notes:
        -----------------------------------------
        1. if marker==None,  the section between header_start and data_start 
        will be treated as header section 
        
        """
        
        if (head_st==None):
            head_st=self.header_start
        
        if (oob.get_ot_type(head_st)==oob.ot_str):
            head_st=self.locate_section(lines, head_st)
            
        
        
        if (head_end==None):
            head_end=self.data_start
            if (oob.get_ot_type(head_end)==oob.ot_str):
                head_end=self.locate_section(lines, head_end)
                
        else:
            
            if (oob.get_ot_type(head_end)==oob.ot_str):
                head_end=self.locate_section(lines, head_end)
                
                ## added one to include these
                
                head_end=head_end+1
                
        
        if (marker==None):
            head_line=lines[head_st:head_end]
        else:
            ln_mark=len(marker)
            head_line=[]
            for line in lines[head_st:head_end+1]:
                line=line.strip()
                if (len(line)>ln_mark):
                    starter=line[0:ln_mark]
                    if (starter==marker):
                        head_line.append(line)
            
                
        return head_line
    
    

    
    
    
    def read_table_from_lines_byvar(self, lines, var_lst,\
                                        lstart=None, lend=None, \
                                        **replaces):
        """ read txt table files into var and return a record table
        
        Inputs:
        -----------------------------------------
        1. lines:<list, t:str>: lines from file 
        2. var_lst:<list, t:io_var_cl_>:list of var to be read
        3. lstar:<int>: start row
        4. lend:<inte>: end row
        5. replaces: <dict>: the words need to be replace before conversion. 
        
        Returns:
        --------------------------------
        1.dtable:<recarray>: data
        
        
        """
        
        colno_lst=self.get_colno_list_byvar(var_lst)
        coltype_lst=self.get_coltype_list_byvar(var_lst)
        colnm_lst=self.get_colnm_list_byvar(var_lst)
        
        
        if (lstart==None):
            lstart=self.data_start
        if (lend==None):
            lend=len(lines)

        if (oob.get_ot_type(lstart)==oob.ot_str):
            lstart=self.locate_section(lines, lstart)
        
            
        self.dtable=lnpm.read_table_from_lines(lines, colnm_lst, colno_lst,\
                                                   coltype_lst, lstart, lend, \
                                                   delim=self.delim, \
                                                   comment_begin=self.comment_begin, **replaces)
        
        return self.dtable
    
    def read_table_from_lines(self, lines, colnm_lst, colno_lst,\
                                  coltype_lst=None,\
                                  lstart=None, lend=None, \
                                  **replaces):
        
        """ read txt lines into  a record table
        
        Inputs:
        -----------------------------------------
        1. lines:<list, t:str>: lines from file 
        2. colnm_lst:<list, t:str_>:list of column name to be read
        3. colno_lst:<list, t:str_>:list of column number to be read
        4. coltype_lst:<list, t:str_>:list of column type to be read
        5. lstar:<int>: start row
        6. lend:<inte>: end row
        7. replaces: <dict>: the words need to be replace before conversion. 

        Returns:
        --------------------------------
        1.dtable:<recarray>: data
        
        Notes:
        
        --------------------------------------------
        1. If one value is not able to converted the required type, it will be replaced by a filling value.  
        2. If coltype is not provided, all columns will be converted to float. 
        
        """
       
        
        if (lstart==None):
            lstart=self.data_start
        if (lend==None):
            lend=len(lines)
            
        if (oob.get_ot_type(lstart)==oob.ot_str):
            lstart=self.locate_section(lines, lstart)
            
        ncol=len(colnm_lst)
        if (coltype_lst==None):
            # set to float by default 
            coltype_lst=ncol*['f8']
        
        self.dtable=lnpm.read_table_from_lines(lines, colnm_lst, colno_lst,\
                                                   coltype_lst, lstart, lend, \
                                                   delim=self.delim, \
                                                   comment_begin=self.comment_begin, \
                                                   **replaces)
        
        
        return self.dtable


    def read_matrix_from_lines(self, lines, colno_lst=None,\
                                   lstart=None, lend=None, **replaces):
        
        """ read txt lines into  a record table
        
        Inputs:
        -----------------------------------------
        1. lines:<list, t:str>: lines from file 
        2. colno_lst:<list, t:str_>:list of column name to be read
        3. lstar:<int>: start row
        4. lend:<inte>: end row
        5. replaces: <dict>: the words need to be replace before conversion. 
        
        Returns:
        --------------------------------
        1.data:<array>: data matrix 
        
        """
        
        
        if (lstart==None):
            lstart=self.data_start
        if (lend==None):
            lend=len(lines)

        ## find colno_lst if necessary
               
        if (oob.get_ot_type(lstart)==oob.ot_str):
            lstart=self.locate_section(lines, lstart)
        
            
        if (colno_lst==None):
            
            for line in lines[lstart:lend]:
                line=line.strip()
                
                if (len(line)>0):
                    for keyname in replaces:
                    ## replace 
                        line=line.replace(keyname, replaces[keyname])
                        
                if (not self.is_comment_line(line)):
                ## decode line into terms
                    terms=lnpm.decode_line(line, self.delim)
                    ncol=len(terms)
                    colno_lst=range(ncol)
                    break
        
        self.data=lnpm.read_matrix_from_lines(lines, colno_lst,\
                                                  lstart, lend, \
                                                  delim=self.delim, \
                                                  comment_begin=self.comment_begin, **replaces)
        

        return self.data
    
    
    def convert_matrix_to_grid(self, colnm_lst, \
                                   colno_lst, zname=None, zval=None):
        """  wrapper for build_grid_from_matrix in gp_grid_m.py 

        Inputs:
        -----------------------------------------------------
        1. colnm_lst:<list, t:str>: axis names
        2. colno_lst:<list, t:integer>: columns of the axis

        3. data:<array, (nrow, ncol)>: data to be gridded' 
        --- if (yname and yval) are given,   data is supposed to be in the column of 
        --- [ax_1, ax_2, ax_3,.., ax_n, y1, y2, y3, ...]
        ---if (yname and yval) are not given, data is supposed to be in the form of 
        --- [ax_1, ax_2, ax_3, ...,ax_n, val]
        
        4. zname:<str>: name of  z (vertical) axis
        5. zval:<array>:  z (vertical) axis 
    

        Returns:
        ----------------------------
        1. grd:<gp_grid_cl>: grid 
        2. data:<array>: gridded data
        
        """
        
        grd, gdata=grid_m.build_grid_from_matrix(colnm_lst, colno_lst, \
                                                     self.data, \
                                                     zname=zname, \
                                                     zval=zval, \
                                                     mask_val=self.mask_val)
        self.grd=grd
        self.gdata=gdata
        
        return self.grd, self.gdata
    
        
    
    def find_colno(self, colnm):
        
        """ find the column no for a given column name 
        Inputs:
        -----------------------------------------
        1. colnm:<str>: column name 
        
        
        Returns:
        ---------------------------------------
        1. col_idx:<integer>: the column with name given by colnm
        
        """
        
        col_idx=int(self.mask_val)
        
        idx=0
        
        for rcolnm in self.colnm_lst:
            if (rcolnm==colnm):
                col_idx=idx
                return col_idx
            idx=idx+1
        
        return col_idx
    
    def get_colno_list_byvar(self, var_lst):
    
        """ find column number for a list of ot_io_var
        
        Inputs:
        ----------------------------------------------------
        1. var_lst:<list, t:io_var_cl_>:list of var to be read
        
        
        Returns:
        -----------------------------------------------------
        1. colno_lst:<list, t:integer>: column number of each var 
        
        """
        
        colno_lst=list()
        nvar=len(var_lst)
        # collect the column name and column no needs to be read
        
        for ivar in range(nvar):
            gvar=var_lst[ivar]
            colnm=gvar.name
            colno=self.find_colno(colnm)
            if (colno<>int(self.mask_val)):
                gvar.set_colno(colno)
            
            colno_lst.append(gvar.colno)
    
        return colno_lst
    
    
    def get_coltype_list_byvar(self, var_lst):
    
        """ find column type for a list of ot_io_var
        
        Inputs:
        ----------------------------------------------------
        1. var_lst:<list, t:io_var_cl_>:list of var to be read
        
        
        Returns:
        -----------------------------------------------------
        1. coltype_lst:<list, t:str>: type for each var 
        
        """
        
        coltype_lst=list()
        nvar=len(var_lst)
        # collect the column name and column no needs to be read
        
        for ivar in range(nvar):
            gvar=var_lst[ivar]
            coltype=gvar.type
            coltype_lst.append(coltype)
            
        return coltype_lst
    

    def get_colnm_list_byvar(self, var_lst):
    
        """ find column names for a list of ot_io_var
        
        Inputs:
        ----------------------------------------------------
        1. var_lst:<list, t:io_var_cl_>:list of var to be read
        
        
        Returns:
        -----------------------------------------------------
        1. colnm_lst:<list, t:str>: names for each var 
        
        """
        
        
        colnm_lst=list()
        nvar=len(var_lst)
        # collect the column name and column no needs to be read
        
        for ivar in range(nvar):
            gvar=var_lst[ivar]
            colnm=gvar.name
            colnm_lst.append(colnm)
        
        return colnm_lst
    

    
    
    def get_colno_list(self, colnm_lst):
    
        """ find column number for given name 

        
        Inputs:
        ----------------------------------------------------
        1. colnm_lst:<list, t:str>:list of column name  to be read
        
        
        Returns:
        -----------------------------------------------------
        1. colno_lst:<list, t:integer>: column no for each column name
        
        
        """
        
        colno_lst=list()
        nvar=len(colnm_lst)
        # collect the column name and column no needs to be read
        
        for ivar in range(nvar):
            colnm=colnm_lst[ivar]
            colno=self.find_colno(colnm)
            print colnm, colno
            
            colno_lst.append(colno)
        
        return colno_lst


    
            
    def get_lon(self, colnm='lon'):
        """get longitude 
        """
        return self.dtable[colnm]
    def get_lat(self, colnm='lat'):
        """get latitude
        """
        return self.dtable[colnm]

    def get_time(self, colnm='time'):
        """get time 
        """
        return self.dtable[colnm]
        
    def get_column_data(self, colnm):
        """get value
        """
        return self.dtable[colnm]

    def set_file_path(self, path):
        """ set path for data file 
        """
        self.set_attr('filepath', path)
    
    def set_filename_format(self, flnm_format):
        """set template for filename """

        self.set_attr('flnmformat', flnm_format)

    
    def construct_filename(self, **keywords):
        """ construction file name use the keyword """
        
        sflnm=self.get_attr('flnmformat')
        if ('filepath' in self.attr_dict):
            filepath=self.get_attr('filepath')
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





#===< TESTS >====          

if (__name__=='__main__'):
    
    
    datapath='/home/lfeng/local_disk/otool_data/oco/aqua_0.25x0.25/'
    flnm=datapath+'aqua_num_nadir_205.dat'
    
    colnm_lst=['lat', 'lon', 'number', 'date', 'sza',  'effective']
    
    header_start=3
    data_start=4

    
    fs=file_desc_cl( header_start, data_start, delim=' ', colnm_lst=colnm_lst)
    
    vlon= ovio.io_var_cl('lon','f',  ['longitude'], None)
    vlon.set_attr_rd_list(['units'])
    vlat= ovio.io_var_cl('lat','f',  ['latititude'], None)
    vlat.set_attr_rd_list(['units'])
    
    vtime=ovio.io_var_cl('time', 'f', ['time'], None)
    vtime.set_attr_rd_list(['units'])

    vdate= ovio.io_var_cl('date','f',  ['date'], None)
    vsza=ovio.io_var_cl('sza','f',  ['date'], None)
    vnum=ovio.io_var_cl('number','f',  ['date'], None)
    
    tables=fs.read_table_from_file_byvar(flnm, [vlon, vlat, vdate, vsza, vnum])
    print tables['lon']
    print tables['lat']
    print tables.dtype.names
    
    # resize 
    numb=tables['number']
    usd_idx=npy.where(numb>0)
    du=tables[usd_idx]
    print 'new-shape', npy.shape(du)
    print 'old-shape', npy.shape(tables)
    del tables
    print 'read from file'

    fl=open(flnm,'r')
    lines=fl.readlines()
    
    fl.close()
    
    table=fs.read_table_from_lines_byvar(lines, [vlon, vlat, vdate, vsza, vnum])
    print len(table)
    print table.dtype.names
    

    data=fs.read_matrix_from_lines(lines, colno_lst=[0, 2,3])
    print npy.shape(data)
    

    
    
    
    
    
    
    
        
    
