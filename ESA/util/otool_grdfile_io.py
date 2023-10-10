""" Class for gridded text files.  
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.12.04
    History: v0.95, 2013.02.19
    

Classes:
=============================================================
1. grdfile_desc_cl: class for ascii file access. 


"""

import numpy as npy
import otool_obj as oob
import otool_var_io as ovio
import gp_grid_m as grid_m
import line_process_m as lnpm 

import message_m as msm

#===< CLASSES >====      
class grdfile_desc_cl:
    
    """ class for description of ascii file  
    
    Members:
    =====================================================
    1. header_start:<int>: header starting row 
    2. data_start:<int>: data starting row 
    3. comment_begin:<str>: beginning character of comment lines
    4. mask_val:<float>: filling for bad or missing data
    5. colnm_lst:<list, t:str>: list of column names in the file
    6. data_lst:<list, t:array>: lits of raw data read from files 
    7. attr_dict:<dict>: attribute dictionarys
    8. grd_lst:<list, t:gp_grid_cl>: list of grids extracted from raw data 
    9. gdata:<list, t:array>:  list of gridded data. 
    10.flnm_lst:<list, t:str>: file names of the stored data set 
    11.ref_lst:<list, t:str/numeric>: values used to identify the data set
    12.cur_flnm:<str>: the name of the file to be read 
    13.cur_ref:<str>: current reference for file 
    ---> the reference can be time step etc
    
    
    

    Functions
    ==================================
    1.  __init__: initialization
    2. is_comment_line: check whether a line is for comments 
    3. find_colno: get column number 
    4. get_colno_list:  find column number for given column name list 
    5. read_matrix_from_lines
    6. convert_matrix_to_grid: convert matrix to grid data
    7. read_header_from_lines:  Read header from lines 
    
    # file 
    8. set_filename:  set the current filename and the current reference ID. 
    9. set_file_path: set the file path 
    10.  set_filename_format: set template for filename
    11. construct_filename: use template and keywords
    
    
    # member manager
    
    12.  get_index:   find the index in the file list for the given reference values
    13. append_data:   add the data set to data list
    14. del_data:  delete data set from the list 
    
    # attributes
    15. set_attr: set attributes
    16. get_attr: get attributes
    
    
    Attributes (reserved words)
    -----------------------------------------------------------
    1. flnm: <str>: full file name 
    2. filepath:  <str>: path of the file name 
    
    3. flnmformat:  <str>: template for filename
    
    """
    
    def __init__(self, filename,\
                     ref,\
                     header_start, \
                     data_start, \
                     comment_begin='#', \
                     delim=',', \
                     mask_val=oob.fill_val, \
                     colnm_lst=[], **keywords):
        """
        Inputs:
        ----------------------------------------------------------
        1. filename:<str>: name the file 
        2. ref:<str/numeric>: reference value 
        
        3. header_start:<integer>: header starting line Number 
        4. data_start:<integer>: header starting line Number
        5. comment_begin:<str>:  comment starting letter 
        6. delim:<str>: separator between items 
        7. mask_val: <float>: fillings for bad or missing values
        8. colnm_lst:<list, t:int>: list of columns to be read
        9. keywords:<dict>: attributes
        
        """
        # file structure
        self.header_start=header_start
        self.data_start=data_start
        self.comment_begin=comment_begin
        self.colnm_lst=colnm_lst
        self.delim=delim
        

        # data storage 
        
        self.cur_flnm=filename
        self.cur_ref=ref
        
        self.flnm_lst=[]
        self.data_lst=[]
        self.grd_lst=[] 
        self.gdata_lst=[]
        self.ref_lst=[]
        self.nset=0
        
        
        self.attr_dict={}
        
        
        if ('ot_type' in self.attr_dict):
            pass
        else:
            self.attr_dict.update({'ot_type':oob.ot_fdesc})

        self.mask_val=mask_val

        # update attributes
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
        
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

    def set_filename(self, flnm, ref):
        """
        set the current filename and the current reference  
        
        
        Inputs:
        ---------------------------------------
        1. yyyy:<int>: year
        2. mm:<int>: month
        3. dd:<int>: dd
        4. flnm:<str>: file name 
                
        """
        
        self.cur_flnm=flnm
        self.cur_ref=ref

    def get_index(self, ref):
        """
        find the index of the given reference values
        
        Inputs:
        ---------------------------------------------
        1. ref:<str/numeric>: reference value 
        
        Returns:
        -------------------------------------------
        1.index:<int>: if found, index of the reference value will be  returned 
        
        """
        
        for idx in range(self.nset):
            if (self.ref_lst[idx]==ref):
                return idx
        
        
        return self.mask_val
    
    def append_data(self, flnm, data, grd, gdata, ref):
        """
        
        add the data set to lists of the class
        
        
        Inputs:
        -------------------------------------------
        1. flnm: <str>:   file name 
        2. data: <array>:  raw data read from the file 
        3. grd: <ctm_grd>: grid extracted from raw data
        4. gdata:<array>:  data gridded over grd 
        5. ref:<str/numeric>: id (such as time etc) for the data set 
        
        
        """
        if (ref in self.ref_lst):
            for idx in self.ref_lst:
                if (ref==ref_lst[idx]):
                    self.flnm_lst[idx]=flnm
                    self.data_lst[idx]=data
                    self.grd_lst[idx]=grd
                    self.gdata_lst[idx]=gdata
                    self.ref_lst[idx]=ref
                    break
        else:
            self.flnm_lst.append(flnm)
            self.data_lst.append(data)
            self.grd_lst.append(grd)
            self.gdata_lst.append(gdata)
            self.ref_lst.append(ref)
            self.nset=self.nset+1
        
    def del_data(self, idx, by_ref=False):
        
        """
        Remove data from the list 
        
        Inputs:
        -----------------------------------------------------
        
        1.idx :<int/obj>: index of the data set to be deleted 
        2. by_ref:<T/F>: if ture , index will be treated as ref values for looking up ref_lst
       
        Returns:
        ---------------------------------------------------
        self.nset:<int>: the length of the current data lst 
        
        """
        
        if (by_ref):
            
            idx=self.find_index(idx)
            if (index==self.mask_val):
                msg='Not found'
                msm.show_err_msg(msg)
                return self.nset
        
            
        del self.flnm_lst[idx]
        
        del self.data_lst[idx]
        
        del self.grd_lst[idx]
        
        del self.ref_lst[idx]
        
        self.nset=self.nset-1
        
        return self.nset
    
    
    def read_header_from_lines(self, lines, head_st=None, head_end=None, marker=None):
        
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
        
        if (len(line)>comm_len):
            starter=line[0:comm_len]
            if (starter==self.comment_begin):
                return True
            
        return False
        
    

    def read_matrix_from_lines(self, lines, colno_lst=None,\
                                   lstart=None, lend=None, **replaces):
        
        """ read txt lines into data matrix 
        
        
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

        ## find lstart if necessary 
        
               
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
        
        data=lnpm.read_matrix_from_lines(lines, colno_lst,\
                                             lstart, lend, \
                                             delim=self.delim, \
                                             comment_begin=self.comment_begin, **replaces)

        
        
        return data
    
    
    def convert_matrix_to_grid(self, data, \
                                   colnm_lst, \
                                   colno_lst, \
                                   zname=None, \
                                   zval=None):
        
        """  wrapper for build_grid_from_matrix in gp_grid_m.py 

        Inputs:
        
        -----------------------------------------------------
        1. colnm_lst:<list, t:str>: (horizontal) axis names
        2. colno_lst:<list, t:integer>: columns of the (horizontal) axis

        3. data:<array, (nrow, ncol)>: data to be gridded' 
        --- if (zname and zval) are given,   each row in data is supposed to be in form of 
        --->:  [ax_1, ax_2, ax_3,.., ax_n, y1, y2, y3, ...]
        ---if (zname and zval) are not given, each row data is supposed to be in the form of 
        --- [ax_1, ax_2, ax_3, ...,ax_n, val]
        
        4. zname:<str>: name of  z (vertical) axis
        5. zval:<array>:  z (vertical) axis 
        

        Returns:
        ----------------------------
        1. grd:<gp_grid_cl>: grid 
        2. data:<array>: gridded data
        
        """
        
        grd, gdata=grid_m.build_grid_from_matrix(colnm_lst, colno_lst, \
                                                     data, \
                                                     zname=zname, \
                                                     zval=zval, \
                                                     mask_val=self.mask_val)
            
        return grd, gdata
    
        
    
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
            
            colno_lst.append(colno)
        
        return colno_lst
    
            
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

    
    def construct_filename(self, **keywords):
        """ construction file name use the keyword 

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
    
    print 'For test, see ESA/instrument/ak_file_io_m.py'
    

    
    
    
    
    
    
    
        
    
