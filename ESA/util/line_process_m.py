"""
    Functions for processing text line
   
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.09.18
    
    
    1.Functions:
    ========================================================
    1. decode_line: split one line into terms
    
    2. decode_keyword_line: split one line into keyword and values 
    
    3. process_terms: select and convert items into required types. 
    
    4. read_table_from_lines: decode lines into a table (recarray)
    
    5. read_matrix_from_lines: decode lines into a matrix (array)
    
"""

import numpy as npy
import otool_obj as oob
import message_m as msm


def decode_line(line, delim):
        
    """ 
    splite one line into terms. 
    
    Inputs:
    -----------------------------------------------------------
    1. line:<str>: line to decode 
    2. delim:<str>: splitor for line 
        
    Returns:
    ---------------------------------------------------------
    out_terms:<list>: terms in the given line 

    
    """
    
    line=line.strip()

    line.replace('\n', '')
    
    # split the line 
    
    terms=line.split(delim)
        
    if (delim==' '):
        out_terms=list()
        for item in terms:
            item=item.strip()
            if (len(item)>0):
                # only non-empty item will be used
                out_terms.append(item)
        
                    
    else:
        # all item will be used 
        
        out_terms=list()
        
        for item in terms:
            item=item.strip()
            out_terms.append(item)
    
    return out_terms



def decode_keyword_line(line, delim_keyword, delim):
        
    """ 
    
    Decode one line into key word and values 
        
    Inputs:
    -----------------------------------------------------------
    1. line:<str>: line to decode 
    2. delim_keyword:<str>: the sperator between keyword and values
    3. delim:<str>: sperator for value line 
    
    Returns:
    ---------------------------------------------------------
    1. keyword:<str>: keyword (left term)
    2. terms:<list, t:str>: value (right term)
    
    Notes:
    
    1. if a line is given as line='A:1, 2,3'
    decode_keyword_line(line, ':', ',') will yield
    keyword='A'
    terms=['1', '2', '3']
    
    """

    
    
    line=line.strip()

    line.replace('\n', '')
    
    keyword, line_left=line.split(delim_keyword)
    terms=decode_line(line_left, delim)
    return keyword, terms


def process_terms_by_guess(terms):
    
    """
    convert items to float if possible
    
    Inputs:
    ----------------------------------------------------------
    1. terms:<list, t:str>: list of terms to be processes
    
    Returns:
    -------------------------------------------------------------
    1. new_terms_lst:<list>: terms in float or string 
    """
    new_terms_lst=list()
    for item in terms:
        try:
            val=float(item)
        except:
            val=item.strip()
        
        new_terms_lst.append(val)

    return new_terms_lst


        
    
    
def process_terms(terms,sel_col, col_type):
        
    """
    select terms, and convert them into certain types
    
    Inputs:
    ----------------------------------------------------------
    1. terms:<list, t:str>: list of terms to be processes
    2. sel_col:<list, t:integer>: terms to be chosen
    3. col_type:<list, t:str>: types to be converted to. 
    
    Returns:
    -------------------------------------------------------------
    1. sel_terms:<list>: selected terms 
    
    Notes:
    -------------------------------------------------------------
    1. If one element can be converted into  desired type, it will be replaced by a filling value. 
    
    
    """
    
    sel_id=0
    sel_terms=list()
    
    # S1 get selected item and convert it to required type
    # print terms[0], len(terms)
    
    
    for col_id in sel_col:
        
        item=terms[col_id]
        stype=col_type[sel_id]
        
        ##T2 convert the term to the required type
        
        if (stype in ['f', 'f4', 'f8', 'f16', 'f32', 'float']):
            ##  float
            try: 
                item=float(item)
            except:
                ## if not possible 
                item=oob.fill_val
            
        elif (stype in ['i', 'i4', 'i8', 'i16', 'i32', 'u8', 'int']):
                
            ## integer
            try: 
                item=int(item)
            except:
                ## if not possible
                item=oob.fill_val_int
        else:
            ## string 
            
            item=item.strip()
            
        sel_terms.append(item)
        sel_id=sel_id+1
    return sel_terms
    



def read_table_from_lines(lines, colnm_lst, colno_lst,\
                              coltype_lst, lstart, lend, \
                              delim=' ', comment_begin='#', \
                              **replaces):
    
    """ read txt lines into a record table
    
    Inputs:
    ------------------------------------------------------
    1. lines:<list, t:str>: lines from file 
    2. colnm_lst:<list, t:str_>:list of column name to be read
    3. colno_lst:<list, t:str_>:list of column number to be read
    4. coltype_lst:<list, t:str_>:list of column type to be converted to
    5. lstar:<int>: start row
    6. lend:<inte>: end row
    7. replaces: <dict>: the words need to be replace before conversion. 

    Returns:
    -------------------------------------------------
    1.dtable:<recarray>: data
        
        
    """
       
    data_lst=[]
    ncol=len(colno_lst)
    comm_len=len(comment_begin)
    
    
    for line in lines[lstart:lend]:
        line=line.strip()
        if (len(line)>0):
            for keyname in replaces:
                ## replace 
                line=line.replace(keyname, replaces[keyname])
            
            line_head=line[0:comm_len]
                
            if (line_head<>comment_begin):
                
                terms=decode_line(line, delim)
                terms=process_terms(terms, colno_lst, coltype_lst)
                
                terms=tuple(terms)
                        
                data_lst.append(terms)
        
        
    sname=tuple(colnm_lst)
    sform=tuple(coltype_lst)
        
    sdtype={'names':sname, 'formats':sform}
    dtable=npy.array(data_lst, dtype=sdtype)
    return dtable
        



def read_matrix_from_lines(lines, colno_lst,\
                               lstart, lend, \
                               delim=' ', comment_begin='#', \
                               **replaces):
    
    
    """ read txt lines into  a matrix 
        
    Inputs:
    ------------------------------------------------------
    1. lines:<list, t:str>: lines from file 
    2. colno_lst:<list, t:str_>:list of column number to be read
    4. lstar:<int>: start row
    5. lend:<inte>: end row
    6. replaces: <dict>: the words need to be replace before conversion. 

    Returns:
    -------------------------------------------------
    1.data:<array>: data matrix 
    
    Notes:
    ----------------------------------------------------
    1. Selected columns will be converted into float. 
    
    """
       
    data_lst=[]
    ncol=len(colno_lst)
    comm_len=len(comment_begin)
    ncol=len(colno_lst)
    
    coltype_lst=ncol*['f8']
    
    
    for line in lines[lstart:lend]:
        line=line.strip()
        if (len(line)>0):
            
            for keyname in replaces:
                ## replace 
                line=line.replace(keyname, replaces[keyname])
            
            line_head=line[0:comm_len]
                
            if (line_head<>comment_begin):
                
                terms=decode_line(line, delim)
                terms=process_terms(terms, colno_lst, coltype_lst)
                
                terms=tuple(terms)
                        
                data_lst.append(terms)
        
    data=npy.array(data_lst) 
    return data



