"""
read configuration file for pre-define flux perturbations 

"""

import numpy as npy 

import ESA.util.line_process_m as lpm

def read_em_conf(flnm):
    """
    Read emission configurations on flux perturbations (Basis functions)  ensemble
    
    Inputs:
    ---------------------------------------------------------------
    
    1. flnm:<str>: name of the text file on flux perturbation ensemble 

    Outputs:
    ----------------------------------------------------------------------
    1. doy:<list, t:int>:  starting doy of the pulse-like perturbations
    2. year:<list, t:int>: starting year of the pulse-like perturbations
    3. npb:<list, t:int>:  number of perturbation functions in each file
    4. emsflnm:<list, t:str>: list of files storing the flux ensemble 

    """
    
    # S1: open config file
    
    try:
        
        fl_cfg=open(flnm, 'r')
    except IOError:
        #br2 if file not exists, exit  
        print 'Can not open '+flnm.strip()
        return None, None
    
    # S2: read in lines
    
    lines=fl_cfg.readlines()
    fl_cfg.close()

    # S3: scan the lines and fill the outputs
    

    doy, year, npb, emsflnm=list(), list(), list(), list()
    
    # #: loop over lines 

    for  line in lines[2:]:
        line=line.strip()
        line.replace('\n', '')
        if (len(line)>0):
            terms=line.split(',')
            year.append(int(terms[0]))
            doy.append(int(terms[1]))
            npb.append(int(terms[2]))
            emsflnm.append(terms[3].strip())

    # loop end lines

    
            
    return doy, year, npb, emsflnm

            
                        
    
    
    
    
