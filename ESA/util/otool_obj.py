"""
  Module for constant and managing function for objects defined in toolkit
  
  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.06.30
  History: v0.95, 2013.01.06
  
  
  1. Parameters
  ===========================================
  # filling values
  1. fill_val_float=-999.0
  2. fill_val=-999.0
  3. fill_val_int=-999
  4. fill_val_str=""
  
  # python types
  
  5. ot_int=1
  6. ot_float=2
  7. ot_string=3
  8. ot_list=4
  9. ot_array=5
  10.ot_slice=6
  
  # ot types 
  11. ot_numeric=9
  12. ot_axis=100
  13. ot_grid=200
  14. ot_field=300
  15. ot_world=400
  16. ot_profile=500
  17. ot_vintpl=600
  18. ot_hinterp=700
  19. ot_slice=800
  20. ot_iovar=900
  21. ot_fdesc=1010
  22. ot_ncfdesc=1010
  23. ot_ak=1030
  24. ot_err=1040
  25. ot_cloud=1050
  26. ot_aod=1060
  
  27. ot_no_ot_type=fill_val_int
  28. ot_bf=1070 
  29. ot_obs_op=1080
  30. ot_stv=1090
  31. ot_menu=1100
  32. ot_def=1200
  33. ot_lc=1300

  
  2. Functions:
  =======================================================
  1.get_ot_type: get type of the obj
  2.check_type : check whether the obj has the type as required 
  3.compare_type: check whether two objects have the same type  
  4.compare_attr: compare attributes of the ot object 
  5.is_attr_matched: compare 
  6. check_ot_shape: check whether obj has the required shape. 
  7. translate_dict_ID: find ID  for given names from a ID dictionary
  
"""


import numpy as npy
#<<< Paramters >>> 

# filling values 
fill_val_float=-999.0
fill_val=-999.0
fill_val_int=-999
fill_val_str=""

# basic type of variables


ot_int=1
ot_float=2
ot_string=3
# short name

ot_str=3

ot_list=4
ot_array=5
ot_slice=6

# type numeric 

ot_numeric=9

# type id (class) defined in tool kit 

ot_axis=100
ot_grid=200
ot_field=300
ot_world=400
ot_profile=500
ot_slice=600

ot_vintpl=700
ot_hinterp=800

ot_iovar=900
ot_fdesc=1010
ot_ncfdesc=1020
ot_ak=1030
ot_err=1040
ot_cloud=1050
ot_aod=1060
ot_bf=1070
ot_obs_op=1080
ot_stv=1090
ot_menu=1100
ot_def=1200
ot_lc=1300



ot_no_ot_type=-1000

#<<< FUNCTIONS >>> 

def check_ot_shape(obj, dims):
    """  
    check_ot_shape: check whether obj has the required shape.

    Inputs:
    ---------------------------------------
    1. obj:<obj>: object which shape will be test. 
    2. dims:<array/list>: the required shape 
    
    Returns:
    -----------------------------------------
    1. matched:<T/F>: True, if the obj shape is matched with the requirements. 
    
    """
    if (obj==None):
        return False
    
    try:
        obs_dims=npy.shape(obj)
    except:
        return False
    
    ndim=len(obs_dims)
    
    if (ndim<len(dims)):
        return False
    ix=0
    
    for nx in dims:
        if (obs_dims[ix]<>nx):
            return False
        ix=ix+1
        
    return True


def get_ot_type(obj):
    
    """ check the object has type as required
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    
    Returns:
    -------------------------------
    1. obj_type:<integer>: type of object tool  
    
    """

    if (obj==None):
        return ot_no_ot_type
    
    
    try:
        # ot object 
        obj_type=obj.get_attr('ot_type')
        
        return obj_type
    
    except AttributeError: 
        # python object 
        
        if (type(obj)==type(1)): # integer 
            return ot_int
        
        elif (type(obj)==type(1.0)): # float
            return ot_float
        
        elif (type(obj)==type("abcd")): # string 
            return ot_string
        
        elif (type(obj)==type([1.0])):  #list
            return ot_list
        
        elif (type(obj)==type(slice(1))): # slice 
            return ot_slice
        
        elif (type(obj)==type(npy.array([0,1,2,3]))): # array
            return ot_array
        
        else:
            return ot_no_ot_type



def check_type(obj, ot_type):
    """ check whether the object has the type as required
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    2. ot_type:<integer>: the expected type of otool object. 
    
    Returns:
    -------------------------------
    1. matched:<T/F>: True, if matched, 
    
    Notes:
    --------------------------------------------------
    1.  when obj is  list or array, we can not check its each element.
    



    """

    
    if (obj==None):
        return False
    
    obj_type=get_ot_type(obj)
    
    
    if (ot_type==ot_numeric): # whether it is a numeric  (see notes 1.)
        

        if (obj_type in [ot_int, ot_float, ot_array]):
            return True
        else:
            return False
    
    if (obj_type==ot_type): # if matched 
        return True
    else:   # if not 
        return False
    

        
    return obj_type




def compare_type(obj, obj2):
    
    """ check the object has type as required
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    2. obj2:<obj>: otool class for comparison
        
    Returns:
    -------------------------------
    1. is_match:<T/F>: True, if matched, 
    
    """
    
    if ((obj==None) or (obj2==None)):
        return False
    
    # check type for obj

    obj_type=get_ot_type(obj)
    if (obj_type==ot_no_ot_type):
        return False
    
    # check type for obj2 

    obj_type2=get_ot_type(obj)
    if (obj_type2==ot_no_ot_type):
        return False
    
    if (obj_type==obj_typ2): # if matched 
        return True
    else:
        return False
    

def is_attr_matched(obj, **keywords):
    
    """ check objecs have the set of attributes set out by keyword_lst
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    2. obj2:<obj>: otool class for comparison
    3. keywords:<dict>: key attributes to be comparied 
    
    Returns:
    -------------------------------
    1. matched:<T/F>: True, if matched. 

    Notes: 
    1. if no criteria (keywords) is given, a True will be return 
    
    
    """
    matched=True
    
    for keyname in keywords:
        keyval=keywords[keyname]
        if (keyval==None):
            pass
        else:
            try:
                att_val=obj.get_attr(keyname)
            except AttributeError: 
                #  print 'error match', keyname
                return False
            
            if (keyval<>att_val):
                # print 'not match', keyval
                return False
    
    return matched


def compare_attr(obj, obj2, attr_lst):
    
    """ check objecs have the set of attributes set out by keyword_lst
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    2. obj2:<obj>: otool class for comparison
    3. attr_lst:<dict>: key attributes to be comparied 
    
    Returns:
    -------------------------------
    1. matched:<T/F>: True, if matched, 
    
    """
    matched=True
    
    for keyname in attr_lst:
        try:
            att_val=obj.get_attr(keyname)
        except AttributeError: 
            return False
        
        try:
            att_val2=obj2.get_attr(keyname)
        except AttributeError: 
            return False
        
        if (att_val<>att_val2):
            return False
    
    return matched


    


def translate_dict_ID(id_dict, names):
    
    """ find ID  for names from a ID dictionary 
    
    Inputs:
    -----------------------------------------
    1.id_dict:<dict>: dictionary in the form of {keyname:id} for ID of keynames (see Note 1)
    2.names:<array/obj>: names to be checked against id_dict
    

    Returns:
    ------------------------------------------
    vid_lst:<array/obj>: id of the given names
    
    
    Notes:
    ---------------------------------------------
    1. ID could be a unqie number or string 
    
    """

    
    if (get_ot_type(names)==ot_str):
        # if names is a single value 
        vid=id_dict[names]
        return vid
    else:
        # if names is a list or array
        
        names=npy.array(names)
        nob=npy.size(names)
        
        
        ##  construct vid_lst array
        keynames=id_dict.keys()
        vid=id_dict[keynames[0]]
        
        vid_lst=nob*[vid]
        vid_lst=npy.array(vid_lst)
        
        set_names=set(names)
        set_names=list(set_names)
        for sname in set_names:
            ## fill iviewmode by each corresponding values
            idx=npy.where(names==sname)
            vid=id_dict[sname]
            vid_lst[idx]=vid
            
    
    return vid_lst

        
    

