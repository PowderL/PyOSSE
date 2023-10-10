# filling values 
fill_val_float=-999.0
fill_val_int=999

# type (class) defined in tool kit 
ot_axis=100
ot_grd=200
ot_field=300
ot_world=400
ot_no_ot_type=-1000


def check_type(obj, ot_type):
    """ check the object has type as required
    Inputs:
    --------------------------------
    1. obj:<obj>: otool class
    2. ot_type:<integer>: the expected type of otool class

    Returns:
    -------------------------------
    1. is_match:<T/F>: True, if matched, 

    """

    if (obj==None):
        return False
    
    try:
        obj_type=obj.get_attr('ot_type')
    Except AttributeError: 
        obj_type=ot_no_ot_type
        
    if (obj_type==ot_grd):
        return True
    else:
        return False
    


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

    if (obj==None):
        return False
    
    try:
        obj_type=obj.get_attr('ot_type')
    Except AttributeError: 
        obj_type=ot_no_ot_type
        
    if (obj_type==ot_grd):
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
    1. matched:<T/F>: True, if matched, 
    
    """
    matched=True
    
    for keyname in keywords:
        keyval=keywords[keyname]
        try:
            att_val=obj.get_attr(keyname)
        Except AttributeError: 
            return False
        
        if (keyval<>att_val):
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
        Except AttributeError: 
            return False
        
        
        try:
            att_val2=obj2.get_attr(keyname)
        Except AttributeError: 
            return False
        
        if (att_val<>att_val2):
            return False
    
    return matched


    


    

        
    

