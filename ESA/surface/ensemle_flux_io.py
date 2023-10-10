"""
 Asistent functions for generating ensemble fluxes 
 

"""

def define_ens_step(yyyy_st, mm_st, dd_st, time_step, nstep, **keywords):
    
    """
    define  time step for ensemble fluxes 
    
    Inputs:
    --------------------------------------------------
    
    1. yyyy_st:<int>: start year 
    2. mm_st:<int>: start month 
    3. dd_st:<int>: start day 
    4. time_step:<int>: time interval for fluxes 
    5. nstep:<int>: numbe of step
    6. keywords:<dict>: dictionary for extra inputs
    ---Reserved entries:
    -->reset_at_new_year: <int>: 1=move doy to doy=1 for each new year 
    -->reset_at_new_month: <int>: 1=move doy to the fist day of each new month 
    
    
    Outputs:
    --------------------------------------------------
    1.ens_yyyy_lst:<list, t:int>: starting year
    2.ens_doy_lst:<list, t:int>: starting doy 
    
    """

    reset_at_new_year=1
    reset_at_new_month=1

    if ('reset_at_new_year' in keywords):
        'reset_at_new_year'=keywords['reset_at_new_year']
    
    
    if ('reset_at_new_month' in keywords):
        'reset_at_new_month'=keywords['reset_at_new_month']
    
    
    
    ens_doy_list=list()
    ens_yyyy_list=list()
    
    # S1: date for first step 

    
    cur_dd=dd_st
    cur_mm=mm_st
    cur_yyyy=yyyy_st
    
    cur_doy=tm.day_of_year(cur_yyyy, cur_mm, cur_dd)
    ens_doy_list.append(cur_doy)
    ens_yyyy_list.append(cur_yyyy)
    ndays_year=tm.days_in_year(cur_yyyy)
    
    # S2: loop over nstep , and end with nstep+1
    
    
    for istep  in range(nstep+1):
        
        new_doy=cur_doy+time_step
        # #c: check year begin
        
        if (new_doy>ndays_year):
            cur_yyyy=cur_yyyy+1
            if (reset_at_new_year):
                new_doy=1
            else:
                new_doy=new_doy-ndays_year
            
            ndays_year=tm.days_in_year(cur_yyyy)

        new_dd, new_mm, new_yyyy=tm.doy_to_time_array(new_doy, cur_yyyy)
        
        # #c: check month begin 
        
        if (new_mm<>cur_mm):
            cur_mm=new_mm
            if (ens_reset_at_new_month):
                new_doy=tm.day_of_year(cur_yyyy, cur_mm, 1)
        cur_doy=new_doy
        
        ens_doy_lst.append(new_doy)
        ens_yyyy_lst.append(cur_yyyy)
    
    
        return ens_doy_lst, ens_yyyy_lst
    
                                   

def create_ens_dim_lst(dim_name_lst):
    """
    construct io avr list for dimensions used in ensemble fluxes 
    
    Inputs:
    ------------------------------------------
    1. dim_name_lst:<str>: name of dimensions used in the file 
    
    Outputs:
    ------------------------------------------
    
    """
    var_lst=[]
    for dname in dim_name_lst:
        var= ovio.io_var_cl(dname,'f',  \
                                [dname], \
                                None)
        
        var_lst=var_lst+[var]

    return var_lst



def create_ens_var_lst(varname_lst, dim_name_lst):
    
    """
    construct io avr list for variables to be stored n ensemble fluxes 
    
    Inputs:
    ------------------------------------------
    1. dim_name_lst:<str>: name of dimensions used in the file 
    
    Outputs:
    ------------------------------------------
    
    """
    var_lst=[]
    nvar=len(varname_lst)
    for ivar in range(nvar):
        for dname in varname:
            var= ovio.io_var_cl(dname,'f',  \
                                    [dname], \
                                None)
        
        var_lst=var_lst+[var]

    return var_lst

