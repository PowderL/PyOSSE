""" module for time conversions.  

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.04.01
    History: v0.95, 2012.08.21
    
    Functions:
    =====================================================
    1. days_in_year: number of days in each year
    2. day_of_year:  day of year for calender date
    3. get_doy: a wrapper for day_of_year
    4. tai85_to_utc:  convert a time in tai85 format to UTC string.
    5. doy_to_utc:    convert day of year to utc string. 
    6. doy_to_time_array: convert day of year to calendar date (yyyy, mm, dd)
    7. utc_to_time_array(utc): convert the utc string to the list of (yyyy, mm, dd, hh, mi, sec)
    8. time_array_to_utc: convert yyyy, mm, dd, hh, mi, sec to utc string
    9. convert yyyy, mm, dd, hh, mi, sec to tau (i.e, hours since 1985.01.01, 00:00)
    10.tau_to_time_array(tau): convert tau (i.e, hours since 1985.01.01, 00:00) to  yyyy, mm, dd, hh, mi, sec 
    11. next_doy: get day of year for current time + n days 
    12. get_ut_loc_time_slot: get indexes for UT time series which fall within the local time slot
    
"""

import numpy as npy
import pylab as pylab


import datetime
# second_leap=33 # there is a 32 seconds in 2000
second_leap=0
TAI0=datetime.datetime(1985, 1,1,0,0,0)
def_yyyy=2004


def days_in_year(yyyy):
   
    """ get the days in a given year for calender date
    
    Input: 
    -------------------------------------------
    1. yyyy, mm, dd:<integer>:year, month and day 
    
    Returns:
    ---------------------------------------------
    1.  days in year
    
    """
    days=365
    if (npy.mod(yyyy, 4)==0):
        days=366
    return days


def day_of_year(yyyy, mm, dd):
   
    """ get the days in a given year for calender date
    
    Input: 
    -------------------------------------------
    1. yyyy, mm, dd:<integer>:year, month and day 
    
    Returns:
    ---------------------------------------------
    1.  days in year
    
    """
        
    date=datetime.date(yyyy, mm, dd)
    date_ref=datetime.date(yyyy, 1, 1)
    delta=date-date_ref
    ndays=delta.days+1
    del date, date_ref, delta
    return ndays


def get_doy(yyyy, mm, dd):
   
    """ get the days in a given year for calender date
    (wrapper for day_of_year)

    Input: 
    ----------------------------
    1. yyyy, mm, dd:<integer>:year, month and day 
   
    Returns:
    ---------------------------------------------
    1.  days in year
    
    """
    return day_of_year(yyyy, mm, dd)


def tai85_to_utc(tai85):

    """ convert a time in tai85 format to UTC string 
    Inputs: 
    ----------------------------------------
    1. tai85:<integer>: Seconds since 1985-01-01 00:00:00
    
    Returns:
    ---------------------------------------------
    1. utc: <string>: time in utc format (yyyy-mm-dd hh:mi:ss)
    
    
    Notes:
    Time_leap is done in a hard (hand) way 
    """
    
    date_new=TAI0+datetime.timedelta(seconds=(tai85-second_leap))
    utc=str(date_new)
    del date_new
    return utc

def utc_to_tai85(utc):
    
    """ convert utc string to tai85 format 
    (i.e., Seconds since 1985-01-01 00:00:00)
    
    Inputs:
    --------------------------------------------
    1. utc: utc time in yyyy-mm-dd hh:mi:ss
    
    Returns:
    ---------------------------------------------
    1. tai85: <integer> seconds since 
    
    """

    t1,t2=utc.split()
    syyyy, smm, sdd=t1.split('-')
    shh, smi, ssec=t2.split(':')
    yyyy=int(syyyy)
    mm=int(smm)
    dd=int(sdd)
    hh=int(shh)
    mi=int(smi)
    fsec=float(ssec)
    sec=int(fsec)
    iv_time=datetime.datetime(yyyy, mm, dd, hh, mi, sec)-TAI0
    tai85=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    return tai85

def doy_to_tai85(doy, sec=0, yyyy=def_yyyy):
    """ convert the day of year to tai85 format 
    (i.e., Seconds since 1985-01-01 00:00:00)
    
    Inputs:
    1. doy:<integer>: Day of year
    2. sec:<integer>: Seconds in the day
    3. yyyy: <integer>: Year
    
    Returns:
    ---------------------------------------------
    1. tai85:  in seconds
    """

    date0=datetime.datetime(yyyy, 1, 1, 0,0,0)
    date0=date0+datetime.timedelta(days=int(doy)-1, seconds=int(sec))
    
    iv_time=date0-TAI0
    tai85=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    del date0
    return tai85

def doy_to_utc(doy, sec=0, yyyy=def_yyyy):
    
    """ convert day of year to utc string
        
    Inputs:
    1. doy:<integer>: Day of year
    2. sec:<integer>: Seconds in the day
    3. yyyy: <integer>: Year
    
    Returns:
    ---------------------------------------------
    1. utc: the time in utc format yyyy-mm-dd hh:mm:ss
    
    """
    
    if (doy>366):
        xdoy=doy
        doy=yyyy
        yyyy=xdoy
    
    date0=datetime.datetime(yyyy, 1, 1, 0,0,0)
    date0=date0+datetime.timedelta(days=int(doy)-1, seconds=int(sec))
    utc=str(date0)
    del date0
    return utc

def doy_to_time_array(doy, yyyy=2005):

    """ convert day of year to calendar date (yyyy, mm, dd)
    
    Inputs:
    1. doy:<integer>: Day of year
    2. yyyy: <integer>: Year
    
    Returns:
    ---------------------------------------------
    1. yyyy, mm, dd:<integer>: year, month and day
    
    """
    if (doy>366):
        xdoy=doy
        doy=yyyy
        yyyy=xdoy
    
    
    utc=doy_to_utc(doy, 0, yyyy)
    yyyy, mm,dd, hh, mi, sec=utc_to_time_array(utc)
    return yyyy, mm, dd

def utc_to_time_array(utc):
    """ convert the utc string to yyyy, mm, dd, hh, mi, sec
    
    Inputs:
    
    1. utc:<string>: the time in utc format yyyy-mm-dd hh:mm:ss
    

    Returns:
    -----------------------------------------
    1. yyyy, mm, dd, hh, mi, sec:<integer>: year, month, day, hour, minute, seconds
    """
    sd, sh=utc.split(' ')
    syyyy, smm, sdd=sd.split('-')
    shh, smi, ssec=sh.split(':')
    yyyy=int(syyyy)
    mm=int(smm)
    dd=int(sdd)
    hh=int(shh)
    mi=int(smi)
    sec=float(ssec)
    return yyyy, mm, dd, hh, mi, sec

def time_array_to_utc(yyyy,mm, dd, hh=0, mi=0, sec=0):
    """ convert yyyy, mm, dd, hh, mi, sec to utc string
 
    Inputs:
    1. yyyy, mm, dd, hh, mi, sec:<integer>: year, month, day, hour, minute and second
        
    Returns:
    -----------------------------------------
    1. utc: the time in utc format yyyy-mm-dd hh:mm:ss
    
    """
    sec=int(sec)
    
    syyyy_mm_dd=r'%4.4d-%2.2d-%2.2d' %(yyyy, mm, dd)
    if (sec>=10):
        shh_mi_sec=r'%2.2d:%2.2d:%5.2f' %(hh, mi, sec)
    else:
        shh_mi_sec=r'%2.2d:%2.2d:%4.2f' %(hh, mi, sec)
    return syyyy_mm_dd+' '+shh_mi_sec

    
def get_tau(yyyy, mm, dd, hh=0, mi=0, sec=0):
    """ convert yyyy, mm, dd, hh, mi, sec to tau (i.e, hours since 1985.01.01, 00:00)
    
    Inputs:
    1. yyyy, mm, dd, hh, mi, sec:<integer>: year, month, day, hour, minute and second
        
    Outputs:
    1. tau:<float>: Hours since 1985.01.01, 00:00
    """
    
    iv_time=datetime.datetime(yyyy, mm, dd, hh, mi, sec)-TAI0
    tau=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    tau=tau/3600.00
    return tau


def tau_to_time_array(tau):
    """ convert tau (i.e, hours since 1985.01.01, 00:00) to  yyyy, mm, dd, hh, mi, sec 
    
    Inputs:
    -----------------------------------------------------------
    1. tau:<float>: Hours since 1985.01.01, 00:00
    
    Returns:
    -----------------------------------------
    1. yyyy, mm, dd, hh, mi, sec:<integer>: year, month, day, hour, minute and second
    
    """
    # change seconds (TAI85)

    tai85=3600.0*tau
    utc=tai85_to_utc(tai85)
    yyyy, mm, dd, hh, mi, sec=utc_to_time_array(utc)
    return yyyy, mm, dd, hh, mi, sec




def next_doy(yyyy_in, doy_in, days=1, return_ymd=False):
    """ get date  for current time + n days later
    
    Inputs:
    -----------------------------------------
    
    1.yyyy_in, doy_in:<integer>: current year and day of year
    2.days:<integer>: day increment
    3.return_ymd:<T/F>: if ture, year, month and day will be returned 
    
    Returns:
    -----------------------------------------
    1.yyyy, doy:<integer>: year and day, If return_ymd==False, 
    
    or 
    1. yyyy, mm, dd:<integer>:year, month and day, If return_ymd==True.
    
    
    
    
    """
    utc=doy_to_utc(doy_in, 0, yyyy_in)
    yyyy, mm, dd, hh, mi, sec=utc_to_time_array(utc)
    date0=datetime.datetime(yyyy, mm, dd, hh,mi,sec)
    date0=date0+datetime.timedelta(days=int(days), seconds=int(sec))
    utc=str(date0)
    yyyy, mm, dd, hh, mi, sec=utc_to_time_array(utc)
    doy=day_of_year(yyyy, mm, dd)
    if (return_ymd):
        return yyyy, mm, dd
    else:
        return yyyy, doy
    
def get_ut_loc_time_slot(lt_st, lt_end, day_time_grid, lon, day_length=24.0*3600):
    
    """ get indexes for UT time series which fall within the local time slot
    
    Inputs:
    1.lt_st, lt_end:<numeric>: local time periods. 
    2.day_time_grid:<array>: UT time series
    3.lon:<array>: longitudes of the locations. 
    4. day_length:<float>: the length of day in seconds or hours. 
    The unit should be consistent with day_time_grid and the local time period.
    
    
    Returns:
    -----------------------------------------
    1. usd_idx:<array>:the index for the UT time series which fall within the local time
    between ls_st and lt_end 
    """
    
    tshift=day_length*lon/360.
    ut_st=lt_st-tshift
    ut_end=lt_end-tshift
    
    if (ut_end<0.0): # whole in earlier day 
        ut_st=ut_st+day_length
        ut_end=ut_end+day_length

        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)
        
    elif (ut_st>day_length): # whole in later day 
        ut_st=ut_st-day_length
        ut_end=ut_end+day_length
    
        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)
    elif (ut_st<0.0): # cross earlier day

        ut_st=ut_st+day_length
        chose1=(day_time_grid>=ut_st) & (day_time_grid<=day_length)
        chose2=(day_time_grid>=0.0) & (day_time_grid<ut_end)
        usd_idx=where(chose1 | chose2)
        usd_idx=squeeze(usd_idx)
    elif (ut_end>day_length): # cross later day
        
        ut_end=ut_end-day_length
        chose1=(day_time_grid>=ut_st) & (day_time_grid<=day_length)
        chose2=(day_time_grid>=0.0) & (day_time_grid<ut_end)
        usd_idx=where(chose1 | chose2)
        usd_idx=squeeze(usd_idx)
    else:
        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)

        
    return usd_idx
