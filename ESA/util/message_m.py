"""
   Functions for warning message processing 
   
   Authors: L. Feng, Edinburgh University
   History: v0.9, 2012.04.30
   History: v0.95, 2012.12.01
    
   
   1. Parameters
   -----------------------------------------------
   # action
   
   1. msm_action_stop=-3
   2. msm_action_wait=-2
   3. msm_action_warn=-1
   4. msm_action_cont=0
   
   # pre-defined error title 
   
   5. msm_no_attr='attribute not found: 
   6. msm_no_var='variable not found: '
   7. msm_wrong_dim='wrong array dimension: '
   8. msm_wrong_type='wrong type: '
   9. msm_wrong_index='Index is not valid:'
   

   2. Functions
   -------------------------------------
   1. show_err_msg  # show error message 
   2. write_log_msg # save error message to log files. 
   
"""


import sys

#<<< Parameters >>> 

msm_action_stop=-3
msm_action_wait=-2
msm_action_warn=-1
msm_action_cont=0

# pre-defined error title 

msm_no_attr='attribute not found: '
msm_no_var='variable not found: '
msm_wrong_dim='wrong array dimension: '
msm_wrong_type='wrong type: '
msm_wrong_index='Index is not valid:'



#<<< Parameters >>> 

def show_err_msg(error_msg, msg_title="", action=msm_action_warn):
    """ show error message 
    
    Inputs:
    ----------------------------------
    1. error_msg:<str>: error message 
    2. msg_title:<str>: title for error message
    3. action:<integer>: action ID
    
    Returns:
    ----------------------------------------
    1. 
    """
    # get module name and the function name 
    
    modulename=sys._getframe(1).f_code.co_filename
    callername=sys._getframe(1).f_code.co_name
    
    print ''

    if (action==msm_action_warn):
        print 'W***     Warning     ***W'
    elif (action==msm_action_stop):
        print 'E***     Error       ***E'
    elif (action==msm_action_wait):
        print 'A***     Action Need    ***A'
    else:
        print 'T***     Troubles    ***T'
        
    print '1. Module:',  modulename
    print '2. Function:', callername

    if (msg_title==""):
        print '3. Error Message:',  error_msg
    else:
        print '3. Error Message:',  msg_title, error_msg
    
    print '***----------------***'
    
    if (action==msm_action_warn): # just warning 
        pass
    elif (action==msm_action_stop): # break
        exit()
        
    elif (action==msm_action_wait): # user to decide whether the process will be continued 
        while (True):
            choice_yn=raw_input('Continue ? (yes/no)')
            choice_yn=choice_yn.lower()
            choice_yn=choice_yn.strip()
            
            if (choice_yn=='yes'):
                break
            elif (choice_yn=='no'):
                exit()
            else:
                pass
    else:
        pass
    
    print ''
    
    

def write_log_msg(logflnm, log_msg):
    """ 
    Save error message into log file 
    
    Inputs:
    ---------------------------------------------------
    1. logflnm:<str>: log file name 
    2. log_msg:<str>: message to be saved
    
    
    """
    
    logfile=open(logflnm, 'a')
    modulename=sys._getframe(1).f_code.co_filename
    callername=sys._getframe(1).f_code.co_name
    
    line_1='# <'+modulename+'> '+'<'+callername+'> '+':'+log_msg+'\n'
    logfile.write(line_1)
    logfile.close()
    

    
