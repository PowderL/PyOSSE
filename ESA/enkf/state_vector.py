"""
    Class for state vector
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. ens_stat_cl: Class for state vector 
    
"""


import ESA.util.otool_gcfile_io as gcfio 
import ESA.util.time_module as tm
import ESA.util.gen_plots as gpl
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import numpy as npy
import pylab as plt 
import etkf_half as et
import etkf_cor as etc

class ens_stat_cl:
    
    """
    Class for state vector used to assimilate observations

    Members:
    -------------------------------------------------------------
    1. attr_dict:<dict>: attributes
    2. step_tag_lst:<list>: ID for run step 
        
    
    # I. variables at current time window
    
    3. wnd_mean_x0:<array, (wnd_nx)>: prior coefficient values for lag window 
    4. wnd_mean_x:<array, (wnd_nx)>: posterior (prior) coefficient values for lag window
    5. wnd_dx:<array, (wnd_nx, wnd_ne)>: perturbations ensemble for coefficients
    6. wnd_inc_m:<array, (wnd_ne>:  increment matrix for ETKF data assimilation 
    7. wnd_xtm: <array, (wnd_nx, wnd_ne)>: transform matrix for ETKF data assimilation 
    8. wnd_xinc: <array, (wnd_nx)>:  increment for ETKF data assimilation
     

    # II: variables for the step 
    
    9.  inc_m:<array, (wnd_nx)>: increment matrix for current step 
    10. xtm: <array, (wnd_nx, wnd_ne)>: transform matrix for current step
    11. xinc: <array, (wnd_nx)>:  increment for ETKF data assimilation
     
    
    # III: list of objects:   
    
    12. x2flux_lst: <list, t:obj>: class for projecting x back to flux
    13. wnd_mean_x_lst:<list, t:array>: list of wnd_mean_x at beginning of each step
    14. wnd_mean_x0_lst:<list, t:array>: list of wnd_mean_x0 at beginning of each step
    15. wnd_xtm_lst:<list, t:array>: list of wnd_xtm at beginning of each step
    16. wnd_dx_lst:<list, t:array>: list of wnd_dx at beginning of each step
    

    # IV: Settings 
        
    17. cur_yyyy:<int>: starting year
    18. cur_mm:<int>: starting month
    19. cur_dd:<int>: starting dd
    20. cur_doy:<int>: starting doy 
    21 .next_yyyy:<int>: next year
    
    22. max_doy:<int>: max doy 
    23. cur_step:<int>: time step
    24. wnd_step:<int>: step inside window
    25. max_step:<int>: maximum step
    26. lag_window:<int>: length of lag_window (max steps inside the window)
    
    # V: time and size
        
        
    27. tau_st_lst:<list, t:float>: starting time
    28. tau_end_lst:<list, t:float>: end time 
    29. x_st_lst:<list, t:int>: starting position of the current period when it becomes the first of the window
    30. x_end_lst: <list, t:int>:  ending position of the current period when it becomes the first of the window
    31. en_st_lst:<list, t:int>: ensemble starting position of the current period when it becomes the first of the window
    32. en_end_lst: <list, t:int>:  ensemble ending position of the current period when it becomes the first of the window
        
    33. nbias:<int>: size of bias correction terms
    34. mod_bias:<array, (nbias)>: first guess of bias correction terms 
    35. bias_err:<array, (nbias)>: uncertainty of bias correction terms 
    . 

    # VI scaling matrix (or spatial correction)
    
    36. wnd_scale_m:<array, (wnd_ne, wnd_ne)>: 
    
    # VII:  # temporal correlation 
        
    37. use_tcor:<T/F>: whether to include temporal correlation 
    38. tcor_len:<int>: length of temporal correlation 
    39. tcor_factor:<float>: scaling factor for temporal correlation 
    40. tcor_inc_m:<wnd_ne, wnd_ne>: increment matrix due to last temporal correlation 
    41. wnd_tcor_inc_m::<wnd_ne, wnd_ne>: increment matrix due to all temporal correlations inside the window
        
    42. tcor_m:<wnd_ne, wnd_ne>: correlation matrix due to last temporal correlation
    43. wnd_tcor_m:<wnd_ne, wnd_ne>: correlation matrix due to  temporal correlations inside the window
    44. wnd_nx, wnd_ne:<int>: size of self.wnd_dx
    
    
    Functions
    ==========================================================
    1. __init__: initialization
    2. join_vector: joint two arrays together 
    3. fifo_vector: push out  first part of long array 
    4. join_matrix: join two matrices along diagonal together 
    5. fifo_matrix: push out first part of a matrix
    6. construct_tcor_matrix:  Create temporal correlation matrix 
    7. add_new_x_to_window:  Add new apriori to assimilation window 
    8. init_x_dx: initializing for first time period 
    9. append_x_dx:  Append prior variable to the end of lag window
    10. do_assim:  assimilate observations and update state vector
    11. vf_construct_aprior: construct x and dx
    12. shift_y_dy: update model mean and model uncertainty in observation space
    13. fifo_x_dx: Append variable to the end, and remove the first one from window
    14. convolve_y_dy:  tranfer  model mean_y and transform dy to reflect changes on mean_x and dx.
    
    
        
    Notes:
    ==============================================================
    1. A sequential data assimilation approach is assumed in this module. y and dy must be kept (by convolve_y_dy_
    must be consistent with x and dx (and wnd_xtm)
    
        

    """
    
    def __init__(self, \
                     yyyy_st,\
                     mm_st, \
                     dd_st,\
                     max_step,\
                     lag_window,\
                     mod_bias=[],\
                     mod_bias_err=[],\
                     use_tcor=False,\
                     tcor_len=30.0*24,\
                     tcor_factor=0.5,\
                     xnorm=1.0,\
                     **keywords):
        """
        
        Inputs:
        =========================================
        1. yyyy_st:<int>: starting year
        2. mm_st:<int>: starting month
        3. dd_st:<int>: starting day
        4. max_step:<int>: maximum step
        5. lag_window:<int>: steps inside the lag window 
        6. mod_bias:<array/list>:  bias to be estimated  
        7. mod_bias_err:<array/list>:  mod_bias_uncertainty 
        8. use_tcor:<T/F>: whether to turn on  temporal correlation 
        9. tcor_len:<float>: temporal correlation length
        10. tcor_factor:<float>: scaling factor for temporal correlation 
        11. keyword:<dict>: attributes
        
        """
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_stv})
        self.step_tag_lst=[]
        
        # S1:  array for data storage at current time window
        
        
        self.wnd_mean_x0=None # prior coef values
        self.wnd_mean_x=None  # posterior coef value 
        self.wnd_dx=None   # prior dx after scaling 
        
        self.wnd_inc_m=None
        self.wnd_xtm=None  # current transform matrix for ETKF approach
        self.wnd_xinc=None
        self.xnorm=xnorm
        
        
        # S2: variables for the step 
        
        self.inc_m=None # current increment matrix  
        self.xtm=None  # analysis X increment
        self.xinc=None  # analysis X increment
        
        # S3: list of class object for projecting back to flux 
        
        self.x2flux_lst=[]
        
        # S2: array for data storage for historical time window
        
        self.wnd_mean_x_lst=[]
        self.wnd_mean_x0_lst=[]
        self.wnd_xtm_lst=[]
        self.wnd_dx_lst=[]
        
        
        
        # setting 
        
        self.cur_yyyy=yyyy_st
        self.cur_mm=mm_st
        self.cur_dd=dd_st
        
        self.cur_doy=tm.day_of_year(self.cur_yyyy, self.cur_mm, self.cur_dd)
        self.max_doy=-1
        
        self.cur_step=0
        self.wnd_step=0
        
        self.max_step=max_step
        self.lag_window=lag_window
        
        #  time and size
        
        
        self.tau_st_lst=[]
        self.tau_end_lst=[]
        
        
        self.x_st_lst=[]
        self.x_end_lst=[]
        
        self.en_st_lst=[]
        self.en_end_lst=[]
        
        self.next_yyyy=self.cur_yyyy
        
        
        self.nbias=len(mod_bias)
        self.mod_bias=npy.array(mod_bias)
        self.bias_err=npy.array(mod_bias_err)
        
        # scaling matrix 
        self.wnd_scale_m=None

        # values for total of window 

        # it will be initialized with zeros when set times 
        
        
        # temporal correlation 
        
        self.use_tcor=use_tcor
        self.tcor_len=tcor_len 
        self.tcor_factor=tcor_factor
        
        self.tcor_inc_m=None  #  
        self.wnd_tcor_inc_m=None  #  
        
        # correlation 
        self.tcor_m=None
        self.wnd_tcor_m=None
        
        self.wnd_nx=0
        self.wnd_ne=0
        
        
    def join_vector(self, x1, x2):
        """
        joint two arrays together 
        
        Inputs:
        ---------------------------------
        1. x1:<array, (old_nx)>: array one
        2. x2:<array, (nx)>: array two 

    4    Returns
        --------------------------------
        1. new_x:<array, (old_nx+nx)>: jointed array
        """
        
        new_x=npy.concatenate((x1, x2))
        
        return new_x
   

    def fifo_vector(self, x1, x2, retain_idx):
        """
        replace first parts of x1  
        
        Inputs:
        ---------------------------------
        1. x1:<array, (old_nx)>: array one
        2. x2:<array, (nx)>: array two 
        3. retain_idx:<array, (nsel)>: index for elements to be retained 
        
        Returns
        --------------------------------
        1. new_x:<array, (nse+nx)>: jointed array
        """
        
        new_x=x1[retain_idx]
        new_x=npy.concatenate((new_x, x2))
        
        return new_x
    


    def join_matrix(self, x1, x2):
        """
        join two matrix along diagonal together 
        
        Inputs:
        ---------------------------------
        1. x1:<array, (old_nx, old_ny)>: matrix one 
       2. x2:<array, (nx,ny)>: matrix  two 

        Returns
        --------------------------------
        1. new_x:<array, (old_nx+nx, old_ny+ny)>: jointed matrix along diagonal 
        
        """
        
        old_nx, old_ny=npy.shape(x1)
        nx, ny=npy.shape(x2)
        
        # #c: the shape (in old_nx, old_ne+ne)
        new_x=npy.zeros([old_nx+nx, old_ny+ny], float)

        # #c: fill the old parts
        new_x[0:old_nx, 0:old_ny]=x1[0:old_nx, 0:old_ny]

        # #c: fill the new parts 
        new_x[old_nx:old_nx+nx, old_ny:old_ny+ny]=x2[0:nx, 0:ny]
        
        
        return new_x
    
    

    def fifo_matrix(self, x1, x2, retain_xid, retain_yid):
        """
        relpace the first part of  matrix  
        
        Inputs:
        ---------------------------------
        1. x1:<array, (old_nx, old_ny)>: matrix one 
        2. x2:<array, (nx,ny)>: matrix  two 
        3. retain_xid:<array, (nsel_x)>: x index for values to be kept 
        4. retain_yid:<array, (nsel_y)>: y index for values to be kept
        

        Returns
        --------------------------------
        1. new_x:<array, (nse_x+nx, nsel_y+ny)>: jointed matrix along diagonal 
        
        """
        # #c: keep the selected parts 
        
        new_x=x1[retain_xid,:]
        new_x=new_x[:, retain_yid]
    
        # #c: join retained new_x with x2

        
        new_x=self.join_matrix(new_x, x2)
        
        return new_x
    
    
    

    
        
    def construct_tcor_matrix(self, dx, tau_st, first_shift=1, use_wnd_inc_m=True):
        
        """
        Create temporal correlation matrix 

        Inputs:
        ---------------------------------
        1. dx:<array, (nx,ne)>: sampling x
        2. tau_st:<float>: starting time 
        3. keep_first:<int>: if first_shift is the starting step inside (previous) window
        4. use_wnd_inc_m:<T/F>:  if True, wnd_inc_m will be considered to be able to reproduce 
        mean_x-mean_x0=self.wnd_dx*wnd_inc_m
        
        returns:
        
        --------------------------------------
        1. tcor_m:<array, (ne, wnd_ne-first_ne)>: temporal correlation matrix 
        2. tcor_scale:<array, (ne, ne)>:  scaling matrix 
        3. tcor_inc_m:<array, (ne)>:  tmporal correlation increment 
        
        
        Notes: 
        
        1. this functions are only approximations 
        Temporal correlation between old (A) and new (B) is  in the form of               
        
        | A  0           |              | A   0 |    | I        0 |
        | B*TS*TC   B*TS |------------> | 0   B |  X | TS*TC   TS |
        TS=tcor_scale: re-scaling  (ne, ne) 
        TC=tcor_m: correlation matrix (ne, old_ne) 
        tcor_inc_m: increment matrix due to correlation:
        tcor_inc_m=TC*inc_m(A). inc_m(A) is the increment matrix for old A 
        tcor_x_inc=B*TS*tcor_inc_m
        
        """
        
        # S1: check size and initialize 
        
        
        nx, ne=npy.shape(dx)
        wnd_nx, wnd_ne=npy.shape(self.wnd_dx)
                
        # c: step for first one inside the window
        

        st_step=self.cur_step-self.wnd_step
        
        first_ne=0
        for ishift in range(0, first_shift):
            enst=self.en_st_lst[st_step]
            ened=self.en_end_lst[st_step]
            first_ne=first_ne+ened-enst
            
        tcor_m=npy.zeros([ne, self.wnd_ne-first_ne], float)
        
        tcor_inc_m=npy.zeros(ne, float)
        tcor_scale=npy.identity(ne)
 
        if (not self.use_tcor):  
            return tcor_m, tcor_scale, tcor_inc_m
        
        
        # #c: index for X element to be kept
        
        tcor_sum=1.0 # it self
        
        print 'st_step', st_step
        print 'wnd_step', self.wnd_step
        
        for istep in range(first_shift, self.wnd_step):
            
            dt=tau_st-self.tau_st_lst[istep+st_step]
            tfactor=dt/self.tcor_len
            tfactor=self.tcor_factor*npy.exp(-tfactor)
            
            # #c: ensemble location 
            
            
            enst=self.en_st_lst[istep+st_step]
            ened=self.en_end_lst[istep+st_step]
            
            new_tcor_m=npy.zeros([ne, ened-enst], float)
            
            # #c: new_tcor in shape of [ne, old_ne]
            # |tfactor,  0      |
            # |  0    , tfactor |
            
            for iex in range(ne):
                if (iex<(ened-enst)):
                    new_tcor_m[iex, iex]=tfactor
            
                
            
            if (istep==first_shift):
                tcor_m=npy.array(new_tcor_m)
            else:
                tcor_m=npy.concatenate((tcor_m, new_tcor_m), axis=1)
                
            tcor_sum=tcor_sum+tfactor*tfactor
        # loop k end
            
        # scaling the value so that 
        tcor_sum=1.0/npy.sqrt(1.0+tcor_sum)
            
        scaled_tcor_m=tcor_sum*tcor_m
        tcor_scale=tcor_sum*tcor_scale
        
        if (self.nbias>0):
            # add non-correlation with bias to tcor_m
            bias_cor=npy.zeros([ne, self.nbias], float)
            scaled_tcor_m=npy.concatenate((bias_cor, scaled_tcor_m), axis=1)
            tcor_m=npy.concatenate((bias_cor, tcor_m), axis=1)
        
        
        if (use_wnd_inc_m):
            
            # #c: as a  approximation xmean-x_mean0=self.wnd_dx*self.wnd_inc_m
            
            wnd_inc_m=self.wnd_inc_m[first_ne:]
            # #c: note that not the scaled_tcor_m
            print npy.shape(tcor_m), npy.shape(wnd_inc_m)
            
            tcor_inc_m=npy.dot(tcor_m, wnd_inc_m)
            
            
        else:
        
            #  #c: It is more accurate when nx=ne 
            
            # solve inc_m for 
            # self.wnd_dx*base_inc_m=self.wnd_mean_x-self.wnd_mean_x0
            
            # this inc_m may be different from self.wnd_inc_m
            
            
            base_dx=self.wnd_dx
            base_dx=base_dx[self.nbias+first_ne:, self.nbias+first_ne:]
            
            
            wnd_xinc=self.wnd_mean_x-self.wnd_mean_x0
            wnd_xinc=wnd_xinc[self.nbias+first_ne:]
            wnd_inc_m=nlg.solve(base_dx, wnd_xinc)
            
            # #c: change to mean_x due to temporal correlation inc_m  
            
            tcor_inc_m=dot(tcor_m, wnd_inc_m)
            
        return scaled_tcor_m, tcor_scale, tcor_inc_m
    
        
    def add_new_x_to_window(self, apr_x, \
                                apr_dx, \
                                tau_st,\
                                tau_end,\
                                step_tag="new",\
                                enlarge_factor=1.0,\
                                scale_m=None, \
                                lag_window=None,\
                                cl_x2flux=None\
                                ):
        
        
        """
        
        Add new apriori to assimilation window 
        
        Inputs:
        
        ============================================
        1. apr_x:<array>: prior x
        2. apr_dx:<array, (nx, ne)>: prior dx
        
        3. tau_st, tau_end:<float>: start and end time 
        
        3. enlarge_factor:<float>: factor used to enlarge xtm 
        4. scale_m:<array, (ne, ne)>: matrix for scaling dx, 
        
        which can be spatial correlation 
        between 'basis functions' 
        
        5. lag_window:<integer>: lag window length 
        
        6. cl_x2flux:<x2flux_cl>: project coefficient back to flux
        
        
        
        Returns:
        ============================================
        
        1. wnd_nx:<int>: size of x inside window
        2. wnd_ne:<int>: size of ensemble inside window 
        
        Notes:
        --------------------------------------------
        1. temporal correlation  is assumed to be in the shape  of:

              |I   0     I |
              |TS*TC, TS*I |
         
         TC is the correlation matrix, TS is scaling matrix. 

         
         
         
        
        """
        
        # S1: size

        
        dx=apr_dx
        mean_x=apr_x
        
        nx,ne=npy.shape(dx)

        # S2: scaling dx by scale_m
        
        if (scale_m==None):
            scale_m=npy.identity(ne)
            
        
        dx=npy.dot(dx, scale_m)
        

        
                       
        if (lag_window<>None):
            self.lag_window=lag_window
        

        if (self.cur_step==0):
        
            
            # T1: keep record

            self.wnd_mean_x0_lst=[]
            self.wnd_mean_x_lst=[]
            self.wnd_dx_lst=[]
            self.wnd_xtm_lst=[]
            self.x2flux_lst=[]
            
            
            # T2: insert new variables
            self.init_x_dx(mean_x, dx, tau_st, tau_end, \
                               enlarge_factor, scale_m)
           

            self.x2flux_lst.append(cl_x2flux)

        elif (self.wnd_step<self.lag_window):
        
            self.wnd_mean_x0_lst.append(self.wnd_mean_x0)
            self.wnd_mean_x_lst.append(self.wnd_mean_x)
            self.wnd_dx_lst.append(self.wnd_dx)
            self.wnd_xtm_lst.append(self.wnd_xtm)
            self.x2flux_lst.append(cl_x2flux)
            
            # #c: append 
            self.append_x_dx(mean_x, dx, tau_st, tau_end, \
                                 enlarge_factor, scale_m)
            

        else:
            self.wnd_mean_x0_lst.append(self.wnd_mean_x0)
            self.wnd_mean_x_lst.append(self.wnd_mean_x)
            self.wnd_dx_lst.append(self.wnd_dx)
            self.wnd_xtm_lst.append(self.wnd_xtm)
            self.x2flux_lst.append(cl_x2flux)
            
            # #c: first in first out 
            self.fifo_x_dx(mean_x, dx, tau_st, \
                               tau_end, enlarge_factor, scale_m)
            
        # set temoral correlation etc
        # self.
        self.tau_st_lst.append(tau_st)
        self.tau_end_lst.append(tau_end)
        self.step_tag_lst.append(step_tag)
        
        # set tcor_inc_m
            
        self.wnd_nx, self.wnd_ne=npy.shape(self.wnd_dx)
        
        
        
        return self.wnd_nx, self.wnd_ne
    
    
    def init_x_dx(self, mean_x, dx, tau_st, tau_end, enlarge_factor, scale_m):
        
        """
        The first time period o be included 
        It will be called when wnd_step<lag_window
        
        1. Inputs:
        -----------------------------------------------------
        1. mean_x:<array, (nx)>: prior mean_x
        2. dx:<array, (nx, ne)>: 'prior' dx
        3. enlarge_factor:<float>: factor used to enlarge xtm 
        4. scale_m:<array, (ne, ne)>: matrix for scaling dx, which can be spatial correlation 
           between regions
           
        """
        
        
        # --------------------------------------------
        # structure of x 
        # 0--nbias-1:  bias  in observations
        # nbias--nbias+nx: the coefficients for regional perturbations in period  1
        # nbias+nx---nbias+nregion+nx: coefficients for regional perturbations 2
        # ------
        # 
        # error covariance
        # prior x 
        # P=dot(dx, transpose(dx))
            
        #   
        #  posteriori x
        #  pdx=dot(x, xtm)
        #  P=dot(pdx, transpose(pdx))
        #  
        #--------------------------------------------------
        
            
        # S1: variable for obs (model) bias 
        


        nx, ne=npy.shape(dx)
        
        # #c: bias 
        
        x_bias=self.mod_bias
        dx_bias=self.bias_err 
        dx_bias=npy.diag(dx_bias)
        
        # #c: window variable sizes 
        
        wnd_nx=self.nbias+nx
        wnd_ne=self.nbias+ne
        
        # #S2: set window variable 
        
        self.wnd_dx=self.join_matrix(dx_bias, dx)
        self.wnd_mean_x0=self.join_vector(self.mod_bias, mean_x)
        
        self.wnd_mean_x=npy.array(self.wnd_mean_x0)
        self.wnd_xtm=npy.identity(wnd_ne,float)
        self.wnd_inc_m=npy.zeros(wnd_ne, float)  # increment matrix for the time window
        self.wnd_xinc=npy.zeros(wnd_nx, float)  #  x increment inside the time window 
        
        
        
        self.cur_step=self.cur_step+1
        self.wnd_step=self.wnd_step+1
        
        
            
        # S4: set scaling matrix:
        
        self.wnd_scale_m=npy.identity(wnd_ne,float)
        
        self.wnd_scale_m[self.nbias:wnd_ne, self.nbias:wnd_ne]=scale_m[0:ne, 0:ne]
        
        
        # S5: time correlation matrices
        
        # #T1: window 
        self.wnd_tcor_m=npy.identity(wnd_ne,float)
        self.wnd_tcor_inc_m=npy.zeros(wnd_ne, float)
        
        # #T2: step one 
        
        self.tcor_m=npy.identity(wnd_ne,float)
        self.tcor_inc_m=npy.zeros(wnd_ne, float)
        
        
        # S6: set step variables to zeros or 1 
        
        
        self.inc_m=npy.zeros(wnd_ne, float)
        self.xinc=npy.zeros(wnd_ne, float)
        self.xtm=npy.identity(wnd_ne)
        
        # S7: step control 
        
        self.x_st_lst.append(self.nbias)
        self.x_end_lst.append(self.nbias+nx)
        self.en_st_lst.append(self.nbias)
        self.en_end_lst.append(self.nbias+ne)
        
        return npy.shape(self.wnd_dx)

                  
    
    def append_x_dx(self, mean_x, dx, tau_st, tau_end, enlarge_factor, scale_m):
        
        """
        Append variable to the end
        It will be called when wnd_step<lag_window
        
        1. Inputs:
        -----------------------------------------------------
        1. mean_x:<array, (nx)>: mean_x
        2. dx:<array, (nx, ne)>: 'prior' dx
        3. tau_st, tau_end:<float>: start and end time 
        4. enlarge_factor:<float>: factor used to enlarge xtm 
        5. scale_m:<array, (ne, ne)>: matrix for scaling dx, which can be spatial correlation 
           between regions
           
        
        Returns:
        -----------------------------------------------------
        1. shape(self.wnd_dx):<array, (wnd_nx, wnd_ne)>: shape of the dx 
        
        
        """


        # S1: check size
        
        # #T1: size for new & old  variables
        
        nx, ne=npy.shape(dx)
        old_nx, old_ne=npy.shape(self.wnd_dx)
        wnd_nx=old_nx+nx
        wnd_ne=old_ne+ne
        
        # #T2:  set variables for current time step 
        # #c: no analysis inc & xtm=I before any bservation is assimilated
        
        inc_m=npy.zeros(ne, float)
        xinc=npy.zeros(nx, float)
        xtm=npy.identity(ne, float)
        
        # #T3: temporal correlation 
        # #:  temporal correlation between old (A) and new (B) is  in the form of               
        # #:  | A  0           |              | A   0 |    | I        0 |
        # #:  | B*TS*TC   B*TS |------------> | 0   B |  X | TS*TC   TS |
        # #:  TS=tcor_scale: re-scaling  (ne, ne) 
        # #:  TC=tcor_m: correlation matrix (ne, old_ne) 
        # #:  tcor_inc_m: increment matrix due to correlation:
        # #:  tcor_inc_m=TC*inc_m(A). inc_m(A) is the increment matrix for old A 
        # #:  tcor_x_inc=B*TS*tcor_inc_m
        
        tcor_m, tcor_scale, tcor_inc_m=self.construct_tcor_matrix(dx, tau_st, 0)
        
        # #T4: set temporal correlation matrices
        # #c: now tcor_m to be enarged to wnd_ne X wnd_ne, and include both TS*TC and TS 
        # ##: current step
        
        self.tcor_m=npy.identity(wnd_ne, float)  # enlarged to wnd_ne X wnd_ne 
        self.tcor_m[old_ne:wnd_ne, 0:old_ne]=tcor_m[0:ne, 0:old_ne]
        self.tcor_m[old_ne:wnd_ne, old_ne:wnd_ne]=tcor_scale[0:ne, 0:ne]
        
        # ##: window 
        self.wnd_tcor_m=self.join_matrix(self.wnd_tcor_m, xtm)
        self.wnd_tcor_m=npy.dot(self.wnd_tcor_m, self.tcor_m)
        
        
        # #c: tcor_inc_m  
        
        # ##: increment matrix due to temporal correlation 
        # ##: current one 
        
        self.tcor_inc_m=npy.zeros(wnd_ne, float)
        self.tcor_inc_m[old_ne:wnd_ne]=tcor_inc_m[0:ne]
        
        # #c: window one
        self.wnd_tcor_inc_m=self.join_vector(self.wnd_tcor_inc_m, tcor_inc_m)
        
        
        # S2: prior 
        
        # #T1: mean 
        self.wnd_mean_x0=self.join_vector(self.wnd_mean_x0, mean_x)
        
        # #T2: dx
        self.wnd_dx=self.join_matrix(self.wnd_dx, dx)
        
        # #c: apply temporal correlation 
        
        self.wnd_dx=npy.dot(self.wnd_dx, self.tcor_m)

        
        # S3: current posterior (which is also for prior at current step)
        
        # #T1: mean 
        print 'wnd_mean_x', npy.shape(self.wnd_mean_x)
        print 'mean_x', npy.shape(mean_x)
        
        self.wnd_mean_x=self.join_vector(self.wnd_mean_x, mean_x)
        
        # #T2: xtm 
        
        self.wnd_xtm=enlarge_factor*self.xtm
        self.wnd_xtm=self.join_matrix(self.wnd_xtm, xtm)
        
        # #T3: increment 
        
        self.wnd_inc_m=self.join_vector(self.wnd_inc_m, inc_m)
        self.wnd_xinc=self.join_vector(self.wnd_xinc, xinc)
        
        tcor_wnd_xinc=npy.dot(self.wnd_dx, self.tcor_inc_m)
        # #c: adding increment due to temporal correlation with previous stat
        
        self.wnd_mean_x=self.wnd_mean_x+tcor_wnd_xinc
        self.wnd_xinc=self.wnd_xinc+tcor_wnd_xinc
        
        
        # S4: 
        
        self.wnd_scale_m=self.join_matrix(self.wnd_scale_m, scale_m)

        
        self.cur_step=self.cur_step+1
        self.wnd_step=self.wnd_step+1
        
        # S5: position 

        # it is the size and position when the current period move to head of window 
        
        self.x_st_lst.append(self.nbias)
        self.x_end_lst.append(self.nbias+nx)
        self.en_st_lst.append(self.nbias)
        self.en_end_lst.append(self.nbias+ne)
        
        # S6: reset step inc_m,  xinc, and xtm
        
        
        self.inc_m=npy.zeros(wnd_ne, float)  # increment matrix for the current step 
        self.xinc=npy.zeros(wnd_nx, float)  #  x increment for the current step 
        self.xtm=npy.identity(wnd_ne)
        
        
        # S9: size 
        return npy.shape(self.wnd_dx)
    
    
    def fifo_x_dx(self, mean_x, dx, tau_st, tau_end, enlarge_factor, scale_m):
        
        """
        Append variable to the end, and remove the first one from window
        It will be called when wnd_step<lag_window
        
        1. Inputs:
        -----------------------------------------------------
        1. mean_x:<array, (nx)>: mean_x
        2. dx:<array, (nx, ne)>: 'prior' dx
        3. enlarge_factor:<float>: factor used to enlarge xtm 
        4. scale_m:<array, (ne, ne)>: matrix for scaling dx, which can be spatial correlation 
           between regions
        5. cl_x2flux:<x2flux_cl>: project coefficient back to flux
        
        """


        # S1: check size
        
        # #c: size for new dx
        nx, ne=npy.shape(dx)
        
            
        old_nx, old_ne=npy.shape(self.wnd_dx)
        

        # #T2:  set variables for current time step 
        # #c: no analysis inc & xtm=I before any bservation is assimilated
        
        inc_m=npy.zeros(ne, float)
        xinc=npy.zeros(nx, float)
        xtm=npy.identity(ne, float)

        
        # #c: position of the oldest coefficient set 
        
        istep=self.cur_step-self.wnd_step
        
        # #c: index for X element to be kept
        
        nxst=self.x_st_lst[istep]
        
        # #c: index for ne to be kept 
        nxed=self.x_end_lst[istep]
        
        retain_x_idx=range(0, nxst)+range(nxed, old_nx)
        retain_x_idx=npy.array(retain_x_idx)
            
        # #c: index for ensemble element to be kept
        
        enst=self.en_st_lst[istep]
        ened=self.en_end_lst[istep]
            
        retain_en_idx=range(0, enst)+range(ened, old_ne)
        retain_en_idx=npy.array(retain_en_idx)
        
        # #c: reset old_nx, old_ne
        
        old_nx=npy.size(retain_x_idx)
        old_ne=npy.size(retain_en_idx)
        
        # #c: window size
        
        wnd_ne=old_ne+ne
        wnd_nx=old_nx+nx
        
        
        # T3: set temporal correlation 
        
        
        # #T3: temporal correlation 
        # #:  temporal correlation between old (A) and new (B) is  in the form of               
        # #:  | A  0           |              | A   0 |    | I        0 |
        # #:  | B*TS*TC   B*TS |------------> | 0   B |  X | TS*TC   TS |
        # #:  TS=tcor_scale: re-scaling  (ne, ne) 
        # #:  TC=tcor_m: correlation matrix (ne, old_ne) 
        # #:  tcor_inc_m: increment matrix due to correlation:
        # #:  tcor_inc_m=TC*inc_m(A). inc_m(A) is the increment matrix for old A 
        # #:  tcor_x_inc=B*TS*tcor_inc_m
        
        tcor_m, tcor_scale, tcor_inc_m=self.construct_tcor_matrix(dx, tau_st, 1)
        
        
        # #T4: set temporal correlation matrices
        # #c: now tcor_m to be enarged to wnd_ne X wnd_ne, and include both TS*TC and TS 
        # ##: current one 
        
        self.tcor_m=npy.identity(wnd_ne, float)  # enlarged to wnd_ne X wnd_ne 
        self.tcor_m[old_ne:wnd_ne, 0:old_ne]=tcor_m[0:ne, 0:old_ne]
        self.tcor_m[old_ne:wnd_ne, old_ne:wnd_ne]=tcor_scale[0:ne, 0:ne]
        
        # ##: window one
        
        self.wnd_tcor_m=self.fifo_matrix(self.wnd_tcor_m, xtm, retain_en_idx, retain_en_idx)
        self.wnd_tcor_m=npy.dot(self.wnd_tcor_m, self.tcor_m)
                
        
        # #c: tcor_inc_m  
        
        # ##: increment matrix due to temporal correlation 
        # ##: current one 
        
        self.tcor_inc_m=npy.zeros(wnd_ne, float)
        self.tcor_inc_m[wnd_ne-ne:wnd_ne]=tcor_inc_m[0:ne]
        
        # #c: window one
        
        self.wnd_tcor_inc_m=self.fifo_vector(self.wnd_tcor_inc_m, tcor_inc_m, retain_en_idx)

         # S2: prior 
        
        # #T1: mean 
        self.wnd_mean_x0=self.fifo_vector(self.wnd_mean_x0, mean_x, retain_x_idx)
        
        # #T2: dx
        self.wnd_dx=self.fifo_matrix(self.wnd_dx, dx, retain_x_idx, retain_en_idx)
        
        # #c: apply temporal correlation 
        
        self.wnd_dx=npy.dot(self.wnd_dx, self.tcor_m)
        
        
        # S3: current posterior (which is also for prior at current step)
        
        # #T1: mean 
        
        self.wnd_mean_x=self.fifo_vector(self.wnd_mean_x, mean_x, retain_x_idx)
        
        # #T2: xtm 
        
        self.wnd_xtm=enlarge_factor*self.xtm
        self.wnd_xtm=self.fifo_matrix(self.wnd_xtm, xtm, retain_en_idx, retain_en_idx)
        
        # #T3: increment 
        
        self.wnd_inc_m=self.fifo_vector(self.wnd_inc_m, inc_m, retain_en_idx)
        self.wnd_xinc=self.fifo_vector(self.wnd_xinc, xinc, retain_x_idx)
        
        tcor_wnd_xinc=npy.dot(self.wnd_dx, self.tcor_inc_m)

        
        self.wnd_mean_x=self.wnd_mean_x+tcor_wnd_xinc
        self.wnd_xinc=self.wnd_xinc+tcor_wnd_xinc
        
        
        # S4: 
        
        self.wnd_scale_m=self.fifo_matrix(self.wnd_scale_m, scale_m, retain_en_idx, retain_en_idx)

        
        self.cur_step=self.cur_step+1
        self.wnd_step=self.wnd_step #  no move 
        
        
        # S5: position 

        # #c: it is the position for first period inside the window
        # it is the size and position for when the current period move to head of window 
        
        self.x_st_lst.append(self.nbias)
        self.x_end_lst.append(self.nbias+nx)
        self.en_st_lst.append(self.nbias)
        self.en_end_lst.append(self.nbias+ne)
        
        # S6: reset step inc_m and wnd_xinc
        
        self.inc_m=npy.zeros(wnd_ne, float)  # increment matrix for the current step 
        self.xinc=npy.zeros(wnd_nx, float)  #  x increment for the current step 
        self.xtm=npy.identity(wnd_ne)
        
        
        # S9: size 
        return npy.shape(self.wnd_dx)

    
    def convolve_y_dy(self, mean_y, dy, use_wnd_xtm=1, use_wnd_inc_m=0):
        
        """
        
        tranfer  model mean_y and transform dy to reflect changes on mean_x and dx. 
        
        
        Inputs:
        -------------------------------------
        1. mean_y: <array, (ny (nobs),) >: mean model observation 
        2. dy:    <array, (ny, ne)>:  aprior dy (before any scaling)
        3. use_wnd_xtm:<T/F>: if ture wnd_xtm will be applied to dy
        (means that no such transform be made before)
        4. use_wnd_inc_m:<T/F>: if ture  (no adjust made to mean_y before)
        

        Returns
        sdy:<array, (ny, ne)>: sdy if the dy matrix transformed 
        
        """
        
        if (use_wnd_xtm==1):
            # #T1: apply scaling 
            sdy=npy.dot(dy, self.wnd_scale_m)
        
            # #T1: apply temporal correlation 
            
            sdy=npy.dot(sdy, self.wnd_t_cor_m)
            
            # #T3: transfer to posterior of previous steps
            
            tdy=npy.dot(sdy, self.wnd_xtm)
            
        
        elif (use_wnd_xtm==2):
            # #: it is assumed scaling factor has been applied to dx already
            sdy=dy
            tdy=npy.dot(dy, self.wnd_dx)

            tdy=npy.dot(tdy, self.wnd_xtm)
            
            
        
        else:
            # it is assumed scaling transfer for previous observation has be done
            sdy=dy
            tdy=npy.dot(dy, self.xtm)
            
            
            
        if (use_wnd_inc_m==1):
            # #T: nothing has been done to shift mean_y due to changes in the fluxes inside the window
            
            yinc_obs=npy.dot(sdy, self.wnd_inc_m)
            yinc_tcor=npy.dot(sdy, self.tcor_inc_m)
            
        elif (use_wnd_inc_m==2):
            # #T: dx has been applied by a scaling factor 
            wnd_inc_m=self.wnd_mean_x-self.wnd_mean_x0
            yinc_obs=npy.dot(sdy, wnd_inc_m)
            yinc_tcor=npy.dot(sdy, self.tcor_inc_m)
            print 'yinc_tor, max, min:', npy.max(yinc_tcor), npy.min(yinc_tcor)
            
        else:
            # #T: only recent assimilation inside the step have not been included
            
            yinc_obs=npy.dot(sdy, self.inc_m)
            yinc_tcor=npy.dot(sdy, self.tcor_inc_m)

        mean_y=mean_y+yinc_obs+yinc_tcor
        
        return mean_y, tdy
    

    
    def do_assim(self, dy, mean_y, yobs, yerr, full_r=None, use_sparse=True):        
        
        """ assimilating observation and update state vector
        
        Inputs:
        ------------------------------------------------
        1. dy:<array, (nobs(ny), ne)>: ensemble for model 
        observation perturbation. (see Note 1)
        2. mean_y:<array, (ny)>: model observations  (see Note 2)
        3. yobs:  <array, (ny)>: observations 
        4. yerr:  <array, (ny)>: observation error^2 
        5. full_r: <array, (ny, ny)>: observation error covariance 
        6. use_sparse:<T/F>:  if True etkf_cor will be used. 
        
        Returns:
        -----------------------------------------------------
        1. xinc:<array, (nx)>: increment from current observations
        2. sum_inc_m:<array, (ne, nx)>: increment matrix (associated with self.wnd_dx)
        3. inc_m:<array, (ne, nx)>: increment matrix (associated with self.wnd_dx*self.xtm(old))
        4. xtm:<array, (ne, ne)>: transform matrix. 
        
        
        Notes:
        --------------------------------------------------
        1. It is assumed shift y_dy has been called so that 
        dy=H(self.wnd_dx*self.wnd_xtm), and mean_y=H(self.wnd_mean_x)
         
        2. self.wnd_dx=self.wnd_x0*self.wnd_scale_m*self.wnd_tcor_m
           self_wnd_dx0 is a quasi-diagonal matrix           
           
        
        """
        
        
        # S1: get current dx and dy applying transform matrix 
       
        # #c: self.wnd_dx have been applied by scale_m
        
        # #1: applied scale_m (i.e., have been applied scale_m and temporal correlation matrix)
        
        tdx=self.wnd_dx
        
        # #2: applied analysis transform matrix 
        
        tdx=npy.dot(tdx, self.wnd_xtm)
        
        tdy=dy
        
        # S2: solve the gain matrix and transform matrix 
        
        xref=1.0
        
        print 'shape ---:', npy.shape(tdy), npy.shape(tdx), npy.shape(mean_y)
        
        
        if (use_sparse):
            # using sparse matrix
            if (full_r==None):
                full_r=npy.diag(yerr)
            
            full_r=full_r.astype(float)
            
            ett=etc.etkf_cor_cl(self.wnd_mean_x, mean_y, tdx, tdy, full_r, \
                                    xref=xref, xnorm=self.xnorm)
        else:
            r=npy.array(yerr)
            ett=et.etkf_cl(self.wnd_mean_x, mean_y, tdx, tdy, r, \
                               xref=xref, xnorm=self.xnorm)
            
        # S3: get transform matrix 
        
        xtm=ett.get_tm() # (on tdx)
                
        # S4: get increment matrix 
        # #c: current increment matrix associated with tdx
        
        inc_m=ett.get_inc_m(yobs)
        
        # #: increment matrix associated with dx
        xtm_inc_m=npy.dot(self.xtm, inc_m)

        # #c: step increment 
        self.inc_m=self.inc_m+xtm_inc_m

        # #c: window increment 

        wnd_xtm_inc_m=npy.dot(self.wnd_xtm, inc_m)
        
        self.wnd_inc_m=self.wnd_inc_m+wnd_xtm_inc_m
        
        
        
        # S5: get increment over x

        xinc=npy.dot(tdx, inc_m)
        
        self.wnd_xinc=self.wnd_xinc+xinc
        self.xinc=self.xinc+xinc
        
        # S6: update class member
        
        self.xtm=npy.dot(self.xtm, xtm)
        self.wnd_xtm=npy.dot(self.wnd_xtm, xtm)
        self.wnd_mean_x=self.wnd_mean_x+xinc
        
        return xinc, xtm_inc_m, xtm, inc_m
    
    
    def vf_construct_aprior(self, nx, ne, cfg=None):
        
        """
        
        Construction prior and error 
        
        Inputs:
        -------------------------------------------------
        1. nx:<int>: size of x 
        2. ne:<int>: size of ensemble 
        3. cfg:<menu>: configuration  
        
        Returns:
        -------------------------------------------------
        1. x:<array, (nx)>: aprior 
        2. dx:<array, (nx, ne)>: 
        """
        if (cfg==None):
            x=npy.zeros(nx, float)
            dx=npy.zeros([nx, ne], float)
            
            for iy in range(ne):
                if (iy>=nx):
                    break
                else:
                    dx[iy, iy]=1.0
    
        
        return x, dx
    


# <<< TEST >>> 

if (__name__=='__main__'):
    
    yyyy_st,mm_st, dd_st=2009,1,1
    
    max_step=20
    lag_window=5
    
    cl_ens= ens_stat_cl(yyyy_st,\
                            mm_st, \
                            dd_st,\
                            max_step,\
                            lag_window,\
                            mod_bias=[],\
                            mod_bias_err=[],\
                            use_tcor=True,\
                            tcor_len=30.0*24.00,\
                            tcor_factor=1.0)
    nx=5
    ne=4
    x,dx=cl_ens.vf_construct_aprior(nx, ne)
    tau_step=30.0*24.00
    tau_st=0.0
    tau_end=tau_st+tau_step
    
    for istep in range(0, 8):
        nx, ne=cl_ens.add_new_x_to_window(x, \
                                              dx, \
                                              tau_st,\
                                              tau_end,\
                                              enlarge_factor=1.0,\
                                              scale_m=None, \
                                              lag_window=5,\
                                              cl_x2flux=None)
        
    
    
        print 'istep, nx, nx:', istep, nx, ne
        
    
        tau_st=tau_end
        tau_end=tau_st+tau_step
    
    plt.pcolor(cl_ens.xtm)
    plt.colorbar()
    plt.show()
    
    
    print cl_ens.tau_st_lst
    print cl_ens.tau_end_lst
    
    t1=tm.get_tau(2009, 1,28)
    t2=tm.get_tau(2009, 1,30)
    print t2-t1
    
