
"""
    Solve Ensemble transform matrix by using  approach LU solver
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.09.11
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. etkf_cl:class for solving ETKF data assimilation equation 
    
    
"""

import numpy as npy
import numpy.linalg as nlg
import sys
from scipy.sparse import csr_matrix
from scipy.sparse.linalg.dsolve import linsolve


class etkf_cor_cl:
    """ class for Ensemble Transform Kalman Filter solver for correlated observations
    
    See Feng et al., ACP 2009; Palmer et al., 2011.
    

    Members:
    -----------------------------------------------------
    
    # control variables 
    1. nx:<int>: size of control variable 
    2. ne:<int>: size of perturbation ensemble 
    3. x:<array, (nx, ne)>: perturbation ensemble 
    4. x_mean:<array, (nx)>: mean value of control variable 
    5. xnorm:<float>: normalized factor 
    
    # model observation 
    6.  y:<array, (nobs, ne)>: ensemble for model observation perturbation 
    7.  ny:<int>: ny=nobs
    8.  y_mean:<array, (nobs)>: model observations for self.mean_x
    9.  sp_y: <array, (nobs, ne)>: sparse matrix of self.y

    # observation errors 

    9.  r:<array, (nobs, nobs)>: observation err covariance 
    10. sp_r:<array, (nobs,)>: sparse matrix of r

    11. sp_b:<array, (nobs, nobs)>: background error covariance in obsevation space 
        sp_b=sp_y*sp_y^T
    12. sp_r:<array, (nobs, nobs)>: sp_r=sp_r+
    
    

    Functions
    ---------------------------------------------------------------
    1.  __init__: initialization
    2.  reset: re-do initialization
    3.  get_posterior: get posterior 
    4.  get_tm: get transform matrix
    6.  get_inc_m: get increment matrix 
    7.  get_x_inc: get analysis increment
    8.  get_y_inc: get analysis increment in observation space
    9.  do_y_transform: transform maxtrix in y space
    
      
    
    """

   
    def __init__(self, x_mean, y_mean, x, y, r, xref=1, xnorm=1.0):
        
        """
        Initialization
        
        Inputs:
        ---------------------------------------------------------
        1. x_mean:<array, (nx)>: mean values of control variables x
        2. Y_mean:<array, (ny)>: mean model observation values.  
        3. x: <array, (nx, ne)>:  ensemble of x perturbations 
        4. y: <array, (nobs, ne)>:  ensemble of model forecast perturbations. ny=nobs 
        5. r: <array, (nobs, nobs)>: observation error
        6. xref:<float>: value used to scale mean_x 
        7. xnorm:<float>: value used to scale ensemble purbations 

        
        """
        # S1: set  x  
        self.nx, self.ne=npy.shape(x)
        self.x=npy.array(x)
        print 'etkf: nx, ne', self.nx, self.ne
        self.x_mean=npy.array(x_mean)
        print 'etkf: x_mean', self.x_mean[0:6]
        
        
        
        
        self.xnorm=xnorm
        self.x=self.x/npy.sqrt(self.xnorm)
        
        # S2: set y 
        self.y=npy.array(y)
        self.y=self.y/npy.sqrt(self.xnorm)
        self.sp_y=csr_matrix(self.y)
        
        self.ny=npy.size(y[:,0])  # ---> nobs 
        self.y_mean=npy.array(y_mean)
        self.yy=None

        # S3: observation error 
        self.r=r
        self.sp_r=csr_matrix(r)

        # #c: Background error covariance in obsevation space
        self.sp_b=self.sp_y*self.sp_y.T
        # #c: S matrix: R+DY*DY^T
        self.sp_s=self.sp_b+self.sp_r
        

    
        
    def reset(self, x_mean, y_mean, x, y, r, xref=1.0, xnorm=1.0, do_debug=False):
        """

        Re-Initialization
        
        Inputs:
        ---------------------------------------------------------
        1. x_mean:<array, (nx)>: mean values of control variables x
        2. Y_mean:<array, (ny)>: mean model observation values.  
        3. x: <array, (nx, ne)>:  ensemble of x perturbations 
        4. y: <array, (nobs, ne)>:  ensemble of model forecast perturbations. ny=nobs 
        5. r: <array, (nobs, nobs)>: observation error
        6. xref:<float>: value used to scale mean_x 
        7. xnorm:<float>: value used to scale ensemble purbations 
        
        """
        
        

        # S1: X 
        self.nx, self.ne=npy.shape(x)
        self.x=npy.array(x)
        print 'etkf: nx, ne', self.nx, self.ne
        self.x_mean=npy.array(x_mean)
        print 'etkf: x_mean', self.x_mean[0:6]

        
        self.xnorm=xnorm
        self.x=self.x[:,:]
        self.x=self.x/npy.sqrt(self.xnorm)
        
        # S2: Y

        max_y, min_y=npy.max(y.flat), npy.min(y.flat)
        max_y=npy.max([abs(max_y), abs(min_y)])
        ulimt=1.0e-8*max_y
        new_y=npy.where(abs(y)<ulimt, 0, y)
        
        self.y=new_y/npy.sqrt(self.xnorm)
        self.sp_y=csr_matrix(self.y)
        
        self.ny=npy.size(y[:,0])  # ---> nobs 
        self.y_mean=npy.array(y_mean)
        
        # S3: observation error 
        self.r=r
        self.sp_r=csr_matrix(r)
        
        # #c: Background error covariance in obsevation space
        self.sp_b=self.sp_y*self.sp_y.T

        # #c: S matrix: R+DY*DY^T
        
        self.sp_s=self.sp_b+self.sp_r
        
        
        
        return 0
    
    
    def get_posterior(self,yobs):
        """ calculate posterior from assimilation of observations 
        Inputs:
        -------------------------------------------------
        1. yobs:<array, nobs>: observations
        
        Returns:
        ---------------------------------------
        1. postprior:<array, nx>: 
        
        Notes:
        -----------------------------------------
        1. xa=xf+k*(yobs-f(x))=xf+dx*inc_m
        
        """

        inc_m=self.get_inc_m(yobs)
        
        x=self.x_mean+npy.dot(self.x, npy.transpose(inc_m))
        return x
    
    def get_tm(self):
        """ the analysis increment on the variables  
        
        Returns:
        ----------------------------------------------
        1. tm :<array, (ne, ne)>: matrix for posterior ensemble: x(f)=x(a)*tm 
        
        Notes: 
        -----------------------------------
        1. tm*tm=[1.0+(y^t*R^-1*y) ]
        
        
        """
        
        
        tm=list()
        
        asf=linsolve.factorized(self.sp_r)
        
        print 'get_tm:', npy.shape(self.r), npy.shape(self.y), self.ne
        
        # S1: solve R*X=Y
        # #c: X=R^-1*Y
        
        for ie in range(self.ne):
            #  t1=linsolve.spsolve(self.sp_r, self.y[:,ie])
            # print 'ie', ie
            
            # print npy.shape(self.y[:,ie])
            tx=self.y[:,ie]
            tx=npy.squeeze(tx)
            # np.save('test_z', tx)
            t1=asf(tx) # self.y[:,ie])
            tm.append(t1)
        tm=npy.array(tm)
        
        print 'get_tm:shape(tm) in R*X=Y:', npy.shape(tm)
        
        # S2: solve TM*TM=(1.0+R*Y*Y^T*R^-1)^-0.5
        

        
        tm=npy.transpose(tm)
        sp_tm=csr_matrix(tm)
        s1=self.sp_y.T*sp_tm
        s1=0.5*(s1+s1.T)
        
        s1=s1.todense()
        s1=npy.array(s1)
        
        s1=npy.identity(self.ne)+s1
        u, w, vt=nlg.svd(s1)
        # tm=u/npy.sqrt(1.0+w)
        tm=u/npy.sqrt(w)
        
        tm=npy.dot(tm, npy.transpose(u))
        
        print 'get_tm: w=', w[0:5]
        print 'get_tm shape(tm):', npy.shape(tm)
        
        
        return tm
    
        
    def get_inc_m(self,yobs):
        
        """ get the analysis increment matrix 
        
        
        Inputs:
        -------------------------------------------------
        1. yobs:<array, (nobs)>: the observations
       
        Returns:     
        ------------------------------------------------
        1. inc_m :<array, (ne,)>:increase matrix 
        
        
        Notes:
        --------------------------------------------------------
        1. inc_m=y^T*(R+Y*Y^t)^-0.5*dy
        
        """
    
        # v_nx,v_ny=npy.shape(self.vy)
        
        dy=yobs-self.y_mean # innovation
        
        # S1: solve the linear equation:  (R+Y*Y^t) X=yobs-self.y_mean
        # so that X=(R+Y*Y^t)^-0.5*dy 
                
        inc_m=linsolve.spsolve(self.sp_s, dy)
        inc_m=npy.array(inc_m)
        
        # S2: inc_m=y^T*X 
        
        inc_m=npy.dot(self.y.T, inc_m)
        
        return inc_m
        # self.ana_npy.array=npy.transpose(ana_npy.array)
        
    def get_x_inc(self, inc_m, x_in=None, mean_x_in=None):
        """ the analysis increment of control variables  
        Inputs:
        ------------------------------------------
        1. inc_m:<array, (ne)>: increment matrix 
        2. x_in:<array, (nx, ne)>: ensemble perturbations
        3. mean_x_in:<array, (nx)>: mean value
        
        Returns:
        ---------------------------
        1. xinc: <array, nx>: the x increments      
        """    
        
        if (x_in==None):
            xinc=npy.dot(self.x, npy.transpose(inc_m))
        else:
            if (mean_x_in==None):
                dx=npy.array(x_in)
            else:
                dx=x_in-mean_x_in[:,newaxis]
            
            dx=dx/npy.sqrt(self.xnorm)
            xinc=npy.dot(dx, npy.transpose(inc_m))
            
        return xinc
    def get_y_inc(self, y_in, inc_m, mean_y=None):
        """ the analysis increment on the variables  
        
        Inputs:
        --------------------------------------------
        1. y_in:<array, (new_ny, ne)>: the 'new' observation perturbation 
        2. inc_m :<array, (ne)>: increment matrix 
        
        Returns:
        -------------------------------------------------------
        1. yinc: <array, new_ny>: y increments from xinc     
        
        """    
        
        if (mean_y==None):
            dy=npy.array(y_in[:,:])
            # -mean_y[:, newaxis]
        else:
            dy=y_in[:,:] -mean_y[:, newaxis]
        
        dy=dy/npy.sqrt(self.xnorm)
        yinc=npy.dot(dy, npy.transpose(inc_m))
        return yinc
    
    def get_y_transform(self, y_in, tm, mean_y=None):
        """ multiply y_in by tm 
        
        
        Inputs:
        --------------------------------------------------------------
        1. y_in:<array>, (new_ny, ne)>: 'new' observation perturbations 
        2. tm  :<array, (ne, ne)>: Transformation matrix
        
        Returns:
        ---------------------------------------------------------
        1. yt:<array, (new_ny, ne)>: yt=y_in*tm
        
        """    
        
        
        if (mean_y==None):
            dy=npy.array(y_in[:,:])
            # -mean_y[:, newaxis]
        else:
            dy=y_in[:,:] -mean_y[:, newaxis]
        
        
        # mean_y=mean(y_in, axis=1)
        # dy=y_in[:,:]-mean_y[:, newaxis]
        yt=npy.dot(dy, tm)
        return yt
        
        
     
        
if (__name__=='__main__'):
    print 'that is a test'
    
    
    
        
        
        
    
        
    
        
        
