"""
    Solve Ensemble transform matrix by using SVD approach
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. etkf_cl:class for solving ETKF data assimilation equation 
    
"""


import numpy as npy
import numpy.linalg as nlg
import sys

# import matrix_tool as mtool

class etkf_cl:
    """ class for Ensemble Transform Kalman Filter solver
    See Feng et al., ACP 2009

    Members:
    -----------------------------------------------------
    
    # control variables 
    1. nx:<int>: size of control variable 
    2. ne:<int>: size of perturbation ensemble 
    3. x:<array, (nx, ne)>: perturbation ensemble 
    4. mean_x:<array, (nx)>: mean value of control variable 
    5. xnorm:<float>: normalized factor 
    
    # model observation 
    6.  y:<array, (nobs, ne)>: ensemble for model observation perturbation 
    7.  ny:<int>: ny=nobs
    8.  mean_y:<array, (nobs)>: model observations for self.mean_x
    9.  r:<array, (nobs, )>: observation err covariance 
    10. sqr:<array, (nobs,)>: sqr=sqrt(r) 
    11. scyt:<array, (ne, nobs): scaled y^t  scyt=transpose(self.y)/self.sqr 
    12. uy:<array, (ne,ne)>: U matrix of svd(scyt) 
    13. vy:<array, (ny, ny)>:V matrix of svd(scyt)
    12. wy:<array, (nwy)>: eigenvalues of scyt. newy=min(ne, ny)
    
    13. md:<array, (nwy)>: md=1.0/(1.0+wy*wy)
    14  n_wy:<int>: size of non-zeros wy 
    15. uw: <array, (ne, n_wy)>: reduced U matrix 
    

    Functions
    ---------------------------------------------------------------
    1.  __init__: initialization
    2.  reset: re-do initialization
    3.  get_posterior: get posterior 
    4.  get_tm: get transform matrix
    5.  get_k: get gain matrix 
    6.  get_inc_m: get increment matrix 
    7.  get_x_inc: get analysis increment
    8.  get_y_inc: get analysis increment in observation space
    9.  do_y_transform: transform maxtrix in y space
    
      
    
    """
    def __init__(self, mean_x, mean_y, x, \
                     y, r, xref=1, xnorm=1.0, do_debug=False):
        """ initialize the class 
        
        Inputs
        ---------------------------------------------------------
        1. mean_x:<array, (nx)>: mean values of control variables x
        2. mean_y:<array, (ny)>: mean model observation values.  
        3. x: <array, (nx, ne)>:  ensemble of x perturbations 
        4. y: <array, (nobs, ne)>:  ensemble of model forecast perturbations. ny=nobs 
        5. r: <array, (nobs)>: observation error covariances given in 1 dimension
        6. xref:<float>: value used to scale mean_x 
        7. xnorm:<float>: value used to scale ensemble purbations 
        
           
        Notes:  
        ----------------------------------------
        1. self.x=(x-mean(x))/sqrt(xnorm) 
        2. self.y=(y-mean(y))/sqrt(xnorm) 
       
        """
        
        self.nx, self.ne=npy.shape(x)
        self.x=npy.array(x)
        
        print 'etkf: ns, ne', self.nx, self.ne
        
        self.mean_x=npy.array(mean_x)
        print 'etkf: shape of mean_x', npy.shape(self.mean_x)
        
        
        self.xnorm=xnorm
        
        # scaling x ensemble 
        
        self.x=self.x/npy.sqrt(self.xnorm)
        
        # #c: error covariance of x 
        
        dev=npy.dot(self.x, npy.transpose(self.x))
        
        self.y=npy.array(y)
        self.ny=npy.size(y[:,0])  # ---> nobs 
        self.mean_y=npy.array(mean_y)
        
        # #c: scaling y ensemble 
        
        
        self.y=self.y/npy.sqrt(self.xnorm)
        self.yy=None
        self.r=r
        # #c: square root of error covariance
        
        self.sqr=npy.sqrt(r) 
        
        # #c:scaling y ensemble:  y^t=r^-0.5*y^t 
        self.scyt=npy.transpose(self.y)/self.sqr 
        
        # S2:  SVD of scaled y^t
        
        # #c: decomposition of tranpose(scaled_y)
        self.uy, self.wy, self.vyt=nlg.svd(self.scyt) 
        if (do_debug):
            
            scx=self.x[:,:]/xref
            dev=dot(scx, transpose(scx))
            plt.pcolor(dev)
            plt.colorbar()
            
            print 'etkf: shape(y)',  shape(self.y)
            print 'etkf: y[0:8,-6]',  self.y[0:8,-6]
            print 'etkf: r',  self.r[0:8]
            print 'etkf: wy', self.wy[0:8]

        # S3: Calculate   (1.0+wy*wy)^-1
        
        # #c:in calculation of  K , we need 1.0/(1.0+w*wt), which is shape of  [nobs, nobs]
        # #c:in calculation of transform matrix, we need 1.0/(1.0+wt*w) of [ne, ne]
        # #c: md is 1d array of size of min(ny, ne)
          
        self.md=1.0/(1.0+self.wy*self.wy) 
        
        
        idx=npy.where(self.wy>0.0)
        idx=npy.squeeze(idx)
        # #c: n_wy is the size of wy >0.0
        
        self.n_wy=npy.size(idx)
        
        
        # #c: should have the shape of [ne, nobs] but  now only have [ne, self.n_wy]

        self.uw=self.uy[:, 0:self.n_wy]*self.wy[0:self.n_wy]  

    def reset(self, mean_x, mean_y, x, y, r, xref=1.0, xnorm=1.0, do_debug=False):
        
        """ re set (reinitialize) the class 
        Input:  
        ------------------------------------------------------------
        1. mean_x:<array, (nx)>: mean values of control variables x
        2. mean_y:<array, (nobs)>: mean model observation values 
        3. x: <array, (nx, ne)>:  ensemble of x perturbations 
        4. y: <array, (nobs, ne)>:  ensemble of model forecast perturbations 
        5. r: <array, (nobs)>: observation error covariances given in 1 dimension
        6. xref:<float>: value used to scale mean_x 
        7. xnorm:<float>: valuse used to scale en

        """
        
        # S1: test whether decomposition of y^t converge
        
        dy=npy.array(y)
        dy=dy/npy.sqrt(xnorm)
        sqr=npy.sqrt(r)
        scyt=npy.transpose(dy)/sqr # the scaled y=r^-0.5*y ==>scaled yt=yt*r^-0.5
        
        try: 
            uy, wy, vyt=nlg.svd(scyt)
        except nlg.LinAlgError:
            return -1
        
        # S2: re-set x 

        self.nx, self.ne=npy.shape(x)
        self.x=npy.array(x)
        self.xnorm=xnorm
        # print self.nx
        self.x=self.x/npy.sqrt(self.xnorm)
        
        # S3: reset y

        self.y=npy.array(y)
        self.ny=npy.size(y[:,0])  # ---> nobs 
        self.mean_y=npy.array(mean_y)
        self.y=self.y/npy.sqrt(self.xnorm)
        self.yy=None
        self.r=r
        self.sqr=npy.sqrt(r)
       
        
        self.scyt=npy.transpose(self.y)/self.sqr # the scaled y=r^-0.5*y ==>scaled yt=yt*r^-0.5
        self.uy, self.wy, self.vyt=nlg.svd(self.scyt) # decomposition of tranpose(scaled_y)
    
        # S4: calculate md 
        
        
        self.md=1.0/(1.0+self.wy*self.wy) #  1d array of min[ny, ne] 
        
        # #c: in calculation of  K , we need 1.0/(1.0+w*wt)  of [nobs, nobs]
        # #c: in calculation of transform matrix, we need 1.0/(1.0+wt*w) of [ne, ne]
        
        if (do_debug):
            # plot(self.wy)
            # plot(self.scyt[0,0:50])
            # plot(self.scyt[1,0:50])
            scx=self.x[:,:]/xref
            dev=dot(scx, transpose(scx))
            plt.pcolor(dev)
            plt.colorbar()
            
        # #c: set self.n_wy 
        
        idx=npy.where(self.wy>0.0)
        idx=npy.squeeze(idx)
    
        self.n_wy=npy.size(idx)
        self.uw=self.uy[:, 0:self.n_wy]*self.wy[0:self.n_wy]  # should have the shape of [ne, nobs] but    now only have [ne, self.n_wy
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
        # S1: get incream matrix 
        inc_m=self.get_inc_m(yobs)

        # S2: return posterior flux 
        x=self.mean_x+npy.dot(self.x, npy.transpose(inc_m))
        
        return x
    
    def get_tm(self):
        """ the analysis increment on the variables  
        
        Returns:
        1. tm :<array, (ne, ne)>: matrix for posterior ensemble: x(f)=x(a)*tm 
        
        Notes: 
        -----------------------------------
        1. tm=U*(1+W*WT)^-0.5)*UT
        
        """    

        inv_s=npy.ones(self.ne, float)
        print 'get_tm: m_wy', self.n_wy
        inv_s[0:self.n_wy]=npy.sqrt(self.md[0:self.n_wy])
        print 'get_tm:  inv_s', inv_s[0:6]
        tm=self.uy*inv_s
        tm=npy.dot(tm, npy.transpose(self.uy))
        
        # ['tm 0:20']
        
        print '-'*8+"etkf: xtm"+'-'*8
        print npy.array2string(tm[-5:-1,-5:-1],precision=5, suppress_small=True)
        
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
        1. inc_m=U*W*(I+WT*W)^-1*VT*R^-0.5*(Yobs-mean_Y)
        """
    
        # v_nx,v_ny=shape(self.vy)
        print npy.shape(yobs), npy.shape(self.mean_y)
        
        dy=yobs-self.mean_y # innovation
        print 'get_inc_m: yobs-mean_y', dy[0:6], max(dy), min(dy)
        
        # print dy
        
        dy=dy/self.sqr
        
        inc_m=npy.dot(self.vyt, npy.transpose(dy))
        
        inc_m=npy.squeeze(inc_m)
        inc_m=npy.transpose(inc_m)
        inv_s=npy.ones(self.ny, float) # filled with one 
        inv_s[0:self.n_wy]=self.md[0:self.n_wy]
        inc_m=inv_s*inc_m
        uw_nx, uw_ny=npy.shape(self.uw) 
        rny=self.ny
        
        if (uw_ny<self.ny):
            rny=uw_ny
        
        inc_m=npy.dot(self.uw[:,0:rny], npy.transpose(inc_m[0:rny]))
        
        return inc_m
        

    def get_k(self):
        
        """ calculate analysis gain k matrix 
        
        Returns: 
        ----------------------------------------------------
        1. k matrix:<array, (nx, nobs)>:  gain matrix
        
        
        2. Notes:
        ------------------------------------------------------
        K=X'*U*W*(I+WT*W)^-1*VT*R^-0.5
        """

        inv_s=ones(self.ny, float) # filled with one 

        inv_s[0:self.n_wy]=self.md[0:self.n_wy]
        
        uw_nx, uw_ny=npy.shape(self.uw)
        
        rny=self.ny
        
        if (uw_ny<rny):
            rny=uw_ny
        
        y_inv=self.uw[:,0:rny]*inv_s[0:rny]
        y_inv_v=npy.dot(y_inv, self.vyt[0:rny, :]) #dimension of [ne, ny]
        k=npy.dot(self.x, y_inv_v) # dimension [nx, ne] X [ne, ny]
        k=k/self.sqr
        
        return k

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
        # mean_y=mean(y_in, axis=1)
        if (mean_y==None):
            dy=npy.array(y_in[:,:])
            
        else:
            dy=y_in[:,:] -mean_y[:, newaxis]
        
        
        dy=dy/npy.sqrt(self.xnorm)
        yinc=npy.dot(dy, npy.transpose(inc_m))
        
        return yinc
    
    def do_y_transform(self, y_in, tm, mean_y=None):
        
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
    import numpy.random as rnd
    
    
    
    
        
        
        
    
        
    
        
        
