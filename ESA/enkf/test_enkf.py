import numpy as npy
import etkf_cor as ec
import etkf_half as et
import numpy.random as rnd


def get_obs(x, nobs):
    dims=npy.shape(x)
    nx=dims[0]
    
    h=npy.zeros([nobs, nx], float)
    
    for iobs in range(nobs):
        if (iobs<nx):
            h[iobs, iobs]=2.0
    
    #     h[iobs, :]=0.5
    
    obs=npy.dot(h, x)
    return obs


if (__name__=='__main__'):
    
    test_etkf_cor=0  # 0: test etkf_half
    use_rnd=False
    nsel=50
    
    nobs=100
    nx=100
    
    
    # true 
    
    xtrue=npy.arange(1, nx+1)
    
    xprior=xtrue+rnd.normal(scale=2, size=nx)
    
    xtrue[4]=10.0
    
    # to observation space
    ytrue=get_obs(xtrue, nobs)
    yprior=get_obs(xprior, nobs)
    
    
    # ensemble 
    
    dx=npy.ones(nx, float)
    dx=npy.diag(dx)
    
    if (use_rnd):
        dx=list()

        for  i in range(nsel):
            sel_x=rnd.normal(scale=2, size=nx)
            dx.append(sel_x)
    
        dx=npy.array(dx)
        dx=dx/npy.sqrt(nsel)
        dx=npy.transpose(dx)
        print 'dx', dx
        

                         
    dy=get_obs(dx, nobs)
    print npy.shape(dx), npy.shape(dy)
    
        
    
    # oerr
    r=npy.ones(nobs, float)
    r=npy.diag(r)
    
    r=0.2*r
    
    if (test_etkf_cor==1):
        
        etk=ec.etkf_cor_cl(xprior, yprior, dx, dy, r)
    else:
        r=npy.diag(r)
        
        etk=et.etkf_cl(xprior, yprior, dx, dy, r)
        
    tm=etk.get_tm()
    inc_m=etk.get_inc_m(ytrue)
    xinc=npy.dot(dx, inc_m)
    xpost=xprior+xinc
    print 'x:'
    print ' '

    print 'true', xtrue
    print 'prior', xprior
    print 'post', xpost
    print ' '

    xpost2=etk.get_posterior(ytrue)
    
    
    print 'post2', xpost2
    
    xinc3=etk.get_x_inc(dx, inc_m)
    xpost3=xprior+xinc3
    print 'post3', xpost3
    
    
    
    print 'y:'
    print ' '

    yinc=npy.dot(dy, inc_m)
    ypost=yprior+yinc
   
    print 'ytrue:', ytrue
    print 'yprior:', yprior
    
    print 'ypost:', ypost
    yinc2=etk.get_y_inc(dy, inc_m)
    ypost2=yprior+yinc2
    print 'ypost2:', ypost2
    
    print ' '
   
    
    
#     print 'tm:'
#     print tm
    
    

    
