from pylab import *
from numpy import *
from Scientific.IO.NetCDF import *
import numpy.linalg as nlg
# import co2_emission_std as co2em
import oco_units as ocunit
import gen_plots as gpl
import geos_chem_def as gcdf     # deflaut settings
import numpy as np
import oco_feedback as ofb

windows_length=5
maxstep=23
# inv_path='./oco_inv_std_prior_20110710/'
# inv_path='./oco_inv_wunch_bias_full/'
# inv_path='./oco_inv_redo_acos_full/'
inv_path='./oco_inv_gv_chosen_20130328/'


# output file

errflnm=inv_path+'inv_err_step.nc'

nreg=144
nocean=44
nland=nreg-nocean

xnorm=1.0
nbias=0
# a priori
all_flux0=list()
all_xcor0=list()
all_x_cor0=list()

# post priori
all_flux=list()
all_x_cor=list()
all_xcor=list()


for istep in range(maxstep+1):
    sel_step=istep+windows_length-1
    if (sel_step>maxstep):
        ist=nreg*(sel_step-maxstep)+nbias
        sel_step=maxstep
        
    else:
        ist=nbias
    
    
    sstep=r'%2.2d' % sel_step
    sflnm=inv_path+'std_oco_assim_res.'+sstep
    print sflnm
    
    ncf=NetCDFFile(sflnm+".nc")
    dx=ncf.variables['dx']
    dx=array(dx)

    
    dx_flux=ncf.variables['dx_flux']
    dx_flux=array(dx_flux)
    dx_flux=dx_flux[ist:ist+nreg]
    dx_flux=dx_flux*ocunit.kg_s_to_GtC_Y
    dx_flux=diag(dx_flux)
    
    # prior uncerntainties

    dx0=dx[ist:ist+nreg, ist:ist+nreg]
    dx_cor0=dot(dx0, dx0)
    
    dx0=dot(dx_flux, dx0)
    xcor0=dot(dx0, transpose(dx0))
    xcor0=xcor0/xnorm

    
    # print shape(dx)
    sum_xtm=ncf.variables['sum_xtm']
    
   
    
  

    # posterior uncertainties
    
    dx=dot(dx, sum_xtm)
    dx=dx[ist:ist+nreg, ist:ist+nreg]
    dx_cor=dot(dx, dx)
    
    dx=dot(dx_flux, dx)
    xcor=dot(dx, transpose(dx))
    xcor=xcor/xnorm
    # debug 
    if (istep==-1):
        
        subplot(2,1,1)
        pcolor(xcor0)
        subplot(2,1,2)
        pcolor(xcor)
        show()
    
    print shape(xcor0)
    print shape(xcor)

    all_x_cor0.append(dx_cor0)
    all_x_cor.append(dx_cor)
    
    all_xcor0.append(xcor0)
    all_xcor.append(xcor)


    
    


# flux

mean_x0=ncf.variables['whole_flux0']
mean_x0=array(mean_x0)
print mean_x0[nland:nreg]

# mean_x0=mean_x0[ist:ist+nreg]
mean_x0=mean_x0*ocunit.kg_s_to_GtC_Y
all_flux0=reshape(mean_x0, [-1, nreg])

print mean_x0[nland:nreg]
print shape(all_flux0)
print all_flux0[1,nland:nreg]
print sum(all_flux0[1, nland:])
sx=all_flux0[:, nland:]
sx=sum(sx, axis=1)

print mean(sx)

# mean_x0=mean_x0[ist:ist+nreg]

mean_x0=ncf.variables['whole_x0']
mean_x0=array(mean_x0)
all_x0=reshape(mean_x0, [-1, nreg])


# jjj=raw_input()

mean_x=ncf.variables['whole_flux']
mean_x=array(mean_x)
# mean_x=mean_x[ist:ist+nreg]
mean_x=mean_x*ocunit.kg_s_to_GtC_Y
all_flux=reshape(mean_x, [-1, nreg])



mean_x=ncf.variables['whole_x']
mean_x=array(mean_x)
# mean_x=mean_x[ist:ist+nreg]
all_x=reshape(mean_x, [-1, nreg])


ncf.close()




    
all_xcor0=array(all_xcor0)
all_xcor=array(all_xcor)

all_x_cor0=array(all_x_cor0)
all_x_cor=array(all_x_cor)

all_flux0=array(all_flux0)
all_flux=array(all_flux)

all_x0=array(all_x0)
all_x=array(all_x)



print shape(all_xcor0)

nt, nx, ny=shape(all_xcor0)


xnt=arange(nt)
xnx=arange(nx)
xny=arange(ny)
dimnames=['nt', 'nx', 'ny']
dimvars=[xnt, xnx, xny]
dimtypes=['i', 'i', 'i']

flux0_info=ofb.geos_varinfo('flux0', 'f', ['nt', 'nx'], all_flux0)
flux_info=ofb.geos_varinfo('flux', 'f', ['nt', 'nx'], all_flux)

err0_info=ofb.geos_varinfo('err0', 'f', ['nt', 'nx', 'ny'], all_xcor0)
err_info=ofb.geos_varinfo('err', 'f', ['nt', 'nx', 'ny'], all_xcor)


x0_info=ofb.geos_varinfo('x0', 'f', ['nt', 'nx'], all_x0)
x_info=ofb.geos_varinfo('x', 'f', ['nt', 'nx'], all_x)

x_err0_info=ofb.geos_varinfo('xerr0', 'f', ['nt', 'nx', 'ny'], all_x_cor0)
x_err_info=ofb.geos_varinfo('xerr', 'f', ['nt', 'nx', 'ny'], all_x_cor)


ofb.ncf_write_by_varinfo(errflnm, dimnames, dimtypes, dimvars, [flux0_info, \
                                                                err0_info, \
                                                                x0_info, \
                                                                x_err0_info, \
                                                                flux_info,err_info, x_info,x_err_info])

