import ESA.util.otool_ncfile_io as ncfio
import ESA.util.geo_constant as gc


import get_region_index as reg_idx
from numpy import *
from pylab import *

step=3
nstep=step+1

sel_reg='nat'

kg_s_to_GtC_Y=gc.mc*3600.0*365*24.0/(1.0e12*gc.mco2)

sstep=r'%2.2d' % step
datapath='./oco_inv/'
filename='oco_assim'+'.'+sstep+'.nc'
filename=datapath+filename
varnames=['whole_flux', 'whole_flux0', 'sum_xtm']

whole_flux, whole_flux0, xtm=ncfio.ncf_read(filename, varnames)
pcolor(xtm)
show()

print shape(whole_flux)
print 'total_flux:', kg_s_to_GtC_Y*sum(whole_flux0[0:144])
whole_flux=reshape(whole_flux, [-1, 144])
whole_flux0=reshape(whole_flux0, [-1, 144])

print 'total_flux2:', kg_s_to_GtC_Y*sum(whole_flux0[0, 0:144])

sel_idx, stitle=reg_idx.get_region_index(sel_reg)

nt3_reg=23
t3_sep=zeros(nt3_reg)

xnorm=1.0

last_nreg=0
for ireg in range(nt3_reg):
    if (ireg==0):
        t3_sep[ireg]=0
        last_nreg=0
    elif (ireg<=11):
        t3_sep[ireg]=last_nreg+9
        last_nreg=last_nreg+9 
    else:
        t3_sep[ireg]=last_nreg+4
        last_nreg=last_nreg+4

print  t3_sep
nreg=144

cur_sep=arange(nreg)

reg_idx_table=searchsorted(t3_sep, cur_sep)
print reg_idx_table

icount=0
sum_reg_idx=list()
for ireg in sel_idx:
    sel_reg_idx=where(reg_idx_table==ireg)
    # print ireg
    
    sel_reg_idx=squeeze(sel_reg_idx)

    if (size(sel_reg_idx)==1):
        print ireg
        sum_reg_idx.append(sel_reg_idx)
    else:
        for usd_reg in sel_reg_idx:
            sum_reg_idx.append(usd_reg)


sum_reg_idx=array(sum_reg_idx)

print 'sum_reg_idx', sum_reg_idx

sel_flux=whole_flux[:, sum_reg_idx]
sel_flux0=whole_flux0[:, sum_reg_idx]

sel_flux=sum(sel_flux, axis=1)
sel_flux0=sum(sel_flux0, axis=1)
sel_flux=kg_s_to_GtC_Y*sel_flux
sel_flux0=kg_s_to_GtC_Y*sel_flux0


xstep=arange(nstep)
title(stitle)
print shape(sel_flux0)
print shape(xstep)
xlabel('Months in 2009')

plot(xstep, sel_flux, 'r')
plot(xstep, sel_flux0, 'b')
plot(xstep, sel_flux0/1.5, 'k')
legend(['post', 'prior', 'true'])
savefig(sel_reg+'.png')
show()







 
