from pylab import *
from numpy import *
import ESA.util.otool_ncfile_io as ofb
import ESA.util.gen_plots as gpl
import ESA.util.time_module as tm
import os

doy_st, doy_end=1,121
doys=range(doy_st, doy_end+1)
sdate_st=r'(Doy: %3.3d to %3.3d)' % (doy_st, doy_end)
fpath='./oco_inv/'

# fpath='./oco_inv_gv_geos_daily_bkup20120408/'



yyyy=2009
ndoy=doy_end-doy_st+1
date_lst=list()
obs_lst=list()
mod_lst=list()
post_lst=list()
nobs_lst=list()

subplot(2,1,1)

doy=doy_st-1
for idoy in range(ndoy):
    doy=doy+1
    if (doy>365):
        doy=doy-365
        yyyy=yyyy+1
    
    yyyy, mm, dd=tm.doy_to_time_array(doy, yyyy)
    sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
    fname=fpath+'daily_obs.'+sdate+".nc"
    print fname
    if (os.path.isfile(fname)):
        varnames=['lon', 'lat', 'obs', 'mod', 'correct_mod']
        lon, lat, obs, mod, post_mod=ofb.ncf_read(fname, varnames)
        
        used_idx=where(obs>0)
        
        used_idx=where((lon>=-180) & (lon<=180) & (lat>=-90) & (lat<90))
        
        # used_idx=squeeze(used_idx)
        if (size(used_idx)>=1):
            
            obs=1.0e6*obs[used_idx]
            mod=1.0e6*mod[used_idx]
            post=1.0e6*post_mod[used_idx]

            
            nx=size(obs)
            x=zeros(nx, float)
            x[:]=idoy
            plot(x, obs, '+', color='b')
            plot(x, mod, '.', color='k')
            date_lst.append(idoy)
            obs_lst.append(mean(obs))
            mod_lst.append(mean(mod))
            post_lst.append(mean(post))
            nobs_lst.append(size(obs))
            print 'nobs:', size(lon[used_idx])
            
            
            
            
            
                            
            
date_lst=array(date_lst)
obs_lst=array(obs_lst)
post_lst=array(post_lst)
nobs_lst=array(nobs_lst)

lobs=plot(date_lst, obs_lst, 'cyan', linewidth=2)
lmod=plot(date_lst, mod_lst, 'r', linewidth=2)
lpost=plot(date_lst, post_lst, 'yellow', linewidth=2)


ylabel('ppm')
xlabel('Days')
ylim([370, 400])


legend([lobs, lmod, lpost], ['obs', 'model', 'Posteriori'])
title('XCO2')
# xlim([0, 360])
title('global observations')

subplot(2,1,2)
plot(date_lst, nobs_lst)
ndate=len(date_lst)

ylabel('Number of Obs')

savefig('global_obs.png')
savefig('global_obs.pdf')
fl1=open('count_obs.dat', 'w')
for idate in range(ndate):
    count_ln=r'%3.3d, %4d' % (date_lst[idate], nobs_lst[idate])
    fl1.write(count_ln+'\n')
fl1.close()

show()



