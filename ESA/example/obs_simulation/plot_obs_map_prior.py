from pylab import *
from numpy import *
import ESA.util.otool_ncfile_io as ofb
# import gen_plots as gpl
# import fig2_gen_plot as gpl
import ESA.util.gen_plots as gpl

import ESA.util.time_module as tm
import os
import ESA.util.gp_axis_m as gax

doy_st, doy_end=1,32
doys=range(doy_st, doy_end)
maxv=390
minv=379
dv=1.0


xx=arange(-180, 181, 2.5)
yy=arange(-90, 91,2)        
ax_lon=gax.gp_axis_cl('lon', xx)
ax_lat=gax.gp_axis_cl('lat', yy)

ngrd_lat=size(yy)
ngrd_lon=size(xx)

count_clear=zeros([ngrd_lon, ngrd_lat], float)
obs_grid=zeros([ngrd_lon, ngrd_lat], float)
mod_grid=zeros([ngrd_lon, ngrd_lat], float)
post_grid=zeros([ngrd_lon, ngrd_lat], float)
show_vob=False
if (show_vob):
    fpath='./dummy_obs/'
    fprefix='aqua'
else:
    fpath='./sat_obs/'
    fprefix='oco_xco2'

icount=0
for doy in doys:
    yyyy, mm, dd=tm.doy_to_time_array(doy, 2009)
    if (show_vob):
        sdate=r'%4.4dD%3.3d' % (yyyy,doy)
    else:
        sdate=r'%4.4d%2.2d%2.2d' % (yyyy,mm, dd)
    
    fname=fpath+fprefix+'.'+sdate+".nc"
    print fname
    
    if (os.path.isfile(fname)):
        varnames=['lon', 'lat', 'obs', 'obs', 'oerr', 'cnt']
        lon, lat, obs, mod, err, cnt=ofb.ncf_read(fname, varnames)
        
        
        lonp1, lonp2, lonwgt=ax_lon.getwgt(lon)
        latp1, latp2, latwgt=ax_lat.getwgt(lat)

        
        for ix in range(size(lonp1)):
            if (lonwgt[ix]>0.5):
                ilon=lonp1[ix]
            else:
                ilon=lonp2[ix]
                
            if (latwgt[ix]>0.5):
                ilat=latp1[ix]
            else:
                ilat=latp2[ix]
                
            count_clear[ilon, ilat]=count_clear[ilon, ilat]+cnt[ix]
            
            obs_grid[ilon, ilat]=obs_grid[ilon, ilat]+1.0e6*obs[ix]*cnt[ix]
            mod_grid[ilon, ilat]=mod_grid[ilon, ilat]+1.0e6*mod[ix]*cnt[ix]
            

figure(1)
used_idx=where(count_clear>0)
used_idx=squeeze(used_idx)
obs_grid[used_idx]=obs_grid[used_idx]/count_clear[used_idx]
mod_grid[used_idx]=mod_grid[used_idx]/count_clear[used_idx] # -1.5
# post_grid[used_idx]=post_grid[used_idx]/count_clear[used_idx]

obs_grid=where(count_clear==0, -999.0, obs_grid)
mod_grid=where(count_clear==0, -999.0, mod_grid)
# post_grid=where(count_clear==0, -999.0, post_grid)



ma_obs = ma.masked_where(transpose(obs_grid) < -990.0, transpose(obs_grid))
ma_mod = ma.masked_where(transpose(mod_grid) < -990.0, transpose(mod_grid))
# ma_post = ma.masked_where(transpose(post_grid) < -990.0, transpose(post_grid))




cx=cm.jet
cx.set_over('r')
cx.set_under('w')
# figure(1)

subplot(2,1,1)

stitle='Observation number'

gpl.plot_map(count_clear, rlon=xx, rlat=yy, use_pcolor=1, linecolor='w', title=stitle, cb_vert=0)



subplot(2,1,2)
stitle='OCO XCO2 in Jan, 2009'

gpl.plot_map(ma_obs.transpose(), rlon=xx, rlat=yy, maxv=maxv, minv=minv,dv=dv,\
                 use_pcolor=1, cmap=cx, title=stitle, cb_vert=0)



# stitle='GEOS-Chem (mod)'

# gpl.plot_map(ma_mod, rlon=xx, rlat=yy, maxv=maxv, minv=minv,dv=dv,use_pcolor=1, cmap=cx, title=stitle)

# figure(2)

savefig('oco_200901.png')

show()



