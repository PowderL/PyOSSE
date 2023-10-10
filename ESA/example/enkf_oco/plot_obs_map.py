from pylab import *
from numpy import *
import ESA.util.otool_ncfile_io as ofb
# import gen_plots as gpl
import fig2_gen_plot as gpl
import ESA.util.time_module as tm
import os
import ESA.util.gp_axis_m as gax

doy_st, doy_end=31,41
doys=range(doy_st, doy_end)
maxv=390
minv=380
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


fpath='./oco_inv_std/'
icount=0
for doy in doys:
    yyyy, mm, dd=tm.doy_to_time_array(doy, 2009)
    sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
    fname=fpath+'daily_obs.'+sdate+".nc"
    if (os.path.isfile(fname)):
        varnames=['lon', 'lat', 'obs', 'mod', 'correct_mod']
        lon, lat, obs, mod, post_mod=ofb.ncf_read(fname, varnames)

        
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
                
            count_clear[ilon, ilat]=count_clear[ilon, ilat]+1
            obs_grid[ilon, ilat]=obs_grid[ilon, ilat]+1.0e6*obs[ix]
            mod_grid[ilon, ilat]=mod_grid[ilon, ilat]+1.0e6*mod[ix]
            post_grid[ilon, ilat]=post_grid[ilon, ilat]+1.0e6*post_mod[ix]

figure(1)
used_idx=where(count_clear>0)
used_idx=squeeze(used_idx)
obs_grid[used_idx]=obs_grid[used_idx]/count_clear[used_idx]
mod_grid[used_idx]=mod_grid[used_idx]/count_clear[used_idx]
post_grid[used_idx]=post_grid[used_idx]/count_clear[used_idx]

obs_grid=where(count_clear==0, -999.0, obs_grid)
mod_grid=where(count_clear==0, -999.0, mod_grid)
post_grid=where(count_clear==0, -999.0, post_grid)



ma_obs = ma.masked_where(transpose(obs_grid) < -990.0, transpose(obs_grid))
ma_mod = ma.masked_where(transpose(mod_grid) < -990.0, transpose(mod_grid))
ma_post = ma.masked_where(transpose(post_grid) < -990.0, transpose(post_grid))




cx=cm.jet
cx.set_over('r')
cx.set_under('w')
subplot(1,1,1)

stitle='Simulated OCO XCO2'

gpl.plot_map(ma_obs, rlon=xx, rlat=yy, maxv=maxv, minv=minv,dv=dv,use_pcolor=1, cbar_vert=0,cmap=cx, title=stitle)
savefig('fig1.png')
figure(2)

subplot(1,1,1)
stitle='Prior Model XCO2 '
gpl.plot_map(ma_mod, rlon=xx, rlat=yy, maxv=maxv, minv=minv,dv=dv,use_pcolor=1, cbar_vert=0,cmap=cx, title=stitle)
savefig('fig2.png')

figure(3)
stitle='Posterior Model XCO2 '
gpl.plot_map(ma_post, rlon=xx, rlat=yy, maxv=maxv, minv=minv,dv=dv,use_pcolor=1, cbar_vert=0, cmap=cx, title=stitle)
savefig('fig3.png')

show()



