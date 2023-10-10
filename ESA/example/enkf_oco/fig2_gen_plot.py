from pylab import *
from numpy import *
from matplotlib import colors
try: 
    from matplotlib.toolkits.basemap import Basemap, shiftgrid
except ImportError:
    from mpl_toolkits.basemap import Basemap, shiftgrid


def plot_map(data, rlon=None, rlat=None,use_pcolor=0,  **keywords):
    
    """ display  the data in a map
    keywords: dict which can include the following keywords  minv=0.0, maxv=0.0, dv=0.0,
    show_map=0, map_proj='cyl', show_lat=0, show_lon=0, lat_0=0.0, lon_0=0.0,
    minlat=-90.0, maxlat=90.0, minlon=0.0, maxlon=360.0, use_log=0, level=level
    """
    
    
    
		
    minv=0.0
    if ('minv' in keywords):
        minv=keywords['minv']
        
    maxv=0.0
        
    if ('maxv' in keywords):
        maxv=keywords['maxv']
        
    dv=0.0
    if ('dv' in keywords):
        dv=keywords['dv']
        
    	
    if ('dv' in keywords):
        dv=keywords['dv']
        
        
        	
    if (maxv>minv):
        rlvl=arange(minv, maxv+dv, dv)
        rlvl[0]=-999.0
        rlvl[size(rlvl)-1]=999.
        
            
    stitle=""
    add_str=""
    if ('title' in keywords):
        add_str=keywords['title']
        add_str=add_str.strip()
        stitle=stitle+' '+add_str
        
    if ('unit' in keywords):
        add_str=keywords['unit']
        add_str=add_str.strip()
        stitle=stitle+'('+ add_str +')'
            
					
    cbar_vert=1
    if ('cbar_vert' in keywords):
        cbar_vert=keywords['cbar_vert']
            
    if (cbar_vert==1):
        orientation='vertical'
    else:
        orientation='horizontal'
            
                	
    show_lat=1
    if ('show_lat' in keywords):
        show_lat=keywords['show_lat']
            
    show_lon=1
    if ('show_lon' in keywords):
        show_lon=keywords['show_lon']
            
    nlon, nlat=shape(data)
            
    vals=array(data)
    
    
    map_proj='cyl'
    if ('map_proj' in keywords):
        map_proj=keywords['map_proj']
        
    lat_0=0.0
    if ('lat_0' in keywords):
        lat_0=keywords['lat_0']
        
    minlat=-90.0
    maxlat=90.0
        
    if ('minlat' in keywords):
        minlat=keywords['minlat']
        
    if ('maxlat' in keywords):
        maxlat=keywords['maxlat']
        
    if (rlat==None):
        dlat=(maxlat-minlat)/nlat
        rlat=arange(minlat,maxlat, dlat)
                        
                        
    minlon=-180
    maxlon=180.0
                        
    if ('minlon' in keywords):
        minlon=keywords['minlon']
        
    if ('maxlon' in keywords):
        maxlon=keywords['maxlon']
			
    if (rlon==None):
        dlon=(maxlon-minlon)/nlon
        rlon=arange(minlon, maxlon, dlon)
                        
    lon_0=0

    do_bdr=0
    if ('do_bdr' in keywords):
        do_bdr=keywords['do_bdr']
    
			
    if ('lon_0' in keywords):
        lon_0=keywords['lon_0']
			
        boundinglat=45

    if ('boundinglat' in keywords):
        boundinglat=keywords['boundinglat']
        
    if (map_proj=='npstere' or map_proj=='spstere'):
        m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
    elif (map_proj=='ortho'):
        m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
    else:
        m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                  urcrnrlon=maxlon, urcrnrlat=maxlat,projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        
    if (rlon[-1]<rlon[0]+360.0):
        rlon=resize(rlon, nlon+1)
        rlon[-1]=rlon[0]+360.0
        vals=squeeze(vals)
        vals=resize(vals, [nlon+1, nlat])
        
    x,y=m(*meshgrid(rlon, rlat))
        
    cmap=cm.Set1
    levels=arange(minv, maxv+1, dv)
    if ('cmap' in keywords):
        cmap=keywords['cmap']
    
    if (maxv>minv):
        if (use_pcolor==1):
            cs0=m.pcolor(x, y, data,  shading='flat', vmin=minv, vmax=maxv, cmap=cmap)
            
        else:
            norm = colors.Normalize(vmin=minv, vmax=maxv)
           
            idx=size(levels)
            print idx
            
            # idx=arange(idx)
            clrs=['w','b', 'g', 'c', 'm', 'y','r','w']
            print clrs
            
            clrs=clrs[0:idx-1]
            print clrs
            if (dv>=1):
                cs0=m.contourf(x, y, transpose(vals), rlvl, levels=levels,vmin=minv, vmax=maxv, colors=clrs,clip_ends = True)
            else:
                cs0=m.contourf(x, y, transpose(vals), rlvl, levels=levels, colors=['b', 'g', 'c','m', 'y','r','w'])
    else:
        if (use_pcolor==1):
            cs0=m.pcolor(x, y, transpose(vals),shading='flat', cmap=cmap)
        else:
            cs0=m.contourf(x, y, transpose(vals))
    
            
    # info(m.drawcoastlines)
    m.drawcoastlines(color='k', linewidth=1.5)
    # m.drawcountries(color=white)
    m.drawmapboundary()
    if (show_lat==1):
        m.drawparallels(arange(minlat, maxlat+30.0, 30.), labels=[1,0, 0,0])
        
    if (show_lon==1):
        m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,1])
        title(stitle, fontsize=12)
    else:
        m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,0])
        title(stitle, fontsize=12)
    
    show_colorbar=1
    if ('cb' in keywords):
        show_colorbar=keywords['cb']
    if (show_colorbar==1):
        colorbar(orientation=orientation, ticks=levels, extend='max')
    if (do_bdr==1):
        lvl=arange(max(vals.flat))
        print shape(x), shape(y), shape(vals)
        # cs2=m.contour(x[0:-1,:], y[0:-1,:], transpose(vals), lvl, colors='k', linewidth=0.5)
    
    return m
def plot_track(rlon, rlat=None,m=None, **keywords):
    show_lat=1
    if ('show_lat' in keywords):
        show_lat=keywords['show_lat']
            
    show_lon=1
    if ('show_lon' in keywords):
        show_lon=keywords['show_lon']

    if (m==None):
        map_proj='cyl'
        if ('map_proj' in keywords):
            map_proj=keywords['map_proj']
        
        lat_0=0.0
        if ('lat_0' in keywords):
            lat_0=keywords['lat_0']
            
        minlat=-90.0
        maxlat=90.0
            
        if ('minlat' in keywords):
            minlat=keywords['minlat']
        
        if ('maxlat' in keywords):
            maxlat=keywords['maxlat']
        
                        
        minlon=-180
        maxlon=180.0
        
        if ('minlon' in keywords):
            minlon=keywords['minlon']
        
        if ('maxlon' in keywords):
            maxlon=keywords['maxlon']
        
                        
        lon_0=0
			
			
        if ('lon_0' in keywords):
            lon_0=keywords['lon_0']
            
        boundinglat=45

        if ('boundinglat' in keywords):
            boundinglat=keywords['boundinglat']
        
        if (map_proj=='npstere' or map_proj=='spstere'):
            m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
        elif (map_proj=='ortho'):
            m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                      urcrnrlon=maxlon, urcrnrlat=maxlat,projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        if (show_lat==1):
            m.drawparallels(arange(minlat, maxlat+30.0, 30.), labels=[1,0, 0,0])
        if (show_lon==1):
            m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,1])
    sgn='+'
    if ('sgn' in keywords):
        sgn=keywords['sgn']
    if ('color' in keywords):
        lcolor=keywords['color']
        if (len(rlon)<100):
            m.plot(rlon, rlat, sgn,color=lcolor, markersize=11)
        else:
            m.plot(rlon, rlat, sgn,color=lcolor, markersize=6)
            
    else:
        m.plot(rlon, rlat, sgn)
    
    m.drawcoastlines(color='w', linewidth=0.5)
    m.drawmapboundary()
    # m.drawcountries(color='k', linewidth=0.5)
    
    if ('title' in keywords):
        stitle=keywords['title']
        title(stitle)
        
    
    return m

def add_text(rlon, rlat,txt, m=None,**keywords):
    show_lat=1
    if ('show_lat' in keywords):
        show_lat=keywords['show_lat']
    rlon=array(rlon)
    rlat=array(rlat)
    
    show_lon=1
    if ('show_lon' in keywords):
        show_lon=keywords['show_lon']

    if (m==None):
        map_proj='cyl'
        if ('map_proj' in keywords):
            map_proj=keywords['map_proj']
        
        lat_0=0.0
        if ('lat_0' in keywords):
            lat_0=keywords['lat_0']
            
        minlat=-90.0
        maxlat=90.0
            
        if ('minlat' in keywords):
            minlat=keywords['minlat']
        
        if ('maxlat' in keywords):
            maxlat=keywords['maxlat']
        
                        
        minlon=-180
        maxlon=180.0
        
        if ('minlon' in keywords):
            minlon=keywords['minlon']
        
        if ('maxlon' in keywords):
            maxlon=keywords['maxlon']
			
                        
        lon_0=0
			
			
        if ('lon_0' in keywords):
            lon_0=keywords['lon_0']
            
        boundinglat=45

        if ('boundinglat' in keywords):
            boundinglat=keywords['boundinglat']
        
        if (map_proj=='npstere' or map_proj=='spstere'):
            m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
        elif (map_proj=='ortho'):
            m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                      urcrnrlon=maxlon, urcrnrlat=maxlat,projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        if (show_lat==1):
            m.drawparallels(arange(minlat, maxlat+30.0, 30.), labels=[1,0, 0,0])
        if (show_lon==1):
            m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,1])
    x,y=m(rlon, rlat)
    for i in range(size(rlon)):
        text(x[i], y[i], txt[i], fontsize=10)
    
        
    if ('title' in keywords):
        stitle=keywords['title']
        title(stitle)
        
    
    return m

def plot_scatter(vals, rlon, rlat=None,m=None, **keywords):
    show_lat=1
    if ('show_lat' in keywords):
        show_lat=keywords['show_lat']
            
    show_lon=1
    if ('show_lon' in keywords):
        show_lon=keywords['show_lon']

    if (m==None):
        map_proj='cyl'
        if ('map_proj' in keywords):
            map_proj=keywords['map_proj']
        
        lat_0=0.0
        if ('lat_0' in keywords):
            lat_0=keywords['lat_0']
            
        minlat=-90.0
        maxlat=90.0
            
        if ('minlat' in keywords):
            minlat=keywords['minlat']
        
        if ('maxlat' in keywords):
            maxlat=keywords['maxlat']
        
                        
        minlon=-180
        maxlon=180.0
        
        if ('minlon' in keywords):
            minlon=keywords['minlon']
        
        if ('maxlon' in keywords):
            maxlon=keywords['maxlon']
			
                        
        lon_0=0
			

        minv=0.0
        if ('minv' in keywords):
            minv=keywords['minv']
        
        maxv=0.0
        
        if ('maxv' in keywords):
            maxv=keywords['maxv']
        dv=0.0
        if ('dv' in keywords):
            dv=keywords['dv']
        

        if ('lon_0' in keywords):
            lon_0=keywords['lon_0']
            
        boundinglat=45

        if ('boundinglat' in keywords):
            boundinglat=keywords['boundinglat']
        
        if (map_proj=='npstere' or map_proj=='spstere'):
            m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
        elif (map_proj=='ortho'):
            m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                      urcrnrlon=maxlon, urcrnrlat=maxlat,projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        if (show_lat==1):
            m.drawparallels(arange(minlat, maxlat+30.0, 30.), labels=[1,0, 0,0])
        if (show_lon==1):
            m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,1])
    sgn='+'
    if ('sgn' in keywords):
        sgn=keywords['sgn']
    x,y=m(rlon, rlat)

    
                  
    # scatter(x, y, c=vals, cmap=cm.jet, vmax=maxv, vmin=minv)
    print max(vals), min(vals)
    s=zeros(size(x))
    s=10.0+s
    h_idx=where(vals>=minv)
    l_idx=where(vals<minv)
    s[h_idx]=30.0
    scatter(x[l_idx],y[l_idx], s=s[l_idx], color='grey', marker='o')
    scatter(x[h_idx],y[h_idx], s=s[h_idx], c=vals[h_idx], cmap=cm.jet, vmax=maxv, vmin=minv, marker='o')
    
    
    m.drawcoastlines(color='b', linewidth=3.0)
    m.drawmapboundary()
    # m.drawcountries(color='k', linewidth=0.5)
    
    if ('title' in keywords):
        stitle=keywords['title']
        title(stitle)
    
    show_colorbar=1
    if ('cb' in keywords):
        show_colorbar=keywords['cb']
    if (show_colorbar==1):
        colorbar() # orientation=orientation)
    return m
    

def plot_track(rlon, rlat=None,m=None, **keywords):
    show_lat=1
    if ('show_lat' in keywords):
        show_lat=keywords['show_lat']
            
    show_lon=1
    if ('show_lon' in keywords):
        show_lon=keywords['show_lon']

    if (m==None):
        map_proj='cyl'
        if ('map_proj' in keywords):
            map_proj=keywords['map_proj']
        
        lat_0=0.0
        if ('lat_0' in keywords):
            lat_0=keywords['lat_0']
            
        minlat=-90.0
        maxlat=90.0
            
        if ('minlat' in keywords):
            minlat=keywords['minlat']
        
        if ('maxlat' in keywords):
            maxlat=keywords['maxlat']
        
                        
        minlon=-180
        maxlon=180.0
        
        if ('minlon' in keywords):
            minlon=keywords['minlon']
        
        if ('maxlon' in keywords):
            maxlon=keywords['maxlon']
			
                        
        lon_0=0
			
			
        if ('lon_0' in keywords):
            lon_0=keywords['lon_0']
            
        boundinglat=45

        if ('boundinglat' in keywords):
            boundinglat=keywords['boundinglat']
        
        if (map_proj=='npstere' or map_proj=='spstere'):
            m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
        elif (map_proj=='ortho'):
            m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                      urcrnrlon=maxlon, urcrnrlat=maxlat,projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        if (show_lat==1):
            m.drawparallels(arange(minlat, maxlat+30.0, 30.), labels=[1,0, 0,0])
        if (show_lon==1):
            m.drawmeridians(arange(minlon,maxlon+60,60.),labels=[0,0,0,1])
    sgn='+'
    if ('sgn' in keywords):
        sgn=keywords['sgn']
    if ('color' in keywords):
        lcolor=keywords['color']
        if (len(rlon)<100):
            m.plot(rlon, rlat, sgn,color=lcolor,  markersize=6)
        else:
            m.plot(rlon, rlat, sgn,color=lcolor, markersize=6)
            
    else:
        m.plot(rlon, rlat, sgn)
    
    m.drawcoastlines(color='w', linewidth=0.5)
    m.drawmapboundary()
    # m.drawcountries(color='k', linewidth=0.5)
    
    if ('title' in keywords):
        stitle=keywords['title']
        title(stitle)
        
    
    return m



if (__name__=="__main__"):
    
    rlon=arange(-180, 180, 10)
    rlat=arange(-90, 90, 10)
    lon_m, lat_m=meshgrid(rlon, rlat)
    data=sin(lon_m*pi/180.0)**2+cos(lat_m*pi/180.0)**2
    data=transpose(data)
    subplot(2,1,1)
    plot_map(data, rlon, rlat, use_pcolor=1, maxv=1.0, minv=-1.0, dv=0.2)
    show()
    
    
    
    
