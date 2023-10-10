"""
Functions for making map plots etc 
1) plot_map: contour plot over a map 
2) plot_track: show track over a map 
3) add_text: show text over a map 
4)  plot_zonal_mean: plot zonal mean data


"""


from pylab import *
from numpy import *
# from matplotlib.toolkits.basemap import Basemap, shiftgrid

try: 
    from matplotlib.toolkits.basemap import Basemap, shiftgrid
except ImportError:
    from mpl_toolkits.basemap import Basemap, shiftgrid


def_map_keywords=['map_proj',\
                      'show_lat',\
                      'show_lon',\
                      'lat_0',\
                      'lon_0', \
                      'minlat',\
                      'maxlat',\
                      'lat_inc',\
                      'minlon',\
                      'maxlon',\
                      'lon_inc',\
                      'lons',\
                      'lats',\
                      'boundinglat',\
                      'drawcountry'
                      'linewidth',\
                      'linecolor'\
                      ]
    
def make_basemap(to_do, **keywords):
    
    """
    
    Inputs:
    ----------------------------------------------------
    1. to_do:<int>: the type of action (reserved for future use)
    
    
    2.keywords: <dict>: extra inputs
    ---reserved keywords:
    --->minv:<float>: minimal values
    --->maxv:<float>: maximal values
    --->dv:<float>: intervals 
    --->show_map:<int>: whether to draw ap or not 
    --->map_proj:<str>: map projection
    --->show_lat:<int>: draw latitude grid
    --->show_lon:<int>: draw longitude grid 
    --->lat_0=0.0:<float>: central latitude
    --->lon_0=0.0:<float> central longitude
    --->minlat:<float>: lower boundary of latitude
    --->maxlat:<float>: upper boundary of latitude
    --->lat_inc:<float>: latitude divide
    --->minlon:<float>: lower boundary of longitude
    --->maxlon:<float>: upper boundary of longitude
    --->lon_inc:<float>: longitude divide
    --->lons:<array>: longitude ticks
    --->lats:<array>: latitude ticks
    --->boundinglat:<float>: bounding latitude 
    --->drawcountry:<int>: if draw_country==1, country boundary will be drawn 
    --->linewidth:<float>: width of coast line and country line 
    --->linecolor:<str>:  linecolor
    

    returns:
    ====================================================
    1. m:<basemap>: class base map 
    
    """

    if (to_do==0):
        print 'do nothing'
        return None
    
    
    
    map_proj='cyl'
    
    show_lat=1
    show_lon=1
    lat_0=0.0
    lon_0=0.0
    
    minlat=-90.0
    maxlat=90.0
    lat_inc=30
    
    minlon=-180.0
    maxlon=180.0
    lon_inc=60.0
    
    boundinglat=45.0
    drawcountry=0

    linewidth=1
    linecolor='k'

    drawcountry=0
    lons=None
    lats=None
    
        
    for keyname in keywords:
        if (keyname=='show_lat'):
            show_lat=keywords['show_lat']
            
        elif (keyname=='show_lon'):
            show_lon=keywords['show_lon']
            
        elif (keyname=='map_proj'):
            map_proj=keywords['map_proj']
            
        elif (keyname=='lat_0'):
            lat_0=keywords['lat_0']
            
        elif (keyname=='minlat'):
            minlat=keywords['minlat']
            
        elif (keyname=='maxlat' in keywords):
            maxlat=keywords['maxlat']
            
        elif (keyname=='lat_inc'):
            lat_inc=keywords['lat_inc']
            
        elif (keyname=='minlon'):
            minlon=keywords['minlon']
            
        elif (keyname=='lon_inc'):
            lon_inc=keywords['lon_inc']
            
        elif (keyname=='maxlon'):
            maxlon=keywords['maxlon']
        
        elif (keyname=='lon_0'):
            lon_0=keywords['lon_0']
            
        elif (keyname=='boundinglat'):
            boundinglat==keywords['boundinglat']
            
        elif (keyname=='linewidth'):
            linewidth==keywords['linewidth']

        elif (keyname=='linecolor'):
            linecolor=keywords['linecolor']
            
        elif (keyname=='drawcountry'):
            drawcountry=keywords['drawcountry']
        
        elif (keyname=='lons'):
            lons=keywords['lons']
 
        elif (keyname=='lats'):
            lats=keywords['lats']
        
        else:
            print 'make_map: unexpected keywords:', keyname
            
            
    
    print 'make_basema--map_proj:', map_proj
        
    if (map_proj=='npstere' or map_proj=='spstere'):
        m=Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
    elif (map_proj=='ortho'):
        m=Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
    else:
        
        
        m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                      urcrnrlon=maxlon, \
                      urcrnrlat=maxlat,\
                      projection=map_proj, \
                      lon_0=lon_0, \
                      lat_0=lat_0, \
                      resolution='l')
        
        
    if (drawcountry==1):
        m.drawcountries(color=linecolor, linewidth=linewidth)
    else:
        m.drawcoastlines(color=linecolor, linewidth=linewidth)
        
        m.drawmapboundary()
        
    if (show_lat==1):
        if (lats==None):
            lats=arange(minlat, maxlat+lat_inc, lat_inc)
        
        m.drawparallels(lats, labels=[1,0, 0,0])
        
    if (show_lon==1):
        if (lons==None):
            lons=arange(minlon,maxlon+lon_inc,lon_inc)
        
        m.drawmeridians(lons,labels=[0,0,0,1])
        
        
    return m
    
    

def plot_map(data, rlon=None, rlat=None, **keywords):
    
    """ display  the data over a map
    Inputs:
    ====================================================
    1. data:<array, (nlon, nlat)>: data over a map
    2. rlon:<array, (nlon)>: longitude
    3. rlat:<array, (nlat)>: latitude
    4. keywords: <dict>: extra inputs
    ---reserved keywords:
    -
    --->map_proj:<str>: map projection
    --->show_lat:<int>: draw latitude grid
    --->show_lon:<int>: draw longitude grid 
    --->lat_0=0.0:<float>: central latitude
    --->lon_0=0.0:<float> central longitude
    --->minlat:<float>: lower boundary of latitude
    --->maxlat:<float>: upper boundary of latitude
    --->lat_inc:<float>: latitude divide
    --->minlon:<float>: lower boundary of longitude
    --->maxlon:<float>: upper boundary of longitude
    --->lon_inc:<float>: longitude divide
    --->boundinglat:<float>: bounding latitude 
    --->drawcountry:<int>: if draw_country==1, country boundary will be drawn 
    --->linewidth:<float>: width of coast line and country line 
    --->linecolor:<str>:  linecolor
    --->cnlevels:<list>: levels for contour
    --->title:<str>: plot title
    --->add_str:<str>: add extra string
    
    --->minv:<float>: minimal values
    --->maxv:<float>: maximal values
    --->dv:<float>: intervals 
    --->use_ma:<int>: If use_ma=1, masked array will be used

    --->show_map:<int>: whether to draw ap or not 
    --->unit:<str>: unit to be shown in the little 
    --->cb:<int>: if cb==1, colorbar will be show
    --->cb_vert:<int>: if cb_vert==0, colorbar will be show at horizontal orientation 
    ---cmap:<cm>: class cm (colormap) 
    ---m: <basemap>: class base map. If m is not provided, 
       a new basemap object will be created
       
    
    Returns:
    ========================================================
    1. m: <Basemap>: class of Basemap 
    
    """
    
    def_keywords=['minv',\
                      'maxv',\
                      'dv',\
                      'use_ma',\
                      'show_map', \
                      'cnlevels',\
                      'add_str',\
                      'title',\
                      'unit',\
                      'cb',\
                      'cb_vert',\
                      'cb_title',\
                      'm',\
                      'cmap']
    
    
    
    m=None
    nlon, nlat=shape(data)
    use_ma=1
    if ('use_ma' in keywords):
        use_ma=keywords['use_ma']
        del keywords['use_ma']
    
    if (use_ma==0):
        
        
        vals=npy.array(data)
    
    else:
        
        vals=data.copy()
        
    minv=0.0
    maxv=0.0
    dv=0.0

    rlvl=None
    stitle=""
    cb_title=""
    
    cbar_vert=1
    cmap=cm.jet
    use_pcolor=0
    cb=1
    
    show_colorbar=1
    
    def_keywords=def_keywords+def_map_keywords
    use_map_keywords={}
    
    for keyname in keywords:
        if (keyname in def_map_keywords):
            use_map_keywords.update({keyname:keywords[keyname]})
            
        else:
            
            if (keyname=='minv'):
                minv=keywords['minv']
                
            elif (keyname=='maxv'):
                maxv=keywords['maxv']
                
            elif (keyname=='dv'):
                dv=keywords['dv']
                
            elif (keyname=='cnlevels'):
                rlvl=keywords['cnlevels']
                
            elif (keyname=='title'):
                stitle=keywords['title']
                
            elif (keyname=='unit'):
                add_str=keywords['unit']
                add_str=add_str.strip()
                stitle=stitle+'('+ add_str +')'
                
            elif (keyname=='use_pcolor'):
                use_pcolor=keywords['use_pcolor']
                
            elif (keyname=='cb'):
                cb=keywords['cb']
                
            elif (keyname=='cb_title'):
                cb_title=keywords['cb_title']
                
            elif (keyname=='cbar_vert'):
                cbar_vert=keywords['cbar_vert']
                
            elif (keyname=='cmap'):
                cmap=keywords['cmap']
                
            elif (keyname=='m'):
                m=keywords['m']
                
                     
            else:
                print 'Warning: plot_map'
                print 'unexpected keywords:', keyname
    

    if (cbar_vert==1):
        orientation='vertical'
    else:
        orientation='horizontal'
    
    if (maxv>minv):
        if (dv==0.0):
            dv=(maxv-minv)/10.0
        
        if (rlvl==None):
            rlvl=arange(minv, maxv+dv, dv)
    
    
    if (rlat==None):
        
        dlat=(maxlat-minlat)/nlat
        rlat=arange(minlat,maxlat, dlat)
    
    if (rlon==None):
        
        dlon=(maxlon-minlon)/nlon
        rlon=arange(minlon, maxlon, dlon)
        
        
    if (m==None):
        m=make_basemap(1, **use_map_keywords)
    
    maxlon=180.0
    minlon=-180.0
    
    if ('maxlon' in use_map_keywords):
        maxlon=use_map_keywords['maxlon']
    

    if ('minlon' in use_map_keywords):
        minlon=use_map_keywords['minlon']
    
    if (maxlon>180.0):
        shift_idx=where(rlon<0.0)

    else:
        shift_idx=where(rlon>180.0)
        
        
    shift_idx=squeeze(shift_idx)
    shift_idx=squeeze(shift_idx)
    nshift=size(shift_idx)
    
    if (nshift>0):
        # exchange longitude
        
        new_vals=zeros(shape(vals), float)
        new_rlon=zeros(nlon, float)
        new_rlon[0:nlon-nshift]=rlon[nshift:nlon]
        new_rlon[nlon-nshift:nlon]=rlon[0:nshift]
        
        if (maxlon>180.0):
            new_rlon[nlon-nshift:nlon]=new_rlon[nlon-nshift:nlon]+360.0

        else:
            new_rlon[0:nlon-nshift]=new_rlon[0:nlon-nshift]-360.0
            
        rlon=new_rlon

        new_vals[0:nlon-nshift, :]=vals[nshift:nlon, :]
        new_vals[nlon-nshift:nlon, :]=vals[0:nshift, :]
        vals=new_vals
    
    if ((maxlon-minlon)>=359.0):
        if (rlon[-1]<rlon[0]+360.0):
            rlon=resize(rlon, nlon+1)
            rlon[-1]=rlon[0]+360.0
            vals=squeeze(vals)
            vals=resize(vals, [nlon+1, nlat])
    
    x,y=m(*meshgrid(rlon, rlat))
    
    if (maxv>minv):
        if (use_pcolor==1):
            cs0=m.pcolormesh(x, y, transpose(vals),  shading='flat', \
                                 vmin=minv, vmax=maxv, cmap=cmap)
            
        else:
            if (rlvl<>None):
                cs0=m.contourf(x, y, transpose(vals), rlvl, cmap=cmap)
            else:
                cs0=m.contourf(x, y, transpose(vals), cmap=cmap)
    
                               
    else:
        if (use_pcolor==1):
            cs0=m.pcolormesh(x, y, transpose(vals),shading='flat', cmap=cmap)
        else:
            cs0=m.contourf(x, y, transpose(vals))
    
    title(stitle)
    
    # info(m.drawcoastlines)
    
    if (cb==1):
        colorbar(orientation=orientation)
        if (len(cb_title)>0): 
            title(cb_title)
    
    return m

def plot_track(rlon, rlat=None,**keywords):

    """ display track over a map
    Inputs:
    ====================================================
    1. rlon:<array, (nlon)>: track longitudes
    3. rlat:<array, (nlat)>: track latitude
    3. keywords: <dict>: extra inputs
    ---reserved keywords:

    --->map_proj:<str>: map projection
    --->show_lat:<int>: draw latitude grid
    --->show_lon:<int>: draw longitude grid 
    --->lat_0=0.0:<float>: central latitude
    --->lon_0=0.0:<float> central longitude
    --->minlat:<float>: lower boundary of latitude
    --->maxlat:<float>: upper boundary of latitude
    --->lat_inc:<float>: latitude divide
    --->minlon:<float>: lower boundary of longitude
    --->maxlon:<float>: upper boundary of longitude
    --->lon_inc:<float>: longitude divide
    --->boundinglat:<float>: bounding latitude 
    --->drawcountry:<int>: if draw_country==1, country boundary will be drawn 
    --->linewidth:<float>: width of coast line and country line 
    --->linecolor:<str>:  linecolor
    --->cnlevels:<list>: levels for contour
    --->title:<str>: plot title
    --->add_str:<str>: add extra string

    --->m:<basemap>: map to be shown. 
       If m==None, a new basemap object will be created
    --->title:<str>: plot title 
    --->unit:<str>: unit of title 
    
    --->cb:<int>: if cb==1, colorbar will be show
    --->cb_vert:<int>: if cb_vert==0, colorbar will be show at horizontal orientation 
    --->cb_title:<str>: title of the colobar
    --->cmap:<cm>: colormap 
    
    --->sgn:<str>: sign used to show the map 
    --->msize:<str>: size of the marker 
    --->color:<list, t:str>: color of marker 
    
    Returns:
    ========================================================
    1. m: <Basemap>: class of Basemap 
    
    """
    
    def_keywords=['m',\
                      'title',\
                      'unit',\
                      'cb',\
                      'cb_vert',\
                      'cb_title',\
                      'sgn',\
                      'msize',\
                      'color']
    

    

    use_map_keywords={}
    
    stitle=""
    cb_title=""
    cb=1
    cbar_vert=1
    cmap=cm.jet
    sgn='o'
    msize=8
    lcolor='k'
    m=None
    
    
    for keyname in keywords:
        if (keyname in def_map_keywords):
            use_map_keywords.update({keyname:keywords[keyname]})
            
        elif(keyname=='title'):
            stitle=keywords['title']
            
        elif (keyname=='unit'):
            add_str=keywords['unit']
            add_str=add_str.strip()
            stitle=stitle+'('+ add_str +')'

            
        elif (keyname=='cb'):
            cb=keywords['cb']
            
        elif (keyname=='cb_title'):
            cb_title=keywords['cb_title']
            
        elif (keyname=='cbar_vert'):
            cbar_vert=keywords['cbar_vert']
            
        elif (keyname=='cmap'):
            cmap=keywords['cmap']
            
        elif (keyname=='m'):
            m=keywords['m']
            
            
        elif (keyname=='sgn'):
            sgn=keywords['sgn']
            
            
        elif (keyname=='msize'):
            msize=keywords['msize']
            
            
        elif (keyname=='color'):
            lcolor=keywords['color']
            
                
        else:
            print 'Warning: plot_track'
            print 'unexpected keywords:', keyname
            

    
    if (m==None):
        
        m=make_basemap(1, **use_map_keywords)
    
    
    m.plot(rlon, rlat, sgn,color=lcolor, markersize=msize)
    
    
    if (cb==1):
        colorbar(orientation=orientation)
        if (len(cb_title)>0): 
            title(cb_title)
            
            
    
    return m

def add_text(rlon, rlat,txt, **keywords):
    
    """
    inputs:
    
    ====================================================
    1. rlon:<array, (nlon)>: track longitudes
    3. rlat:<array, (nlat)>: track latitude
    3. txt:<list, t:str>: text to be show 
    

    4. keywords: <dict>: extra inputs
    ---reserved keywords:
    
    --->map_proj:<str>: map projection
    --->show_lat:<int>: draw latitude grid
    --->show_lon:<int>: draw longitude grid 
    --->lat_0=0.0:<float>: central latitude
    --->lon_0=0.0:<float> central longitude
    --->minlat:<float>: lower boundary of latitude
    --->maxlat:<float>: upper boundary of latitude
    --->lat_inc:<float>: latitude divide
    --->minlon:<float>: lower boundary of longitude
    --->maxlon:<float>: upper boundary of longitude
    --->lon_inc:<float>: longitude divide
    --->boundinglat:<float>: bounding latitude 
    --->drawcountry:<int>: if draw_country==1, country boundary will be drawn 
    --->linewidth:<float>: width of coast line and country line 
    --->linecolor:<str>:  linecolor
    --->cnlevels:<list>: levels for contour
    --->title:<str>: plot title
    --->add_str:<str>: add extra string

    --->m:<basemap>: map to be shown. 
       If m==None, a new basemap object will be created
       
   
    --->title:<str>: plot title 
    --->unit:<str>: unit of title 
    
    --->color:<list, t:str>: color of marker 
    
    Returns:
    ========================================================
    1. m: <Basemap>: class of Basemap 
    
    """
    
    def_keywords=['m',\
                  'title',\
                  'unit']


    

    use_map_keywords={}
    
    stitle=""
    m=None
    
    
    for keyname in keywords:
        if (keyname in def_map_keywords):
            use_map_keywords.update({keyname:keywords[keyname]})

            
        elif (keyname=='title'):
            stitle=keywords['title']
            
        elif (keyname=='unit'):
            add_str=keywords['unit']
            add_str=add_str.strip()
            stitle=stitle+'('+ add_str +')'
                

        else:
            
            print 'Warning: add_text'
            print 'unexpected keywords:', keyname
            
    

    
    if (m==None):
        
        m=make_basemap(1, **use_map_keywords)
    
            
    x,y=m(rlon, rlat)
    
    for i in range(size(rlon)):
        text(x[i], y[i], txt[i], fontsize=10)
    
        
    if (len(stitle))>0:
        title(stitle)
    
    return m




def plot_zonal_mean(vals, pres_lvl, rlat,use_pcolor=0,  **keywords):
    
    """ plot zone mean 

    Inputs:
    -----------------------------------------
    1. vals:<array, (nlat, nz)>: zonal mean model or observation data
    2. pres_lvl:<array, nz>: pressure level
    3. rlat:<array, (nlat)>: latitude
    4. use_pcolor:<int>: if use_pcolor=1, pcolor instead of contourf will be used 
    5. keywords:<dict>: options
    ---reserved keywords:
    -->use_log:<int>: use_log=1, pressure will be shown  in log10 
    --->minv:<float>: minimal values
    --->maxv:<float>: maximal values
    --->dv:<float>: intervals 
    --->title:<str>: plot title 
    --->unit:<str>: unit of title 
    
    --->cb:<int>: if cb==1, colorbar will be show
    --->cb_vert:<int>: if cb_vert==0, colorbar will be show at horizontal orientation 
    --->cb_title:<str>: title of the colobar
    --->cmap:<cm>: colormap 
    
    --->bottom_pres:<float>: pressure at bottom level
    --->top_pres:<float>: pressure at top level 
    

    """
    
    
    use_log=0
    minv=0.0
    maxv=0.0
    dv=0.0
    stitle=""
    
    cb_title=""
    cb=1
    cbar_vert=1
    cmap=cm.jet
    
    bottom_pres=max(pres_lvl)
    top_pres=min(pres_lvl)
    
    
    for keyname in keywords:
        
        if (keyname=='use_log'):
            use_log=keywords['use_log']
        
        elif (keyname=='minv'):
            minv=keywords['minv']
        
        elif (keyname=='maxv'):
            maxv=keywords['maxv']
            
        elif (keyname=='dv'):
            
            dv=keywords['dv']
            
        elif (keyname=='title'):
            stitle=keywords['title']
            
        elif (keyname=='unit'):
            
            add_str=keywords['unit']
            add_str=add_str.strip()
            stitle=stitle+'('+ add_str +')'

        elif (keyname=='bottom_pres'):
            bottom_pres=keywords['bottom_pres']

        elif (keyname=='top_pres'):
            top_pres=keywords['top_pres']

            
        elif (keyname=='cb'):
            cb=keywords['cb']
            
        elif (keyname=='cb_title'):
            cb_title=keywords['cb_title']
            
        elif (keyname=='cbar_vert'):
            cbar_vert=keywords['cbar_vert']
            
        elif (keyname=='cmap'):
            cmap=keywords['cmap']

        else:

            
            print 'Warning: plot_zonal_mean'
            print 'unexpected keywords:', keyname

            
    if (use_log==1):

        cur_ax=gca()
        cur_ax.set_yscale('log')
    

    if (dv==0.0):
      	dv=(maxv-minv)/10.0

    if (maxv>minv):
        rlvl=arange(minv, maxv+dv, dv)
#        rlvl[0]=-999.0
#        rlvl[size(rlvl)-1]=999.
    
    vals=squeeze(vals)
    
    nx, ny=shape(vals)
    
    y=array(pres_lvl)
    if (rlat<>None):
        x=array(rlat)
    else:
        x=arange(-90.0, 90.1, 180./(nx-1.0))
    

    if (maxv>minv):
        if (use_pcolor==1):
            cs0=pcolor(x, y, transpose(vals),  shading='flat', \
                           vmin=minv, vmax=maxv, cmap=cmap)
            
        else:
            cs0=contourf(x, y, transpose(vals), rlvl, cmap=cmap)
    else:
        if (use_pcolor==1):
            cs0=pcolor(x, y, transpose(vals),shading='flat', cmap=cmap)
        else:
            cs0=contourf(x, y, transpose(vals))
    
    if (cbar_vert==1):
        orientation='vertical'
    else:
        orientation='horizontal'


    if (len(stitle)>0):
        title(stitle)
        
    if (cb==1):
        colorbar(orientation=orientation)
    
        if (len(cb_title)>0):
            title(cb_title)
    
    
    ylim([bottom_pres, top_pres])
    
    xmin=x[0]
    xmax=x[-1]
    xmin=int(xmin+0.5)
    xmax=int(xmax+0.5)
    
    xlim([xmin, xmax])
    
    
    xlabel('latitude')
    ylabel('Pressure (hPa)')
    
    return cs0

def show_plot():
    show()
    

if (__name__=="__main__"):
    
    rlon=arange(-180, 180, 10)
    rlat=arange(-90, 90, 10)
    lon_m, lat_m=meshgrid(rlon, rlat)
    data=sin(lon_m*pi/180.0)**2+cos(lat_m*pi/180.0)**2
    data=transpose(data)
    subplot(2,1,1)
    plot_map(data, rlon, rlat, use_pcolor=1, maxv=1.0, minv=-1.0, dv=0.2)
    show()
    
    
    
    
