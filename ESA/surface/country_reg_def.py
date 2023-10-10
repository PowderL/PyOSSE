""" this module is used to define the sub regions based on
    country map and Transcom region definition

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.21
    
    
    Classes:
    ===============================================
    1. country_region: class for basis functions 
    
"""

import ESA.util.gen_plots as gpl

from pylab import *

from numpy import *
import numpy as npy

import Scientific.IO.NetCDF as nfio


# spatial resolution of the sub-regions for snow,  land or ocean t3 regions 

flnm_flux_map_t3='region_flux_ml.nc'

  

sp_reg_dict={}



flnm_sub_reg_flux_out='cell_flux_ml.nc'



#  definition of the class

class country_region:
    """
    class for country region
    
    
    """
    def __init__(self, lon, lat, xid, yid, reg_id):
        self.lon=lon
        self.lat=lat
        self.xid=xid
        self.yid=yid
        self.reg_id=reg_id
        self.nsel=size(xid)
        

class country_map:
    """
    class for country region
    
    
    """
    
    def __init__(self, sel_country, cn_map):
        
        """
        Initialization
        
        """
        nf=nfio.NetCDFFile(cn_map, 'r')
        dig_map=nf.variables['map']
        dig_map=npy.array(dig_map)
        cn_id_list=nf.variables['country_id']
        atr_line=getattr(cn_id_list, 'country_name')
        country_name_list=atr_line.split(';')
        ientry=-1
        
        for cn in country_name_list:
            cn=cn.strip()
            if (cn==sel_country):
                ientry=ientry+1
                break
    
            ientry=ientry+1
        
        cn_id_list=npy.array(cn_id_list)
        cn_id=cn_id_list[ientry]
        dig_map=dig_map.astype(integer)
        # print 'cn_id', ientry, cn_id
        
        
        sel_cells=npy.where(dig_map==cn_id)
        sel_cells=npy.squeeze(sel_cells)
        idx_lon=sel_cells[0]
        idx_lat=sel_cells[1]
        self.cn_xid=array(idx_lon)
        self.cn_yid=array(idx_lat)

        self.sel_country=sel_country
        self.country_id=cn_id
        self.cn_map=dig_map
        self.cn_map[sel_cells]=cn_id
        
        
        map_lon=nf.variables['longitude']
        map_lon=npy.array(map_lon)
        map_lat=nf.variables['latitude']
        map_lat=npy.array(map_lat)
        nf.close()
        
        self.map=dig_map
        self.lon=map_lon
        self.lat=map_lat
        
        self.country_name_list=country_name_list
        self.country_id_list=cn_id_list

        self.nsel_cl=0
        self.cell_xid=None
        self.cell_yid=None
        self.sel_reg=list()

        self.nsel_reg=0
        
        
        
    def sel_cell(self,  lon_l, lon_r, \
                     lat_l, lat_r,in_map=None, in_pid=None):
        
        sel_cell_xid=list()
        sel_cell_yid=list()
        idx_lon=self.cn_xid
        idx_lat=self.cn_yid
        ncn_cell=npy.size(idx_lon)
        map_lon=self.lon
        map_lat=self.lat

        
        if (in_map==None):
            
            sel_map=npy.zeros(shape(self.map))
        else:
            sel_map=in_map

        if (in_pid==None):
            sel_pid=list()
        else:
            sel_pid=in_pid

        
        cid=npy.max(sel_map.flat)+1
        
        nsel=0
        for isel in range(ncn_cell):
            ix=idx_lon[isel]
            iy=idx_lat[isel]
            # print ix, iy
            
            if ((map_lon[ix]>=lon_l) & (map_lon[ix]<=lon_r) &\
                (map_lat[iy]>=lat_l) & (map_lat[iy]<=lat_r)):
                
                sel_cell_xid.append(ix)
                sel_cell_yid.append(iy)
                sel_map[ix, iy]=cid
                sel_pid.append(self.country_id)
                cid=cid+1
        
                nsel=nsel+1
                
                
        
        self.nsel_cl=nsel
        self.cell_xid=sel_cell_xid
        self.cell_yid=sel_cell_yid
        # gpl.plot_map(sel_map, map_lon, map_lat, use_pcolor=1)
        # show()
        return sel_map, sel_pid
    
    

    def divid_country(self, xst, xend, yst, yend, dx, dy, in_map=None):
        
        if (in_map==None):
            cur_map=in_map
        else:
            cur_map=self.cn_map

            
        ix=npy.argmin(abs(cmap2.lon-xst))
        wnd_xst=cmap2.lon[ix]
        idx=npy.argmin(abs(cmap2.lon-xend))
        wnd_xend=cmap2.lon[ix]

        # x coordinate
        wnd_x=npy.arange(wnd_xst, wnd_xend+dx, dx)
        
    
        ix=npy.argmin(abs(cmap2.lat-yst))
        wnd_yst=cmap2.lat[ix]
        idx=npy.argmin(abs(cmap2.lat-yend))
        wnd_yend=npy.cmap2.lon[ix]


        # y coordinate
        
        wnd_y=npy.arange(wnd_yst, wnd_yend+dy, dy)
        
        
        idx=npy.where(cur_map>0)
        cell_x=idx[0]
        cell_y=idx[1]

        nx=npy.size(wnd_x)
        ny=npy.size(wnd_y)
        ncell=npy.size(cell_x)
        
        wnd_cnt=npy.zeros([nx, ny])
        
        
        x=self.lon[cell_x]
        y=self.lat[cell_y]
        
        div_map=zeros(shape(cur_map))
        ncnt=0

        for ic in range(ncell):
            if ((x[ic]>=wnd_xst) & (x[ic]<=wnd_xend) & (y[ic]>=wnd_yst) & (y[ic]<=wnd_yend)):                
                ix=npy.argmin(npy.abs(wnd_x-x[ic]))
                iy=npy.argmin(npy.abs(wnd_y-y[ic]))
            
                if (wnd_cnt[ix, iy]==0):
                    ncnt=nct+1
                    div_map[x[ic], y[ic]]=ncnt
                    wnd_cnt[ix,iy]=ncnt
                else:
                    div_map[x[ic], y[ic]]=wnd_cnt[ix,iy]
                    
        return div_map
    
            
    def sel_cell(self,  lon_l, lon_r, lat_l, lat_r,in_map=None, in_pid=None):
        
        sel_cell_xid=list()
        sel_cell_yid=list()
        idx_lon=self.cn_xid
        idx_lat=self.cn_yid
        ncn_cell=npy.size(idx_lon)
        map_lon=self.lon
        map_lat=self.lat

        
        if (in_map==None):
            
            sel_map=zeros(shape(self.map))
        else:
            sel_map=in_map

        if (in_pid==None):
            sel_pid=list()
        else:
            sel_pid=in_pid

        
        cid=max(sel_map.flat)+1
        
        nsel=0
        for isel in range(ncn_cell):
            ix=idx_lon[isel]
            iy=idx_lat[isel]
            # print ix, iy
            
            if ((map_lon[ix]>=lon_l) & (map_lon[ix]<=lon_r) &\
                (map_lat[iy]>=lat_l) & (map_lat[iy]<=lat_r)):
                
                sel_cell_xid.append(ix)
                sel_cell_yid.append(iy)
                sel_map[ix, iy]=cid
                sel_pid.append(self.country_id)
                cid=cid+1
        
                nsel=nsel+1
                
                
        
        self.nsel_cl=nsel
        self.cell_xid=sel_cell_xid
        self.cell_yid=sel_cell_yid
        # gpl.plot_map(sel_map, map_lon, map_lat, use_pcolor=1)
        # show()
        return sel_map, sel_pid
    
    

    
    def div_region(self, lon_l_lst, lon_r_lst, lat_l_lst, lat_r_lst, in_map=None, in_pid=None):
        nreg=len(lon_l_lst)
        
        idx_lon=self.cn_xid
        idx_lat=self.cn_yid
        ncn_cell=size(idx_lon)
        map_lon=self.lon
        map_lat=self.lat
        reg_id=1
        sel_reg_lst=list()
        if (in_map==None):
            sel_map=zeros(shape(self.map))
            
        else:
            sel_map=array(in_map)
        
        reg_id=max(sel_map.flat)+1

        if (in_pid==None):
            sel_pid=list()
        
        else:
            sel_pid=in_pid
        
                          
        for ireg in range(nreg):

            lon_l=lon_l_lst[ireg]
            lon_r=lon_r_lst[ireg]
            lat_l=lat_l_lst[ireg]
            lat_r=lat_r_lst[ireg]
            sel_cell_xid=list()
            sel_cell_yid=list()
            
            nsel=0
                
            for isel in range(ncn_cell):
                ix=idx_lon[isel]
                iy=idx_lat[isel]
                # print ix, iy
                
                if ((map_lon[ix]>=lon_l) & (map_lon[ix]<=lon_r) &\
                    (map_lat[iy]>=lat_l) & (map_lat[iy]<=lat_r) & (sel_map[ix,iy]==0)):
                    nsel=nsel+1
                    sel_cell_xid.append(ix)
                    sel_cell_yid.append(iy)
                    sel_map[ix, iy]=reg_id
                    
            
            print 'nsel', nsel, lon_l, lat_l, lon_r, lat_r
            
            if (nsel>0):
                clreg=country_region(map_lon, map_lat, sel_cell_xid, sel_cell_yid, reg_id)
                reg_id=reg_id+1
                sel_reg_lst.append(clreg)
                sel_pid.append(self.country_id)
            
        self.sel_reg=sel_reg_lst
        self.nsel_reg=reg_id-1
        return sel_map, sel_pid
    
    


    
def divide_uk_ireland(do_debug=False):

    cmap1=country_map('United Kingdom', 'country_map_0.5x0.666.nc')
    cmap2=country_map('Ireland', 'country_map_0.5x0.666.nc')
    cmap3=country_map('France', 'country_map_0.5x0.666.nc')
    cmap4=country_map('Netherlands', 'country_map_0.5x0.666.nc')
    cmap5=country_map('Belgium', 'country_map_0.5x0.666.nc')
    cmap6=country_map('Germany', 'country_map_0.5x0.666.nc')
    
    
    lat_l_lst=[53,     53,   53,  56]
    lon_l_lst=[-10.0, -5.3, -3.4, -10]
    
    lat_r_lst=[56,    56,   56, 70]
    lon_r_lst=[-5.3,-3.4,   5,  5]

    sel_map, pid=cmap1.sel_cell(-15, 15, 30, 53)
    print max(sel_map.flat), '0'
    sel_map, pid=cmap1.div_region(lon_l_lst, lon_r_lst, \

                                 lat_l_lst, lat_r_lst, sel_map, pid)
    
    lat_l_lst=[48]
    lon_l_lst=[-40]
    
    lat_r_lst=[70]
    lon_r_lst=[40]
    
    print max(sel_map.flat), '1'
    sel_map, pid=cmap2.div_region(lon_l_lst, lon_r_lst, \
                                  lat_l_lst, lat_r_lst, sel_map, pid)



    print max(sel_map.flat), '2'
    

    sel_map, pid=cmap3.div_region(lon_l_lst, lon_r_lst, \
                                  lat_l_lst, lat_r_lst, sel_map, pid)



    print max(sel_map.flat), '3'

    sel_map, pid=cmap4.div_region(lon_l_lst, lon_r_lst, \
                                  lat_l_lst, lat_r_lst, sel_map, pid)
    


    print max(sel_map.flat), '4'


    sel_map, pid=cmap5.div_region(lon_l_lst, lon_r_lst, \
                                  lat_l_lst, lat_r_lst, sel_map, pid)
    


    print max(sel_map.flat), '5'

    sel_map, pid=cmap6.div_region(lon_l_lst, lon_r_lst, \
                                  lat_l_lst, lat_r_lst, sel_map, pid)
    


    print max(sel_map.flat), '6'
    
    
    if (do_debug==True):
        make_change=True
        if (make_change):
            sel_map=where(sel_map==68, 100, sel_map)
            sel_map=where(sel_map==69, 300, sel_map)
            sel_map=where(sel_map==70, 10, sel_map)
            sel_map=where(sel_map==71, 250, sel_map)
            sel_map=where(sel_map==72, 180, sel_map)
            sel_map=where(sel_map==73, 100, sel_map)
            sel_map=where(sel_map==74, 220, sel_map)
            # sel_map=where(sel_map==0, -99.0
            sel_map=where(sel_map>0,sel_map+10,sel_map)
            sel_map=where(sel_map<78,4*sel_map,sel_map)
            sel_map=where(sel_map==0, -1, sel_map)
        else:
            sel_map=where(sel_map>0,2*sel_map,sel_map)
            sel_map=where(sel_map==0, -1, sel_map)
        
        cx=cm.Paired
        # cx=cm.Set1
        cx=cm.jet
        cx.set_over('r')
        cx.set_under('w')
        
        m=gpl.plot_map(sel_map, cmap1.lon, cmap1.lat, maxlat=65, minlat=45, maxlon=5, minlon=-15, use_pcolor=1,cb=0, cmap=cx, minv=10, maxv=320,dv=2)
        # m=gpl.plot_map(sel_map, cmap1.lon, cmap1.lat, maxlat=65, minlat=45, maxlon=5, minlon=-15, use_pcolor=1,cb=0)

        # m.plot([-15, -15, 10, 10, -15], [45, 65, 70, 40, 40], linewidth='1.5', color='w')

        
        savefig('variables_map.png')
        show()
        # gpl.plot_map(sel_map, cmap.lon, cmap.lat, use_pcolor=1)
        # show()
        
    print len(pid)
    
    return sel_map, pid




        
if (__name__=='__main__'):
    
    
        
    divide_uk_ireland(do_debug=True)
    print 'OK'
    

    
