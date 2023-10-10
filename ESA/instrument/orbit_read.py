from pylab import *
from numpy import *
import gen_plots as gpl
import time_module as  tm
import geos_chem_def as gcdf
import decode_line as decl
import numpy.random as rnd
import read_ecmwf_cld as rcld # read in cloud faction data
from Scientific.IO.NetCDF import *
import compute_model_grid as cgrd
import geo_constant as gc


try: 
    from matplotlib.toolkits.basemap import Basemap, shiftgrid
except ImportError:
    from mpl_toolkits.basemap import Basemap, shiftgrid


# default_pixelsize_oco= 0.25*pi*10.0 * 10.0 # km x km (diameter==10 km)

default_pixelsize_oco= 1.25 * 2.26 # km x km

grid_s0=gc.ra*gc.rb*(0.25*pi)**2

def penalty_function (pixel_size):
    # print pixel_size
    pf=(26.098*pixel_size**(-0.485)+10.18)/(26.098+10.18)
    # print pf
    return pf


class orbit_data:
    def __init__(self, yyyy, doy, viewmode, flnm=None, hh=13.5, \
                 ecmwf_cld_flnm='cld_ecmwf'):

        self.doy=doy
        self.yyyy=yyyy
        self.viewmode=viewmode
        
        # get_model_lat(model_res=None,nlat=0):
        
        # self.regrid_lon=arange(-180.0, 180.0, 1.0)
        # self.regrid_lat=arange(-90.0, 91.0, 1.0)
        self.regrid_lat=cgrd.get_model_lat(model_res='1x1')
        self.regrid_lon=cgrd.get_model_lon(model_res='1x1')
        
        self.nregrid_lon=size(self.regrid_lon)
        self.nregrid_lat=size(self.regrid_lat)
        self.regrid_clr_count=zeros([self.nregrid_lon, self.nregrid_lat], float)
        
        
        
        yyyy, mm, dd=tm.doy_to_time_array(doy, yyyy)
        syyyymm=r'%4.4d%2.2d' % (yyyy, mm)
        cld_flnm=ecmwf_cld_flnm+"."+syyyymm+".nc"
        cld_fl=NetCDFFile(cld_flnm, "r")
        tcc=None
        cld_lon, cld_lat, cur_tcc, tcc=rcld.get_cld_map(cld_fl, yyyy, mm,dd,hh, tcc)
        cld_fl.close()
        
        
        
        self.cnt_all_obs=0
        
        scenes_count=0
        obs_count=0
        
        if (flnm==None):
            self.data=None
        else:
            forb=open(flnm, 'r')
            line1=forb.readline()
            line2=forb.readline()
            terms=decl.decode_line(line2)
            rec_doy=int(terms[2])
            line3=forb.readline()
            lines=forb.readlines()
            forb.close()
            data=list()
            line_cnt=0
            for line in lines:
                line=line.strip()
                if (len(line)>0):
                    line=line.replace('-999', ' -999')
                    
                    terms=decl.decode_line(line)
                    terms=array(terms)
                    terms=terms.astype(float)
                    
                    nsounding=terms[2]
                    cflag=1
                    
                    if (nsounding>0):
                        self.cnt_all_obs= self.cnt_all_obs+ 4*nsounding
                        
                        obs_count=obs_count+4*nsounding
                        
                        scenes_count=scenes_count+1
                        if (terms[1]>180.0):
                            terms[1]=terms[1]-360.0
                            
                        orbit_doy=int(terms[3])
                        # check throught clound map
                        
                        # if (orbit_doy<rec_doy):
                        #    print line_cnt, terms[0], terms[1], terms[3], rec_doy
                        #    print line
                        
                        terms[3]=24.0*3600*(terms[3]-rec_doy+1)
                        olon, olat=terms[1], terms[0]
                        
                        dlon=abs(cld_lon-olon)
                        idx_lon=argmin(dlon)
                        dlat=abs(cld_lat-olat)
                        idx_lat=argmin(dlat)
                        cld_prob=cur_tcc[idx_lon, idx_lat]
                        if (viewmode=='nadir'):
                            pixel_size=default_pixelsize_oco
                            pf= penalty_function(default_pixelsize_oco)
                        else:
                            pixel_size=default_pixelsize_oco/abs(cos(terms[4]*pi/180.0)) #
                            pf= penalty_function(pixel_size)
                        
                        # pcld=(1.0-terms[-1])**(4*nsounding) # sqrt(nsounding)
                        # pcld=(1.0-terms[-1])**sqrt(nsounding)
                        # nclr=1.0*(1.0-cld_prob)*nsounding*pf
                        # nclr=4*nsounding
                        # print pf, cld_prob
                        
                        # zzzz=raw_input()
                        
                        cflag=1
                        rnd_cf=rnd.random(4*nsounding)
                        s_grid=cos(olat*pi/180.0)*grid_s0
                        s_grid=(1.0-cld_prob)*s_grid
                        s_pixel=4*nsounding*pixel_size
                        
                        
                        alignment_factor=(1.0-cld_prob)*(1.0-3.0*cld_prob*(1.0-cld_prob))
                        # the standard one
                        alignment_factor=0.5
                        if (s_grid>s_pixel):
                           
                        
                            clr_std=alignment_factor*(1.0-cld_prob)*pf
                        else:
                            alignment_factor= alignment_factor*s_pixel/s_grid
                            clr_std=alignment_factor*pf*(1.0-cld_prob)
                        
                        
                        clr_idx=where(rnd_cf<=clr_std)
                        
                        nclr=size(clr_idx)
                        if (nclr>0):
                            nold_clr=self.check_data_use(terms[1], terms[0], nclr)
                            if (nold_clr>=1):
                                cflag=1
                            else:
                                cflag=0
                                
                    
                    if (cflag==0):
                        # vals=[terms[3], terms[0], terms[1], terms[4], nsounding, terms[-1]]
                        
                        # print 'cflag, nclr', cflag, nclr
                        
                        vals=[terms[3], terms[0], terms[1], terms[4], nclr, terms[-1]]
                        vals=array(vals)
                        data.append(vals)
                        # print terms[0]
                        # testing
                        # if ((abs(terms[1]-60.0)<2.0) and (abs(terms[0]+45.7)<2.0)):
                        #    print 'local_time-blank', terms[3]/(24.0*3600)
                        #    txxxx=raw_input()
                            
                            
                        # if ((abs(terms[1]-60.0)<1.0) and (abs(terms[0]-45.7)<1.0)):
                        #    print 'local_time North', terms[3]/(24.0*3600)
                        #    txxxx=raw_input()
                line_cnt=line_cnt+1
            
            
            data=array(data)
            print 'clear data-shape', shape(data)
            print 'clear data-type', type(data)
            rtime=data[:,0]
            
            ord_idx=argsort(rtime)
            new_data=data[ord_idx, :]
            data=new_data
            self.data=array(data)
            print 'clear sum ()ll_clear)', sum(self.regrid_clr_count.flat)
            print 'direct clear sum 2 (i.e data shape)', sum(data[:, 4])
            
            self.recount_clr_data()
            # print 'driect sum 3', sum(data[:, 4])
            
            # print 'shape of data:', shape(data)
            # print 'scenes_count:', scenes_count
            print 'raw observation count:',  obs_count
            
    def check_data_use(self, lon, lat, nclr):
        lon_idx=searchsorted(self.regrid_lon, lon)
        lat_idx=searchsorted(self.regrid_lat, lat)

        if (lon_idx>=self.nregrid_lon):
            lon_idx=self.nregrid_lon-1

        if (lat_idx>=self.nregrid_lat):
            lat_idx=self.nregrid_lat-1
        
        old_nclr=self.regrid_clr_count[lon_idx, lat_idx]
        self.regrid_clr_count[lon_idx, lat_idx]=self.regrid_clr_count[lon_idx, lat_idx]+nclr
        return old_nclr
    
    def recount_clr_data(self):
        nobs=size(self.data[:,0])
        
        for iobs in range(nobs):
            lon_idx=searchsorted(self.regrid_lon, self.data[iobs, 2])
            lat_idx=searchsorted(self.regrid_lat, self.data[iobs, 1])
            
            if (lon_idx>=self.nregrid_lon):
                lon_idx=self.nregrid_lon-1
            if (lat_idx>=self.nregrid_lat):
                lat_idx=self.nregrid_lat-1
                
            
            self.data[iobs, 2]=self.regrid_lon[lon_idx]
            self.data[iobs, 1]=self.regrid_lat[lat_idx]
            
            
            old_nclr=self.regrid_clr_count[lon_idx, lat_idx]
            self.data[iobs, 4]=old_nclr
        
        
    def update_data(self, flnm=None, data=None):
        if (flnm==None):
            self.data=array(data)
        else:
            forb=open(flnm, 'r')
            line=forb.readline()
            lines=forb.readlines()
            data=list()
            for line in lines:
                terms=line.split()
                terms=array(terms)
                terms=terms.astype(float)
                data.append(terms)

    def show_orbit(self):
        data=self.data
        lon=data[:,2]
        lat=data[:,1]
        rtime=data[:,0]
        # rtime=rtime-rtime[0]
        ord_idx=argsort(rtime)
        lon=lon[ord_idx]
        lat=lat[ord_idx]
        rtime=rtime[ord_idx]
        
        maxlat=90.0
        minlat=-90.0
        maxlon=180.0
        minlon=-180.0
        lat_0=0.0
        lon_0=0.0
        ion()
        ax=subplot(2,1,1)
        m=Basemap(llcrnrlon=minlon, llcrnrlat=minlat, \
                  urcrnrlon=maxlon, urcrnrlat=maxlat,projection='cyl', lon_0=lon_0, lat_0=lat_0,resolution='l')
        m.drawcoastlines(linewidth=0.5)
        m.drawmapboundary()
        parallels = arange(-90.,90,30.)
        m.drawparallels(parallels,labels=[1,0,0,0])
        meridians = arange(-180.,180.,60.)
        m.drawmeridians(meridians)
        line1,=m.plot(lon, lat, 'b+')
        sec=int(rtime[4])
        tsec=r'Sec: %3.3d' % sec
        tt=text(0, -100,tsec) 
        nstep=size(lon)-5
        maxp=5
        for istep in range(nstep):
            maxp=maxp+1
            line1.set_data(lon[0:maxp], lat[0:maxp])
            sec=int(rtime[maxp-1])
            tsec=r'Sec: %3.3d' % sec
            tt.set_text(tsec)
            draw()

        
               
            
               
    
class oco_orbit:
    def __init__(self, yyyy_st, doy_st, yyyy_end,doy_end, \
                 viewtype="aqua",\
                 viewmode='nadir', \
                 datapath=gcdf.orbit_path):
        
        def_year=2006
        doy=doy_st
        yyyy=yyyy_st
        sdate=r'%4.4dD%3.3d' % (def_year, doy_st)
        
        self.data=list()
        self.doy=list()
        self.yyyy=list()
        doy=doy_st
        self.cnt_all_obs=0
        
        while True:
            sdate=r'%4.4dD%3.3d' % (def_year, doy)
            sdoy=r'%3.3d' % doy
            flnm=datapath+"/"+viewtype+"_num_"+viewmode+"_"+str(doy)+".dat"
            print flnm
            
            orbd=orbit_data(yyyy, doy, viewmode, flnm)
            self.cnt_all_obs=orbd.cnt_all_obs
            print sdoy,  self.cnt_all_obs
            
            self.data.append(orbd)
            self.doy.append(doy)
            self.yyyy.append(yyyy)
            yyyy, doy=tm.next_doy(yyyy, doy)
            # print yyyy, doy
            
            if (doy>doy_end and yyyy>=yyyy_end):
                break
        self.doy=array(self.doy)
        self.yyyy=array(self.yyyy)
        
            
            
    def print_orb_info(self):
        print '='*30+'Orbit information'+'='*30
        print  'number of date:', size(self.doy)
        
    def get_orbit(self, yyyy, doy):
        found=False
        for i in range(len(self.doy)):
            if (self.doy[i]==doy and self.yyyy[i]==yyyy):
                found=True
                break
        if (found):
            return self.data[i]
        else:
            return None
            
        
if (__name__=='__main__'):
    from pylab import *
    yyyy=2003
    doy=1
    oco=oco_orbit(yyyy, doy, yyyy, 2)
    sel_orbit=oco.get_orbit(2003, 1)
    # sel_orbit.show_orbit()
    counts=sel_orbit.data[:,4]
    print max(counts)
    
    # hist[counts]
    show()
    

            
            
                  
