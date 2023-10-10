from pylab import *
from numpy import *
import time_module as tm
import gp_axis as gax
import oco_feedback as ofb
import oco_units as ocunit
import geos_chem_def as gcdf 
import surface_co2_jacobian as co2s
import os

black_stn_list={'wkt':[31.32, -97.33],\
                'bao':[40.05,-105.00],\
                'wgc':[38.26, -121.49],\
                'nha':[42.95, -70.63],\
                'thd':[41.05, -124.15],\
                'esp':[49.58, -126.37],\
                'wsa':[49.93, -60.02],\
                'lef':[45.95, -90.27],\
                'dnd':[48.38, -99.00]\
                }


black_stn_list={'wkt':[31.32, -97.33],\
                'bao':[40.05,-105.00],\
                'wgc':[38.26, -121.49],\
                'nha':[42.95, -70.63],\
                'dnd':[48.38, -99.00],\
                'wbi':[41.72, -91.35],\
                'sgp':[36.80,-97.50], \
                'lef':[45.95, -90.27],\
                'amt':[45.03, -68.68],\
                'sct':[33.41,  -81.83],\
                'fsd':[49.88, -81.57],\
                'key':[25.67, -80.16]\
                }



class surface_gp:
    """ observation operator for satellite observations """
    
    def __init__(self, viewtype, obs_datapath,id, mean_update, hm_update):
        self.id=id
        self.viewtype=viewtype
        self.obs_datapath=obs_datapath
        self.stn_lat, self.stn_lon, self.stn_z=None, None, None
        self.stn_name=None
        self.grid_time=None
        self.mean_update=mean_update
        self.hm_update=hm_update
        
        self.otime, self.olon, self.olat, self.oflag=None, None, None, None
        self.oz, self.osf, self.obs, self.oerr=None, None, None, None
        self.obs_platform="surface"
        
        self.nobs=0
               
        
    
    def read_obs(self, yyyy, doy):
        yyyy, mm,dd=tm.doy_to_time_array(doy, yyyy)
        sdate=r'%4.4dD%3.3d' % (yyyy, doy)
        ncflnm=self.obs_datapath+"/"+self.viewtype+"."+sdate+".nc"
        print 'read obs from:', ncflnm
        
        if (not os.path.isfile(ncflnm)):
            self.nobs=0
            self.stn_lat=[]
            self.stn_lon=[]
            self.oz=[]
            self.obs=[]
            self.olat=[]
            self.olon=[]
            self.oerr=[]
            
            
            return
        
        varnames=['stn', 'grid_time', 'stn_lat', 'stn_lon', 'stn_z', 'altitude', \
                  'latitude', 'longitude', 'obs', 'sq_obs','err', 'time', 'flag']

        
        
        self.stn_name, self.grid_time, \
                       self.stn_lat, self.stn_lon, self.stn_z, \
                       self.oz, self.olat, \
                       self.olon, self.obs, \
                       self.sq_obs, self.oerr, \
                       self.otime, self.oflag=ofb.ncf_read(ncflnm, varnames)
        
        self.osf=ones(shape(self.stn_lat))
        self.nobs=size(self.obs)
        
        # nobs=size(olat)
        
        # used_idx=where(oerr>0.0)
        # used_idx=squeeze(used_idx)
        
        # self.olat=array(olat[used_idx])
        # self.olon=array(olon[used_idx])
        # self.oz=array(oz[used_idx])
        # self.otime=array(otime[used_idx])
        # self.obs=array(obs[used_idx])
        # self.osf=array(osf[used_idx])
        # self.oerr=sqrt(2)*1.0e6*array(oerr[used_idx])
        # self.odoy=doy
        # self.nobs=size(self.obs)
        
        
        
    def get_stn_obs_mean(self, lt_st, lt_end):
        """ if we plan to use time mean instead of hourly data
        """
        
        nstn=size(self.stn_lon)
        stn_obs=zeros(nstn, float)
        
        stn_obs[:]=-999.0
        stn_err=zeros(nstn, float)
        stn_err[:]=-999.0

        stn_otime=zeros(nstn, float)
        stn_otime[:]=-999.0
        stn_flag=zeros(nstn)
        stn_flag[:]=4

        
        
        for istn in range(nstn):
            sel_lon=self.stn_lon[istn]
            usd_idx=tm.get_ut_time_slot(lt_st, lt_end, self.grid_time, sel_lon)
            
            if (size(usd_idx)>0):
                print usd_idx
                
                sel_obs=self.obs[usd_idx, istn]
                sel_sq_obs=self.sq_obs[usd_idx, istn]
                sel_err=self.oerr[usd_idx, istn]
                sel_flag=self.oflag[usd_idx, istn]
                sel_otime=self.otime[usd_idx, istn]
                sel_oz=self.oz[usd_idx, istn]
                
                # sub_idx=where((sel_obs>0) & (sel_flag==0) & (sel_oz<=1.0))
                sub_idx=where((sel_obs>0) & (sel_flag==0) & (sel_oz<=3.5))
                
                sub_idx=squeeze(sub_idx)
                
                sel_nobs=size(sub_idx)
                # print 'sel_nobs', sel_nobs
                if (sel_nobs>0):
                    
                    # sub_idx=squeeze(sub_idx)
                    # print sub_idx
                    
                    # print sel_obs[sub_idx]
                    sel_obs=mean(sel_obs[sub_idx])
                    sel_otime=mean(sel_otime[sub_idx])
                    sel_err=sum(sel_err[sub_idx]**2)
                    sel_sq_obs=sum(sel_sq_obs[sub_idx]**2)
                    sel_sq_obs=sel_sq_obs-sel_nobs*sel_obs**2
                    
                    
                    print 'istn, sel_obs, sel_sq_obs, sel_err', istn, \
                          sel_obs, sel_sq_obs, sel_err
                    
                    
                    # chose the maximal one
                    
                    # sel_err=max(sel_err, sel_sq_obs)
                    
                    
                    # if (sel_sq_obs>0.0):
                    
                        
                    sel_err=sqrt(sel_err/sel_nobs)
                    stn_obs[istn]=sel_obs
                    stn_err[istn]=sel_err
                    stn_flag[istn]=0
                    stn_otime[istn]=sel_otime

        
        for bk_stn in black_stn_list.keys():
            blat, blon=black_stn_list[bk_stn]
            black_id=where(abs((self.stn_lat-blat)<=1) & (abs(self.stn_lon-blon)<=1))
            black_id=squeeze(black_id)
            stn_err[black_id]=-999.0
        
            
        return stn_obs, stn_err, stn_otime, stn_flag
    
                
                    
    
    def read_sel_xgp(self, yyyy, mm, dd, \
                     est, eend, \
                     sel_em,\
                     datapath,\
                     mod_offset=0.0e-6,\
                     do_update=True, \
                     do_save=True):
        
        """ read selected tracer xgp values from files in datapath """
        
        # import co2_jacobian as ddcl
        out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        
        daily_obsflnm=gcdf.inv_path+"/"+self.viewtype+"_obs_mean"+"."+out_sdate+".nc"
        
        if (do_update):
            mod_xgp=co2s.co2jacob_sel(yyyy, mm, dd, \
                                      self.stn_lon, self.stn_lat, self.stn_z,\
                                      [est], [eend],\
                                      sel_em,\
                                      mod_offset=mod_offset,\
                                      datapath=datapath,\
                                      do_debug_jco2=False)
            if (do_save):
                nstn=size(self.stn_lon)
                dimnames=['nstn']
                dimtypes=['i']
                xno=range(nstn)
                dimvars=[xno]
                
                mod_info=ofb.geos_varinfo('mod', 'f', ['nstn'], mod_xgp)
                
                
                
                lon_info=ofb.geos_varinfo('lon', 'f', ['nstn'], self.stn_lon)
                lat_info=ofb.geos_varinfo('lat', 'f', ['nstn'], self.stn_lat)
                z_info=ofb.geos_varinfo('z', 'f', ['nstn'], self.stn_z)

                ofb.ncf_write_by_varinfo(daily_obsflnm, dimnames, dimtypes, dimvars, [lon_info,\
                                                                                      lat_info,\
                                                                                      z_info,\
                                                                                      mod_info])
                
        else:
            
            varnames=['mod']
            mod_xgp=ofb.ncf_read(daily_obsflnm, varnames)
            mod_xgp=squeeze(mod_xgp)
            
        return mod_xgp
    
               
        
    
    def read_simulated_obs(self, yyyy, mm, dd, \
                           datapath):
        
        """ read selected tracer xgp values from files in datapath """
        
        # import co2_jacobian as ddcl

        out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        
        daily_obsflnm=datapath+"/"+self.viewtype+"_obs_mean"+"."+out_sdate+".nc"
        
        
        varnames=['lon', 'lat', 'z', 'mod']
        lon, lat, z, mod_xgp=ofb.ncf_read(daily_obsflnm, varnames)
        mod_xgp=squeeze(mod_xgp)
        return mod_xgp
    
    
    
    
    def read_jacobian(self, yyyy, mm, dd,
                      em_st_list, em_end_list,\
                      em_yst_list,\
                      em_step_list,\
                      datapath,\
                      mod_offset=0.0e-6,\
                      do_update=True, \
                      do_save=True):
        
        out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        daily_obsflnm=gcdf.inv_path+"/"+self.viewtype+"_obs_hm"+"."+out_sdate+".nc"


        out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        daily_obsflnm=gcdf.inv_path+"/"+self.viewtype+"_obs_hm"+"."+out_sdate+".nc"
        
        if (do_update):
            hm=co2s.co2jacob(yyyy, mm, dd, \
                             self.stn_lon, self.stn_lat, \
                             self.stn_z,\
                             em_st_list, em_end_list, \
                             em_yst_list,\
                             em_step_list,\
                             mod_offset=mod_offset,\
                             datapath=datapath,\
                             do_debug_jco2=False)
            
            
            

            if (do_save):
                sav_nobs, sav_ne=shape(hm)
                dimnames=['nstn', 'ne']
                dimtypes=['i', 'i']
                xnx, xne=arange(sav_nobs), arange(sav_ne)
                dimvars=[xnx, xne]
                
                hm_info=ofb.geos_varinfo('hm', 'f', ['nstn', 'ne'], hm)
                lon_info=ofb.geos_varinfo('lon', 'f', ['nstn'], self.stn_lon)
                lat_info=ofb.geos_varinfo('lat', 'f', ['nstn'], self.stn_lat)
                z_info=ofb.geos_varinfo('z', 'f', ['nstn'], self.stn_z)
                
                ofb.ncf_write_by_varinfo(daily_obsflnm, dimnames, \
                                         dimtypes, dimvars, [lon_info,\
                                                             lat_info,\
                                                             z_info,\
                                                             hm_info])
                
        else:
            varnames=['hm']
            print daily_obsflnm
            
            hm=ofb.ncf_read(daily_obsflnm, varnames)
            hm=squeeze(hm)


        return hm
    
        
    

    

            
                
                
