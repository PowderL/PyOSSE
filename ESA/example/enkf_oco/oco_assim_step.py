# print 'try to import matrix tool'
# import matrix_tool as mtool

# print 'after  matrix tool'

import numpy as npy

import ESA.util.time_module as tm
import ESA.enkf.state_vector as stv_c
import ESA.enkf.assim_def_m as defm
import ESA.enkf.assim_daily_obs as assim_daily 
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.util.gen_plots as gpl

import prior_co2_flux as co2prior



#  parameters to countrol the run
class oco_assim_step_cl:
    def __init__(self, \
                     yyyy, mm, dd,\
                     menu_flnm,\
                     temporal_resolution,\
                     err_scale=1.0,\
                     start_step=0, \
                     rerun_path='./'\
                     ):
        
        # S1: read in configuration value 
        self.root_menu=assim_daily.assim_read_config(menu_flnm)
        
        # S2: create objects used in data assimilation 
        
        # #c: the classes included in self.cl_def
        # 1. cl_obs: class for reading observation data
        # 2. cl_fc_prof: class for sampling single tracer model run
        # 3. cl_enr_prof: class for sampling 
        # 4. cl_sat_opr: class for observation operator (sampling models)
        # 5. cl_x2flux: class for projecting coefficients to fluxes
        # 6. cl_run_desc: class for description outputs from ensemble runs 
        # 7. cl_stv: state vector
    
        self.cl_def=assim_daily.assim_create_objects(yyyy, mm, dd, self.root_menu)
        # #T:  make a reference to these objects
        
        self.cl_obs=self.cl_def.cl_obs
        self.cl_fc_prof=self.cl_def.cl_fc_prof
        self.cl_enr_prof=self.cl_def.cl_enr_prof
        self.cl_sat_opr=self.cl_def.cl_sat_opr
        self.cl_x2flux=self.cl_def.cl_x2flux
        self.cl_run_desc=self.cl_def.cl_run_desc
        self.cl_stv=self.cl_def.cl_stv
        
        
        # #c: time step for ensemble run and forecast run 
        
        
        self.em_yyyy_st, self.em_doy_st, \
            self.em_yyyy_end, self.em_doy_end= self.cl_run_desc.extract_current_time_step(yyyy, mm, dd)
        
        # #c: forecast run step 
        
        self.fc_doy_st=self.em_doy_st[0]
        if (npy.size(self.em_doy_st)>1):
            self.fc_doy_end=self.em_doy_st[1]
        else:
            self.fc_doy_end=self.fc_doy_st
            
        
        self.fc_yyyy_st=self.em_yyyy_st[0]
        
        if (npy.size(self.em_yyyy_st)>1):
            self.fc_yyyy_end=self.em_yyyy_st[1]
        else:
            self.fc_yyyy_end=self.fc_doy_st

        
        
        # #T: create the CTM control for forecast 
        
        
        self.step=start_step
        
        self.cl_ctm=self.root_menu['ctm.fclass'](1,self.fc_doy_st, self.fc_yyyy_st, \
                                                     self.fc_doy_end, self.fc_yyyy_end,\
                                                     self.root_menu['ctm.runpath'], \
                                                     self.root_menu['ctm.datapath'], \
                                                     self.root_menu['ctm.tra_st'],\
                                                     self.root_menu['ctm.tra_end'],\
                                                     self.root_menu['ctm.finput_gen'],\
                                                     self.root_menu['ctm.ext'],\
                                                     self.root_menu['ctm.runtype'],\
                                                     self.root_menu['ctm.runscript'])
        
                                                   
        
        # #c: class for BASIS function

        self.cl_bf=self.root_menu['bf.fclass'](self.root_menu['bf.path'], \
                                                   self.root_menu['bf.flnm'], \
                                                   yyyy, mm, dd,  \
                                                   varname_lst=self.root_menu['bf.varname_lst'],\
                                                   varname_dict=self.root_menu['bf.varname_dict'],\
                                                   fopen=self.root_menu['bf.fopen'], \
                                                   fread=self.root_menu['bf.fread'], \
                                                   fget=self.root_menu['bf.fget'], \
                                                   fclose=self.root_menu['bf.fclose'])
        
        
        
        
        # S2: assim_create_objects(yyyy, mm, dd, root_menu
        
        self.first_run=True
        self.err_scale=err_scale
        self.do_retry=False
        
        # S3: projection to regional flux 
        
        self.whole_flux0=None
        self.whole_flux=None
        self.whole_dx_flux=None
        
        # #T: setup time steps for running GEOS-Chem with prior fluxes
        # #T: but the 
        
        
        self.cur_mm=mm
        self.cur_dd=dd
        self.cur_yyyy=yyyy
        
        days_of_year=365
        
        
        if (npy.mod(self.cur_yyyy,4)==0):
            days_of_year=366
            
        self.max_doy=days_of_year
            

        self.inv_time_step=temporal_resolution
        
        
        
        self.cur_doy=tm.day_of_year(self.cur_yyyy, self.cur_mm, self.cur_dd)
        # #c: cur step is the time step for whole time period
        self.cur_step=start_step
        # #c: cur_end_step is the time step within the lag window
        self.cur_wnd_step=0
        self.step_end_day=1000*self.fc_yyyy_end+self.fc_doy_end
        self.lag_window=self.root_menu['inv.lag_window']
        
        
        # #c: diag
        self.output_daily_obs=self.root_menu['diag.save_obs']
        
    
            
                
            
    
    def do_one_step(self, update_step=3, rerun_step=6, hm_update_date='20080101'):

        """ do one assimilation and re-run 
        the model to the next step  """

        self.cur_doy=self.cur_doy
        cur_doy=self.cur_doy
        # #c: if it is a new year, add the list to table
        
        yyyy, mm, dd=tm.doy_to_time_array(self.cur_doy, self.cur_yyyy)
        
        nsets=self.cl_run_desc.extract_coverage(yyyy, mm, dd)
        
        self.em_yyyy_st, self.em_doy_st, \
            self.em_yyyy_end, self.em_doy_end= self.cl_run_desc.extract_current_time_step(yyyy, mm, dd)
        
        
        istep=self.cur_step
        
        # #c: fill the time for forecast run 
        
        self.fc_doy_st=self.em_doy_st[istep]
            
        if (npy.size(self.em_doy_st)>istep+1):
            self.fc_doy_end=self.em_doy_st[istep+1]
        else:
            # then last one 
            self.fc_doy_end=self.em_doy_st[-1]
            
        
        self.fc_yyyy_st=self.em_yyyy_st[istep]
            
        if (npy.size(self.em_yyyy_st)>istep+1):
            self.fc_yyyy_end=self.em_yyyy_st[istep+1]
        else:
            self.fc_yyyy_end=self.em_yyyy_st[-1]
            
            
        self.step_end_day=1000*self.fc_yyyy_end+self.fc_doy_end
        
        
        if (self.cur_step>=rerun_step):
            
            # run the model 
            self.cl_ctm.rerun_geos_chem(self.cur_step,self.fc_doy_st, self.fc_yyyy_st, \
                                            self.fc_doy_end, self.fc_yyyy_end)
            
        else:
            # just set time 
            
            self.cl_ctm.set_time(self.cur_step,self.fc_doy_st, self.fc_yyyy_st, \
                                     self.fc_doy_end, self.fc_yyyy_end)
                      
        
        
        # # set state vector 
        
        
        self.set_stv(yyyy, mm, dd)
        
        # loop over the 
        print 'do_one_step: loop for assimilating daily observations:', self.cur_doy, self.step_end_day
        
        
        for idoy in range(self.cur_doy, self.cur_doy+self.inv_time_step):
            
            
            do_mean_update=True
            # #c: reset cur_yyyy, and max doy 
            if (idoy>self.max_doy):
                self.cur_yyyy=self.cur_yyyy+1
                idoy=1
                days_of_year=365
        
                if (npy.mod(self.cur_yyyy,4)==0):
                    days_of_year=366
                    
                self.max_doy=days_of_year
            
            
            self.cur_doy=idoy
            yyyy, mm,dd=tm.doy_to_time_array(self.cur_doy, self.cur_yyyy)
            full_yyyydoy=self.cur_yyyy*1000+idoy
            
            if (full_yyyydoy>=self.step_end_day):
                break
   
            
            # #c: loading coverage by ensemble runs 
            
            nsets=self.cl_run_desc.extract_coverage(yyyy, mm, dd)
            enr_yst_lst=self.cl_run_desc.st_year_lst
            enr_dst_lst=self.cl_run_desc.st_doy_lst
            enr_step_lst=self.cl_run_desc.step_lst
            enr_em_st_lst=self.cl_run_desc.st_no_lst
            enr_em_end_lst=self.cl_run_desc.end_no_lst

            
            print ' ---------- do_one_step ----------------' 
            print ' time:', yyyy, mm, dd
            
            self.cur_yyyy=yyyy
            self.cur_mm=mm
            self.cur_dd=dd
            
            
            full_yyyydoy=yyyy*1000+idoy
            
            out_sdate=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
                
            do_mean_update=True
            mean_save=True
            do_hm_update=False
            if (hm_update_date==None):
                do_hm_update=True
            elif (out_sdate>hm_update_date):
                do_hm_update=True
                
            
                    
            do_obs_write=False
                
            # S2: read in observations
    
            nobs=self.cl_sat_opr.read_obs(yyyy, mm, dd, viwemode=self.root_menu['sat_obs.viewmode'], \
                                              viewtype=self.root_menu['sat_obs.viewtype'])
            
            if (nobs>0):
                
            # S3: read mean upate 
    
                mod_xgp=self.cl_sat_opr.read_sel_xgp(yyyy, mm, dd, \
                                                         gpname_lst=self.root_menu['sat_opr.gpname_lst'],\
                                                         sel_idx=self.root_menu['sat_opr.sel_idx'],\
                                                         sel_traname=self.root_menu['sat_opr.traname'],\
                                                         do_update=do_mean_update, \
                                                         do_save=mean_save,\
                                                         xgp_flnm=self.root_menu['sat_opr.xgp_flnm'])
                
            # S4: read in hm 
    
                
                
                
                
                hm=self.cl_sat_opr.read_jacobian(yyyy, mm, dd, \
                                                     enr_yst_lst,\
                                                     enr_dst_lst,\
                                                     enr_step_lst,\
                                                     enr_em_st_lst,\
                                                     enr_em_end_lst,\
                                                     gpname=self.root_menu['sat_opr.traname'],\
                                                     gpname_lst=self.root_menu['sat_opr.gpname_lst'],\
                                                     append_data=False,\
                                                     do_update=do_hm_update, \
                                                     do_save=True,\
                                                     hm_flnm=self.root_menu['sat_opr.hm_flnm'])
                
                flux_pb=self.root_menu['inv.flux_pb']
                if (flux_pb>0.0):
                    add_mod_obs=npy.dot(hm, self.cl_stv.wnd_mean_x0)
                    mod_xgp=mod_xgp+add_mod_obs
                
    
                # S5: shift y and dy as a result of previous assimilations. 
    
        
                mean_y, dy=self.cl_stv.convolve_y_dy(mod_xgp, hm,  use_wnd_xtm=self.root_menu['stv.use_wnd_xtm'], \
                                                         use_wnd_inc_m=self.root_menu['stv.use_wnd_inc_m'])
                
                # S6: assimilation 
                
                yobs=self.cl_sat_opr.obs
                yerr=self.cl_sat_opr.oerr
                mod_err=self.root_menu['inv.mod_err']
                
                yerr=yerr*yerr+(1.0e-6*mod_err)*(1.0e-6*mod_err)
                
                xinc, xtm_inc_m, xtm, inc_m=self.cl_stv.do_assim(dy, mean_y, yobs, yerr)     
            
                
                if (self.output_daily_obs):
                    self.dump_daily_data(xinc, inc_m, mod_xgp, hm, mean_y)
                

            # if: nobs end
        # loop: idoy 
        
            
        self.do_assim_update()
        
            
        
        
        
        print 'prepare_restart'
        
        self.prepare_restart()
        # move to next step
        
        self.cur_step=self.cur_step+1
        
        return self.cur_step
    
    

    
    def dump_daily_data(self, xinc, inc_m, mod_xgp,hm,  mean_y):
        
        """
        save data to netcdf file 
        
        """
        
        sav_nobs, sav_ne=npy.shape(hm)
        dimnames=['nobs', 'ne','nx']
        dimtypes=['i', 'i', 'i']
        
        daily_obsflnm=self.root_menu['diag.obsflnm']
        invpath=self.root_menu['diag.path']
        print self.cur_yyyy, self.cur_mm, self.cur_dd
        
        sdate=r'%4.4d%2.2d%2.2d' % (self.cur_yyyy, self.cur_mm, self.cur_dd)
        daily_obsflnm=daily_obsflnm.replace('XDATEX', sdate)
        daily_obsflnm=invpath+daily_obsflnm
        
        xnobs, xne, xnx=npy.arange(sav_nobs), npy.arange(sav_ne), npy.arange(self.cl_stv.wnd_nx)
        
        nobs_info=ovar.io_var_cl('nobs', 'f', ['nobs'], xnobs, \
                                     varattr={'desc':'obs id'})
        
        ne_info=ovar.io_var_cl('ne', 'f', ['ne'], xne, \
                                   varattr={'desc':'ne'})
        
        
        nx_info=ovar.io_var_cl('nx', 'f', ['nx'], xnx, \
                                   varattr={'desc':'nx'})
                
        hm_info=ovar.io_var_cl('hm', 'f', ['nobs', 'ne'], hm, \
                                   varattr={'desc':'jacobian'})
        
        mean_x_info=ovar.io_var_cl('mean_x', 'f', ['nx'], self.cl_stv.wnd_mean_x, \
                                       varattr={'desc':'mean_x'})
        
        
        xinc_info=ovar.io_var_cl('xinc', 'f', ['nx'], xinc, \
                                     varattr={'desc':'analysis increment'})
        
        inc_m_info=ovar.io_var_cl('inc_m', 'f', ['ne'], inc_m, \
                                      varattr={'desc':'inncrement matrix'})
        
        
        mod_info=ovar.io_var_cl('mod', 'f', ['nobs'], mod_xgp, \
                                            varattr={'desc':'a prior model value', 'units':'ppm'})
        
        y_info=ovar.io_var_cl('correct_mod', 'f', ['nobs'], mean_y, \
                                  varattr={'desc':'corrected model value', 'units':'ppm'})
        
        
        obs=self.cl_sat_opr.obs
        oerr=self.cl_sat_opr.oerr

        obs_info=ovar.io_var_cl('obs', 'f', ['nobs'], obs, \
                                    varattr={'desc':'observation', 'units':'ppm'})
        
        oerr_info=ovar.io_var_cl('err', 'f', ['nobs'], oerr, \
                                     varattr={'desc':'observation error', 'units':'ppm'})

        olon=self.cl_sat_opr.olon
        
        
        olon_info=ovar.io_var_cl('lon', 'f', ['nobs'], olon,\
                                     varattr={"desc":"observation longitude",\
                                                  "units":"degrees_east",\
                                                  'standard_name':'longitude',\
                                                  'long_name':'Longitude, positive east'})
        
        
        olat=self.cl_sat_opr.olat
        olat_info=ovar.io_var_cl('lat', 'f', ['nobs'], olat,\
                                     varattr={"desc":"observation latitude",\
                                                  "units":"degrees_north",\
                                                  'standard_name':'latitude',\
                                                  'long_name':'Latitude, positive north'})
        

        otime=self.cl_sat_opr.otime
        otime_info=ovar.io_var_cl('time', 'f', ['nobs'], otime)
        
        ncfio.ncf_save_var(daily_obsflnm,  [hm_info, \
                                                obs_info, \
                                                olon_info,\
                                                olat_info,\
                                                oerr_info,\
                                                otime_info,\
                                                mod_info,\
                                                y_info,\
                                                mean_x_info, xinc_info, inc_m_info], \
                               [nobs_info, ne_info, nx_info],\
                               create_new=True)
        

        
            
            
                
        
    def set_stv(self,yyyy, mm, dd):
        
        """
        
        create and maintain stv 
        Inputs:
        ---------------------------------
        1. yyyy, mm, dd:<date>
        
        """
        import pylab as plb
        # S1: read BF (i.e, perturbation ensemble)
        
        bf_read_lst=self.root_menu['bf.rdname_lst']
        lon, lat, reg_map,  sf_area, reg_err=co2prior.read_reg_bf(self.cl_bf, yyyy, mm, dd, \
                                                                      varname_lst=bf_read_lst)
        
        # S2: 
        
        
        datapath=self.root_menu['prior.datapath']
        fctm=self.root_menu['prior.fctm']
        nmext=self.root_menu['prior.ext']
        ctm_date0=r'%4.4d%2.2d%2.2d' % (yyyy, mm, dd)
        ftracerinfo=self.root_menu['prior.ftracerinfo']
        fdiaginfo=self.root_menu['prior.fdiaginfo']
        
        tau0=tm.get_tau(yyyy, mm, dd)
        
        ems_tau0, ems_tau1, reg_src=co2prior.get_gc_prior_flux(yyyy, mm, dd,\
                                                                   datapath,\
                                                                   fctm,\
                                                                   nmext,\
                                                                   ctm_date0, \
                                                                   ftracerinfo,\
                                                                   fdiaginfo,\
                                                                   tau0, \
                                                                   reg_map,\
                                                                   sf_area)
        
        
        print npy.max(reg_src)
        print npy.min(reg_src)
        
        
        enlarge_factor=1.0
        ocean_err_factor=1.0 # 0.8 #   .0/2.0
        
        nland=100
        nx=npy.size(reg_err[:])
       
        # dx=identity(nx, float) # diag(reg_err[:])
        
        dx=npy.ones(nx, float)
        nocean=nx-nland
        
        
        # #T1: set x
        x=npy.zeros(nx, float)


        flux_pb=self.root_menu['inv.flux_pb']
        if (flux_pb>0.0):
            x=flux_pb*reg_src/reg_err
        
        # #T2: set dx

        dx[1:nland]=0.5*abs(reg_src[1:nland])/reg_err[1:nland]
        
        dx[nland:]=0.5*abs(reg_src[nland:])/reg_err[nland:]
        
        # #c: change dx to matrix 
        dx=npy.diag(dx)
        
        # #c: set first one to small amount
        dx[0, 0]=0.1
        
            
        cur_tag=r'ST%3.3f' % (self.cur_step)
        
        self.cl_stv.add_new_x_to_window(x, \
                                            dx, \
                                            ems_tau0,\
                                            ems_tau1,\
                                            step_tag=cur_tag,\
                                            enlarge_factor=1.0,\
                                            lag_window=self.lag_window,\
                                            cl_x2flux=None)
        
        
        
        
        # S4: keep history records
        
        if (self.whole_flux0==None):
            self.whole_flux0=reg_src+x*reg_err
        else:
            new_flux0=reg_src+x*reg_err
            self.whole_flux0=npy.concatenate((self.whole_flux0, new_flux0))
        
        if (self.whole_flux==None):
            self.whole_flux=reg_src+x*reg_err
        else:
            new_flux0=reg_src+x*reg_err
            self.whole_flux=npy.concatenate((self.whole_flux, new_flux0))
        
        if (self.whole_dx_flux==None):
            self.whole_dx_flux=reg_err
        else:
            self.whole_dx_flux=npy.concatenate((self.whole_dx_flux, reg_err))
            
        print 'set_stv: new size:', self.cl_stv.wnd_nx, self.cl_stv.wnd_ne
        return True
    
    
    def do_assim_update(self):
        """ upate the system using xtm and inc_m
        
        """
        # #! note that the size can be different, if bias is included
        
        xinc=self.cl_stv.wnd_mean_x-self.cl_stv.wnd_mean_x0
        xinc=xinc[self.cl_stv.nbias:]
        whole_nx=npy.size(self.whole_flux)
        wnd_nx=npy.size(xinc)
        pos0=whole_nx-wnd_nx
        flux_inc=xinc*self.whole_dx_flux[pos0:]
        self.whole_flux[pos0:]=self.whole_flux0[pos0:]+flux_inc
        

        return flux_inc
    
    
    
    def prepare_restart(self):
        """
        prepare restart files 

        """
        
        if (self.cl_stv.wnd_step>=self.cl_stv.lag_window):
        
            yyyy, mm, dd=tm.doy_to_time_array(self.cl_ctm.doy_st, self.cl_ctm.yyyy_st)
            nsets=self.cl_run_desc.extract_coverage(yyyy, mm, dd)
            
            enr_yst_lst=self.cl_run_desc.st_year_lst
            enr_dst_lst=self.cl_run_desc.st_doy_lst
            enr_step_lst=self.cl_run_desc.step_lst
            enr_em_st_lst=self.cl_run_desc.st_no_lst
            enr_em_end_lst=self.cl_run_desc.end_no_lst
            
            usd_lst=list()
            mdoy=npy.min(enr_dst_lst)
            myyyy=npy.min(enr_yst_lst)
        
            for iset in range(nsets):
                # #c: 'oldest one will be used and removed' 
            
                if ((enr_dst_lst[iset]==mdoy) and (enr_yst_lst[iset]==myyyy)):
                    usd_lst.append(iset)
        
        
            
            xinc=self.cl_stv.wnd_mean_x-self.cl_stv.wnd_mean_x0
            # #c: without bias term 
            
            # #c: xinc has a size of window,  larger than size of the time period. 
            # #c: however, usd_set_lst will be used to cut it off in transfer_restart_file
            
            xinc=xinc[self.cl_stv.nbias:]
                
            newdata=self.cl_ctm.transfer_restart_file(enr_em_st_lst, \
                                                          enr_em_end_lst, \
                                                          enr_yst_lst, \
                                                          enr_step_lst,\
                                                          usd_set_lst,\
                                                          xinc,\
                                                          self.root_menu['enr_prof.path'],\
                                                          self.root_menu['enr_prof.restart_flnm'])
        
            
        # add the adjust ment to restart file 
        
            self.cl_ctm.prepare_next_restart_file(newdata, resflnm=None, action=1)
        
        else:
            newdata=None
            
            self.cl_ctm.prepare_next_restart_file(newdata, resflnm=None, action=0)
                
                
            
        
        

if (__name__=='__main__'):
    
    print 'see run_job.py'
    
        
    
    







        
    
    
