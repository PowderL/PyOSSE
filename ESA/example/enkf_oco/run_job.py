#! /geos/u23/epd7/bin/python2.7

import numpy as npy
import ESA.util.time_module as tm
import ESA.enkf.state_vector as stv_c
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio

import oco_assim_step as ost

yyyy, mm, dd=2009, 1,1
menu_flnm='osse_def.cfg'
temporal_resolution=32
ocs=ost.oco_assim_step_cl(yyyy, mm, dd,\
                     menu_flnm,\
                     temporal_resolution,\
                     err_scale=1.0,\
                     start_step=0, \
                     rerun_path='./'\
                     )
        


for i in range(0, 5):
    ocs.do_one_step(update_step=0, rerun_step=0,  hm_update_date='20081231')

    
    # #c: outputs
    
    xne=npy.arange(ocs.cl_stv.wnd_ne)
    xnx=npy.arange(ocs.cl_stv.wnd_nx)
    xwhole=npy.arange(npy.size(ocs.whole_flux))
    

    sstep=r'%2.2d' % (i)

    datapath=ocs.root_menu['inv.path']
    flnm=ocs.root_menu['inv.flnm']
    resflnm=datapath+flnm+'.'+sstep+'.nc'
    
    
    
    xne_info=ovar.io_var_cl('ne', 'f', ['ne'], xne, \
                                     varattr={'desc':'xne'})
    
    xnx_info=ovar.io_var_cl('nx', 'f', ['nx'], xnx, \
                                     varattr={'desc':'xnx'})
    
    xwhole_info=ovar.io_var_cl('nwhole', 'f', ['nwhole'], xwhole, \
                                   varattr={'desc':'xwhole'})
    
    
    
    x0_info=ovar.io_var_cl('x0', 'f', ['nx'], ocs.cl_stv.wnd_mean_x0)
    x_info=ovar.io_var_cl('x', 'f', ['nx'], ocs.cl_stv.wnd_mean_x)
    xinc_info=ovar.io_var_cl('xinc', 'f', ['nx'], ocs.cl_stv.wnd_xinc)
    
    xtm_info=ovar.io_var_cl('sum_xtm', 'f', ['ne', 'ne'], ocs.cl_stv.wnd_xtm)
    dx_info=ovar.io_var_cl('dx', 'f', ['nx', 'ne'], ocs.cl_stv.wnd_dx)
    
    
    # flux

    whole_flux0_info=ovar.io_var_cl('whole_flux0', 'f', ['nwhole'], ocs.whole_flux0)
    whole_flux_info=ovar.io_var_cl('whole_flux', 'f', ['nwhole'], ocs.whole_flux)
    whole_dx_flux_info=ovar.io_var_cl('whole_dx_flux', 'f', ['nwhole'], ocs.whole_dx_flux)
    
    
    ncfio.ncf_save_var(resflnm, [dx_info, xtm_info, x_info, \
                                     x0_info, xinc_info, \
                                     whole_flux_info, whole_flux0_info,\
                                     whole_dx_flux_info], \
                           [xne_info, xnx_info, xwhole_info],create_new=True) 
    
    
    

                                                                
