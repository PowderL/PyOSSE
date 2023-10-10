# util 

pydoc -w ESA.util.bpch2_rw_smp
pydoc -w ESA.util.gc_ts_file_m
pydoc -w ESA.util.gen_plots
pydoc -w ESA.util.geo_constant
pydoc -w ESA.util.gp_axis_m
pydoc -w ESA.util.gp_data_m
pydoc -w ESA.util.gp_grid_m
pydoc -w ESA.util.horizontal_interp_m
pydoc -w ESA.util.__init__
pydoc -w ESA.util.line_process_m
pydoc -w ESA.util.message_m

pydoc -w ESA.util.otool_descfile_io
pydoc -w ESA.util.otool_gcfile_io
pydoc -w ESA.util.otool_grdfile_io
pydoc -w ESA.util.otool_menu_m
pydoc -w ESA.util.otool_ncfile_io
pydoc -w ESA.util.otool_obj
pydoc -w ESA.util.otool_txtfile_io
pydoc -w ESA.util.otool_var_io
pydoc -w ESA.util.pres_m
pydoc -w ESA.util.time_module
pydoc -w ESA.util.vertical_interp_m
# share library 
pydoc -w ESA.util.bpch2_rw_py
pydoc -w ESA.util.flib
pydoc -w ESA.util.flux_regrid
pydoc -w ESA.util.great_circle_distance
pydoc -w ESA.util.interpolate_f
pydoc -w ESA.util.process_nf_array
pydoc -w ESA.util.read_data
pydoc -w ESA.util.sample_model_field
pydoc -w ESA.util.sigma_pres_mod
pydoc -w ESA.util.vertical_column
pydoc -w ESA.util.vertical_profile


# atmosphere


pydoc -w ESA.atmosphere.compute_gc_grid
pydoc -w ESA.atmosphere.ctm_field_2d
pydoc -w ESA.atmosphere.ctm_field_m
pydoc -w ESA.atmosphere.ctm_grid_2d
pydoc -w ESA.atmosphere.ctm_grid_3d
pydoc -w ESA.atmosphere.ctm_profile_m
pydoc -w ESA.atmosphere.ctm_slice_m
pydoc -w ESA.atmosphere.ctm_world_m
pydoc -w ESA.atmosphere.enr_slice_m
pydoc -w ESA.atmosphere.gc_grid_2d
pydoc -w ESA.atmosphere.gc_grid_3d
pydoc -w ESA.atmosphere.__init__
pydoc -w ESA.atmosphere.ts_slice_m



# surface 

pydoc -w ESA.surface.convert_netcdf_flux_bpch
pydoc -w ESA.surface.convert_t3_flux
pydoc -w ESA.surface.convert_t3_map
pydoc -w ESA.surface.country_reg_def

pydoc -w ESA.surface.divide_region
pydoc -w ESA.surface.ensemble_flux_io
pydoc -w ESA.surface.extract_layer_flux
pydoc -w ESA.surface.extract_layer_flux_t3
pydoc -w ESA.surface.gen_ensemble_flux
pydoc -w ESA.surface.__init__
pydoc -w ESA.surface.sample_reg_flux
pydoc -w ESA.surface.spatial_cor
pydoc -w ESA.surface.split_t3_region
pydoc -w ESA.surface.svd_region_flux
pydoc -w ESA.surface.svd_t3_region
pydoc -w ESA.surface.write_flux_bpch

# instrument

pydoc -w ESA.instrument.ak_file_m
pydoc -w ESA.instrument.ak_m
pydoc -w ESA.instrument.aod_file_m
pydoc -w ESA.instrument.aod_m
pydoc -w ESA.instrument.cloud_file_ecmwf
pydoc -w ESA.instrument.cloud_file_m
pydoc -w ESA.instrument.cloud_m
pydoc -w ESA.instrument.err_file_m
pydoc -w ESA.instrument.err_m
pydoc -w ESA.instrument.orbit_file_m
pydoc -w ESA.instrument.orbit_m
pydoc -w ESA.instrument.read_ecmwf_cld
pydoc -w ESA.instrument.satellite_m

# enkf
pydoc -w ESA.enkf.assim_daily_obs
pydoc -w ESA.enkf.assim_def_m
pydoc -w ESA.enkf.construct_state_vector
pydoc -w ESA.enkf.etkf_cor
pydoc -w ESA.enkf.etkf_half
pydoc -w ESA.enkf.__init__
pydoc -w ESA.enkf.read_em_cfg
pydoc -w ESA.enkf.run_desc_file_m
pydoc -w ESA.enkf.run_desc_m
pydoc -w ESA.enkf.state_vector
pydoc -w ESA.enkf.x2flux_m

mv *.html ./doc
