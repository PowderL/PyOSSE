"""
This directory containes modules for representation of atmospheric compositions  at model grids. 

Modules:
----------------------------------
1. gp_grid_m.py:  general class for high-dimension grids
2. ctm_grid_m.py: general class for 3D or 4D (longitude, latitude, z, time) 
grid used to grid atmospheric compositions. 
3. gc_grid_m.py:  derived class from ctm_grid_cl for grid GEOS-Chem models 
4. gc_field_m.py: derived class from gp_field_cl for GEOS-Chem components
5. bpch_io_bufr.py: access to model bpch2 data set


Imported modules and shared libraries
-----------------------------------------
A. /general/
1.numpy
2.pylab

B ./util/
1. time_module: time conversion
2. geo_constants: constants
3. gp_axis_m.py:  general class for axis 
4. h_interpolate_m: horizontal interpolation
5. vert_interpolate_m: vertical interpolation and integration
7.gp_grid_m:general class for grid definitions 
8.gp_field_m:general class for gridded data
 

C./IO/
1. netCDF_gen.py:  code for  reading /writing netCDF files. 
2. bpch2_rw_py.so (bpch2_rw_v2.py.f90): code for reading /writing bpch files

"""
