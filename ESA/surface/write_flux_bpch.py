""" Functions for write file to GEOS-Chem BPCH2 file 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.10.12
    History: v0.95, 2013.02.01
    

    Functions:
    ==========================================================
    1. open_flux_bpch2_file: open file to write
    2. write_flux_record: save flux data into bpch2 file
    3. close_flux_bpch2_file: close bpch2 file 
    4. write_ens_flux_ml: save ensemble fluxes to bpch2 file 
    5. write_ens_map_ml: save map for ensemble fluxes to bpch2 file 
    
    
                    
   
    """
# import module  

import pylab as plb  # plot package
import numpy as npy

import ESA.util.bpch2_rw_py as bp  
import ESA.util.time_module as tm      # time format conversion 
import ESA.util.gen_plots as gpl       # interface to map and contour plots 
import ESA.atmosphere.compute_gc_grid as cmgrd 



    
def open_flux_bpch2_file(flnm, funit, title='flux'):
    """
    open the flux file to write 
    
    Inputs:
    -------------------------------------------
    1. flnm:<str>: file name
    2. funit:<int>: file unit
    3. title:<str>: title for flux data
    
    Returns:
    --------------------------------------------
    1. funit:<int>: file unit
    """
    
    stat=bp.open_bpch2_for_write(funit,full_flnm, title)
    
    if (stat<>0):
        print 'wrong in open flux '+full_flnm
    
    return funit  

def write_flux_record(funit, flux, \
                          lon, lat, \
                          tau0, tau1,\
                          ntracer, \
                          modelname,\
                          category='CO2_FLUX',\
                          unit="mole/s", \
                          ifirst=1,\
                          jfirst=1,\
                          lfirst=1,\
                          centre180=None,\
                          lonres=None,\
                          latres=None\
                          ):
    
    


    """
    Save emission to bpch2 file
    
    Inputs:
    ============================================
    1.funit:<int>: file unit
    2.flux:<array, (nlon, nlat, [nlayer])>: emission flux
    3. lon:<array, (nlon)>: longitude  grid
    4. lat:<array, (nlat)>: latitude  grid
    5. tau0:<float>:  starting time 
    6. tau1:<float>:  end time 
    7. ntracer:<int>: tracer id
    8. modelname:<str>: the source of the data
    9. category:<str>: category of the data 
    10. unit:<str>: unit of the flux
    11. ifirst: <int>: first longitude location 
    12. jfirst: <int>: first latitude location
    13. lfirst: <int>: first vertical location 
    14. lonres: <float>: longitude resolution 
    14. latres: <float>: latitude resolution 
    15. centre180:<int>: central location of longitude grid 
    
    
    Returns:
    ========================================
    1. stat:<int>: status 
    
    """
    
    
     
    # S1: information on grid guessed from lon and lat input 
    
    
    if (lonres==None):
        lonres=lon[2]-lon[1] 
    
    if (latres==None):
        latres=lat[2]-lat[1]
    
    if (centre180==None):
        # guess centre180
        if ((lon[0]>=-190) & (lon[-1]<=190)):
            centre180=1
        else:
            centre180=0
            
    # S2: half polar 
    
    halfpolar=bp.get_halfpolar()
    
    
    dims=npy.shape(flux)
    reserved="9999999"
    
    # write the data  
    
    stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                   lonres,latres,halfpolar,centre180,\
                                   ntracer,unit,tau0,tau1,\
                                   ifirst,jfirst,lfirst,flux)
    return stat



    
def close_flux_bpch2_file(funit):
    stat=bp.close_bpch2_file(funit)
    return stat



def write_ens_map_ml(flnm, \
                         regmap, \
                         lon, lat, nz,\
                         modelname='geos5',\
                         category='LANDMAP'):
    
    
    """
    Save map for ensemble  to bpch2 file 
    
    Inputs:
    ============================================
    1. flnm:<str>: file name to be save 
    2. regmap:<array, (nlon, nlat, [nlayer])>: emission map
    3. lon:<array, (nlon)>: longitude  grid
    4. lat:<array, (nlat)>: latitude  grid
    5. nz:<int>: nlayer 
    6. category:<str>: category of the data 
    
    
    """
    
    
    
    # S1: write output bpch2 file
    ext1=modelname
    nx=npy.size(lon)
    ny=npy.size(lat)
    ext2=cmgrd.get_model_resolution(nx, ny)
    
    
    if ('666' in ext2):
        # special treatment for 0.5x0.666 data
        ext2=ext2.replace('.', '')
    
    full_flnm=flnm+"."+ext1+"."+ext2
    new_regmap=1.0*regmap
    funit=39
    
    title='geos5'
    
    
    stat=bp.open_bpch2_for_write(funit,full_flnm, title)
    

    # S2: set up grid from lon and lat inputs
    
    
    lonres=lon[2]-lon[1] # the first one is different from the old-one
    latres=lat[2]-lat[1]
    
    if ((lon[0]<-190) or (lon[-1]>190)):
        centre180=0
    else:
        centre180=1
    
        
    
    
    halfpolar=bp.get_halfpolar()
    
    
    ifirst, jfirst, lfirst=1, 1, 1
    
    reserved="9999999"
    centre180=1
    modelname=bp.get_modelname()
    ntracer=1
    tau0=0.0
    tau1=24.0
    unit='unitless'
    # S3: save data to output file 
    stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                   lonres,latres,halfpolar,centre180,\
                                   ntracer,unit,tau0,tau1,\
                                   ifirst,jfirst,lfirst,new_regmap)
    
    # S4: close file 
    
    bp.close_bpch2_file(funit)


def write_ens_flux_ml(fname, regflux, \
                          lon, lat, nz,\
                          tau0, tau1,\
                          ntracer, \
                          modelname='geos5',\
                          category='CO2_FLUX',\
                          unit="mole/s"):



    """
    Save map for ensemble  to bpch2 file 
    
    Inputs:
    ============================================
    1. flnm:<str>: file name to be save 
    2. regflux:<array, (nlon, nlat, [nlayer])>: emission ensemble 
    3. lon:<array, (nlon)>: longitude  grid
    4. lat:<array, (nlat)>: latitude  grid
    5. nz:<int>: nlayer 
    6. category:<str>: category of the data 
    7. unit:<str>: unit of the data 

    Returns:
    ============================
    1. stat:<int>: status 
    
    """
    
    
    
    # S1: open output file 
    ext1=modelname
    nx=npy.size(lon)
    ny=npy.size(lat)
    ext2=cmgrd.get_model_resolution(nx, ny)
    
    
    
    if ('666' in ext2):
        # special treatment for 0.5x0.666 data
        ext2=ext2.replace('.', '')
    

    full_flnm=fname+"."+ext1+"."+ext2

    funit=41
    title=modelname
    
    stat=bp.open_bpch2_for_write(funit,full_flnm, title)
    
    # S2: set up grip

    lonres=lon[2]-lon[1] # the first one is different from the old-one
    latres=lat[2]-lat[1]
    
    if ((lon[0]<0) & (lon[-1]>0)):
        centre180=1
    else:
        centre180=0
    
        
    halfpolar=bp.get_halfpolar()
    ifirst, jfirst, lfirst=1, 1, 1
    
    reserved="9999999"
    
    # S3: save data 
    
    stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                   lonres,latres,halfpolar,centre180,\
                                   ntracer,unit,tau0,tau1,\
                                   ifirst,jfirst,lfirst,regflux)
    
    # S4: close file 
    
    bp.close_bpch2_file(funit)

    return stat

    


        
    






























    
              
    
        
