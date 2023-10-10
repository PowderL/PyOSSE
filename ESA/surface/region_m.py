"""
    class for stored region information 
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.03.24
    
    class: 
    ===========================================
    region_cl: class for regional information 
        
"""
class region_bf_cl:
    """
    container for regional flux  
    """
    
    def __init__(self, name, stype, \
                     reg_id,\
                     pid,\
                     lon, lat, area, \
                     region_map,  \
                     flux, \
                     err,\
                     err_div=[], \
                     div_map=[],\
                     **keywords):
        
        
        """
        Inputs:
        --------------------------------------
        1.  name:<str>: Name of the region 
        2.  stype:<str>:  Surface type 
        3.  reg_id:<str>: Region ID  
        4.  pid:<str>:    Parent ID 
        5.  lon:<array>:  longitude 
        6.  lat:<array>:  latitude 
        7.  area:<array>:  area 
        8.  region_map:<array>:  regional map 
        9.  flux:<array>:        fluxes
        10. err:<array>:         
        11. err_div:<list>: # 
        12. div_map:<list>: 
        **keywords
        """

        
        name=name
        self.stype=stype
        self.id=reg_id
        self.pid=pid
        
        self.lon=self.lon
        self.lat=self.lat 
        self.map=region_map
        self.area=area
        self.flux=flux
        self.err=err
        self.err_div=err_div
        self.div_map=div_map
        
        
        """
        
        self.name=name
        self.stype=stype
        self.id=reg_id
        self.pid=pid
        
        self.lon=self.lon
        self.lat=self.lat 
        self.map=region_map
        self.area=area
        self.flux=flux
        self.err=err
        self.err_div=err_div
        self.div_map=div_map
        

        
    def add_err_div(self, err_div_lst, err_map_lst):
        """ """
        
        self.err_div=self.err_div+err_div_lst
        self.div_map=self.div_map+err_div_lst
        
        
        

def read_bf_to_region_list(path, flnm, \
                               varnames=def_reg_flux_varnames, \
                               varname_dict={}, \
                               yyyy, mm, dd,\
                               **keywords):
    """
    
    # read BF from netCDF file into regional list 
    
    Members:
    ===================================================
    1. path:<str>: directory
    2. flnm: <str>: file name
    3. varnames:<>: number of variables in flux file 
    4. varname_dict:<dict>: dictionary translate varname to names in file 
        
    5. keywords:<dict>: extra inputs
    
    
    """
    # S1:  construct fdesc
    
    varnm_lst=list(varnames)
    nvar=len(varnm_lst)
    vattr=['units']
    attr_read_lst=nvar*[vattr]
    
    var_lis=vario.make_var_read_list(varnm_lst, attr_read_lst)
    
    # #c: fdesc is for file access 
    
    fdesc=ncfio.ncfile_desc_cl(var_lst=var_lst, \
                                   varname_dict=varname_dict,\
                                   read_to_array=True)
    
    fdesc.set_file_path(self.datapath)
    fdesc.set_filename_format(self.flnm)
    
    doy=tm.day_of_year(yyyy, mm, dd)

    # S2: read map and flux  into netcdf file 
    
    full_flnm=fdesc.construct_filename(XYYYYX=yyyy, XMMX=mm, XDDX=dd, XDOYX=doy)
    
    fdesc.read_var_from_file(full_flnm)
    lon=fdesc.get_lon()
    lat=fdesc.get_lat()
    parent_id=fdesc.get_data('parent_id')
    pid=fdesc.get_data('pid')
    
    parent_flux=fdesc.get_data('parent_flux')
    parent_map=fdesc.get_data('parent_map')
    parent_stype=fdesc.get_data('parent_stype')
    area=fdesc.get_data('area')
    nreg=len(parent_id)
    # creat parent 
    
    parent_name=ncfio.get_attr(full_flnm, '', 'parent_name')
    
    
    if (parent_name==None):
        parent_name=nreg*['']
    
    else:
        parent_name=parent_name.split(';')
        
    
    # S3: create a region member for each sub regions
    
    reg_id_lst=[]
    reg_cl_lst=[]
    
    for ireg in nreg:
    
        stype=parent_stype[ireg]
        reg_map=parent_map[:,:,ireg]
        flux=parent_flux[:,:,ireg]
        name=parent_name[ireg]
        reg_id=parent_id[ireg]
        
        region_bf_cl(name, stype, \
                         reg_id,\
                         pid,\
                         lon, lat, area, \
                         reg_map,  \
                         flux, \
                         err_div=[],\
                         err_map=[]\
                         )
        
        
        reg_cl_lst.append(reg_id)
        reg_id_lst.append(reg_id)
        
        
        
        
    
                 

