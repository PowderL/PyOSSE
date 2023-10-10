""" class for reading and writing binary punch files
    This is class is stand-alone version for otool_gcfile_io
    
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2012.09.30
    
    Classes:
    =========================================
    1. bpgrid : grid for bpch2 data 
    2. bpch2_data : container for bpch2 data
    3. tracer_info: class for tracerinfo file 
    4.  diag_info:  class for diaginfo file 
    5. bpch2_file_rw: class for functions to read and write BPCH2 data
    
      
"""

from numpy import *

import bpch2_rw_py 
import ESA.util.time_module as tm
import ESA.util.geo_constant as gc
import ESA.util.gp_axis_m as gax
import ESA.atmosphere.compute_gc_grid as cmg


bpch2_fill_val=-999
bpch2_fill_str_val="???"
bpch2_fill_flt_val=-999.0
def_tracerinfo_file='tracerinfo.dat'
def_diaginfo_file='diaginfo.dat'

class bpgrid:
    
    """ the class for the grid used to store GEO-CHEM binary punch data

    Members:
    --------------------------------------------------
    1. ix:<int>: size of longitude
    2. jx:<int>: size of latitude
    3. lx:<int>: size of model level
    4. mod_res:<str>:model resolution  ('4x5', '2x2.5' etc)
    5. ifirst:<int>: starting longitude
    6. jfirst:<int>: starting latitude
    7. lfirst:<int>: starting vertical level
    8. halfpolar:<int>:  flag for half polar data set
    9. centre180:<int>: location of longitude centre
    10. glb_ix:<int>: global longitude size
    11. glb_jx:<int>: global latitude size
    12. glb_lx:<int>: global vertical size
    13. lonres:<float>: longitude resolution
    14. latres:<float>: latitude resolution
    15. lons:<array>: longitude grid
    16. lats:<array>: latitude grid 
    
    Functions:
    -----------------------------------------------------------
    1. __init__: initialization
    2. get_lat: Get the latitude grid 
    3. get_lon:  Get the longitude grid 
    4. get: Get vertical level
    
    

    """    
    def __init__(self, ix, jx, lx, \
		     ifirst=1, jfirst=1, lfirst=1,\
		     halfpolar=0, centre180=0, \
		     lonres=None, latres=None):
        """ Initialize the class 
        Inputs:
        ============================================================
        1. ix, jx, lx: 	<int>:	the sizes of lon, lat, z 
        2. ifirst, jfirst, lfirst:<int>:  the start point for the data to store 
        3. halfpolar, centre180:<int>: 0 or 1 the type of the longitude grid
        4. lonres, latres: <float>: the resolution of the lat and lon 
        
        Notes:
        =============================================================
        1. The user can include the details on the vertical grid as well.    	
        
        """
        # S1: model information 
        
        self.ix=ix
        self.jx=jx
        self.lx=lx
        self.mod_res=""
        self.ifirst=ifirst
        self.jfirst=jfirst
        self.lfirst=lfirst
        self.halfpolar=halfpolar
        self.centre180=centre180
        self.glb_ix=ix
        self.glb_jx=jx
        self.glb_lx=lx
        
        
        # S2: set resolution 
        
        # #T: longitude 
        
        if (lonres==None):
            self.lonres=360.0/(self.glb_ix)
        else:
            self.lonres=lonres
            self.glb_ix=int(360.0/lonres)
            
        if (self.glb_ix>=539):
            self.glb_ix=540
            
        # #T: latitude

        if (latres==None):
            self.latres=180.0/(self.glb_jx-1)
        else:
            self.latres=latres
            self.glb_jx=int(180.0/latres)+1
            if (self.glb_jx>361):
                self.glb_jx=361
                
                
        self.mod_res=cmg.get_model_resolution(self.glb_ix, self.glb_jx)
        
        # S3: set model grid 
        
        self.lats=cmg.get_model_lat(self.mod_res, self.glb_jx)
        self.lons=cmg.get_model_lon(self.mod_res, self.glb_ix)
        
        self.zs=self.get_z()
		
    def get_lat(self):
        """
        Get the latitude grid 
        
        Return:
        ================================
        1. sel_lat:<array>: latitude grid 
        
        """
        # the middle point
        
        lat=cmg.get_model_lat(self.mod_res, self.glb_jx)
        sel_lat=lat[self.jfirst-1:self.jfirst+self.jx-1]
            
        return sel_lat
        
           
    def get_lon(self):
            
        """
        Get the longitude grid 
        
        Return:
        ================================
        1. sel_lon:<array>: longitude grid
        
        """
            
        if ((self.ix+self.ifirst-1)<=self.glb_ix):
            sel_lon=self.lons[self.ifirst-1:self.ifirst+self.ix-1]
        else:
            
            sel_lon1=self.lons[self.ifirst-1:self.glb_ix]
            sel_lon2=self.lons[0:self.ifirst-1+self.ix-self.glb_ix]
            sel_lon=concatenate((sel_lon1, sel_lon2))
            
            
        return sel_lon
	
			
			
    def get_z(self):
        """
        Get vertical level 
        
        Return:
        ================================
        1. level:<array>: vertical level
        """
            
        
        level=arange(self.lx)
        return level

		
                 
        
class bpch2_data:
    """
    container for bpch2 data set 
    
    Members:
    ------------------------------------------------
    1. ntracer:<integer>: the order number of the tracer in GEOS-CHEM 
    2. category:  <str>: category of the tracer 
    3. unit:      <str>: unit of the data 
    4. bpgrid:  <bpgrid>:class for data grid 
    5  attr: <dict>: attributes

    Functions
    ----------------------------------------------------------
    1. __init__: initialization 
    2. write_to_bpch2_file_window: write window data set into one bpch2 file
    3. write_to_bpch2_file: write data set into one bpch2 file
    4. regrid: regrid data 
    5. def set_attr: set one attribute to the data set 
    6. get_attr: get attributes of the data set  
    7. update_data: update data 
    
    """
    def __init__(self, ntracer, category, unit, bpgrid, data, **keywords):
        
        """ 
        Initialization 
        
        Inputs:
        ------------------------------------------------------
        1. ntracer:<integer>: the order number of the tracer in GEOS-CHEM 
        2. category:  <str>: category of the tracer 
        3. unit:      <str>: unit of the data 
        4. bpgrid:  <bpgrid>:class for data grid 
        5  keyword: <dict>: additional information
        which given in form of  tau0=a, tau1=b, modelname='md', reserved='a'
		"""
        self.ntracer=ntracer  # the tracer no in model
        self.grid=bpgrid
        self.data=data
        self.category=category
        self.unit=unit
        
        self.attr={}
        if (len(keywords)>0):
            self.attr.update(keywords)
        
    def write_to_bpch2_file_window(self, funit, ist, jst, nx, ny):
        """ write window data set into one bpch2 file 
        
        Inputs:
        ==============================================
        1. funit: <int>: the unit number of the bpch2 file 
        to be consistent with  GEOS-Chem Fortran index,
        2) ist:<int>: starting poistion of longitude
        3) jst:<int: starting poistion of latitude 
        4. nx: <int>: size of longitude 
        5. ny: <int>: size of latitude
        
        Returns:
        -------------------------------------------
        1. stats:<int>: IO status
        
        """
	        
        # S1: information on model grid
        
        ix=nx
        jx=ny	
        lx=self.grid.lx
        
        
        lonres=self.grid.lonres
        latres=self.grid.latres
        
        # ifirst=self.grid.ifirst
        # jfirst=self.grid.jfirst
        # lfirst=self.grid.lfirst
        
        ifirst=ist
        jfirst=jst
        
        lfirst=self.grid.lfirst
        
        halfpolar=self.grid.halfpolar
        centre180=self.grid.centre180
        
        ntracer=self.ntracer

        # S2: data slice to be saved 
        
        data=self.data
        data=data[ifirst-1:ifirst-1+ix, jfirst-1:jfirst-1+jx, :]
        
        # S3: information on tracer 
        
        category=self.category
        unit=self.unit
        tau0=0.0
        tau1=0.0
        modelname=" "
        reserved=" "
        modelname=" "	
        
        if ('modelname' in self.attr):
            modelname=self.attr['modelname']
        # S5: time and unit 
        
        if ('tau0' in self.attr):
            tau0=self.attr['tau0']
            
        if ('tau1' in self.attr):
            tau1=self.attr['tau1']
		
            
        if ('reserved' in self.attr):
            reserved=self.attr['reserved']
                
        reserved=reserved.strip()
        
        if (len(reserved)==0):
            reserved="000000"
        

        traunit=self.attr['traunit']

        # S6: save data 
        
        stat = bpch2_rw_py.write_bpch2_data(funit,modelname,category,reserved, \
                                                lonres,latres,halfpolar,centre180,\
                                                ntracer,unit,tau0,tau1,\
                                                ifirst,jfirst,lfirst,data)
        return stat
            

    def write_to_bpch2_file(self, funit):
        
        """ write data into one open bpch2 file 
        
        Inputs:
        =============================================
        
        1. funit: <int>: unit number of the bpch2 file 
        
        Returns:
        -----------------------------------------------
        1. stat:<int>: IO status 
        
        """
        
        # S1: model grid 
        
        ix=self.grid.ix
        jx=self.grid.jx	
        lx=self.grid.lx
        
        lonres=self.grid.lonres
        latres=self.grid.latres
        ifirst=self.grid.ifirst
        jfirst=self.grid.jfirst
        lfirst=self.grid.lfirst
        halfpolar=self.grid.halfpolar
        centre180=self.grid.centre180
    
        # S2: information on tracer 
        
        ntracer=self.ntracer
        data=self.data
        category=self.category
        unit=self.unit
        
    
        modelname=" "
        reserved=" "
        modelname=" "	
	

        if ('modelname' in self.attr):
            modelname=self.attr['modelname']


        # S3: time and unit 

        tau0, tau1=0.0, 0.0
        
        if ('tau0' in self.attr):
            tau0=self.attr['tau0']
            
        if ('tau1' in self.attr):
            tau1=self.attr['tau1']
            
        if ('reserved' in self.attr):
            reserved=self.attr['reserved']
            
        
        reserved=reserved.strip()
        
        if (len(reserved)==0):
            reserved="000000"
        

        
        traunit=self.attr['traunit']
        # print unit, traunit
        
        # S4: save data 
        
        stat = bpch2_rw_py.write_bpch2_data(funit,modelname,category,reserved, \
                                                lonres,latres,halfpolar,centre180,\
                                                ntracer,unit,tau0,tau1,\
                                                ifirst,jfirst,lfirst,data)
        return stat
    
	
    def regrid_data(self, new_lon, new_lat):
        """ regrid data 
        
        Inputs:
        ------------------------------------------------------
        1. new_lon:<array>: new longitude grid 
        2. new_lat:<array>: new latitude grid 
        
        
        """
        
        ix=self.grid.ix
        jx=self.grid.jx	
        lx=self.grid.lx
		
        lonres=self.grid.lonres
        latres=self.grid.latres
        ifirst=self.grid.ifirst
        jfirst=self.grid.jfirst
        lfirst=self.grid.lfirst
        halfpolar=self.grid.halfpolar
        centre180=self.grid.centre180
        
        ntracer=self.ntracer
        data=self.data
        category=self.category
        
        lat=self.grid.get_lat()
        lon=self.grid.get_lon()
        ax_lat=gax.gp_axis('lat', lat)
        ax_lon=gax.gp_axis('lat', lon)
        lonp1, lonp2, lonw=ax_lon.getwgt(new_lon)
        latp1, latp2, latw=ax_lat.getwgt(new_lat)
        new_ix=size(new_lon)
        new_jx=size(new_lat)
        new_lx=lx
        
        new_data1=lonw[:,newaxis, newaxis]*data[lonp1, :, :]+\
            (1.0-lonw[:,newaxis, newaxis])*data[lonp2,:,:]
        new_data=latw[newaxis, :, newaxis]*new_data1[:,latp1, :]+\
            (1.0-latw[newaxis, :, newaxis])*new_data1[:,latp2, :]
        
        self.data=array(new_data)
        new_grid=bpgrid(new_ix, new_jx, new_lx)
        
        self.grid=new_grid
	
		
    def set_attr(self, attr_name, value):
        """ set one attribute to the data set 
        
        Inputs:
        ----------------------------------
        1. attr_name:<str>:name of the attribute
        2. value: <obj>: value of the attribute 
        
        """	      
        self.attr.update({attr_name:value})
        
    def get_attr(self, attr_names):
        """ get attributes of the data set  
        Inputs:
        --------------------------------------------
        1. attr_names:<str/list>: name of the attributes to be fetched. 
        
        """        
        if (attr_names==None):
            return self.attr	
		
        at=list()
        ats=list()
        if (type(attr_names)==type(at)):
            at=attr_names
        else:
            at.append(attr_names)
            
        for attr_name in at:
            if (attr_name in self.attr):
                ats.append(self.attr[attr_name])
            else:
                print 'no attribute '+attr_name +'in data '+self.category 
                return None
			
        
        return ats
	
    def update_data(self,data):
        """
        update the data set 
        
        Inputs:
        -------------------------------------
        1. data:<array>: data
        
        """
        
        if (self.data<>None):
            del self.data
        
        self.data=data
	
	
    def display(self, dirc, **keywords):
        
        """ display the data over a map 
        
        
        Inputs:
        --------------------------------------------
        dirc:<int>: action type (reserved)
        
        keywords:<dict>: options for display
        
        """
	import gen_plots as gpl
	
		
        rlat=self.grid.get_lat()
        rlon=self.grid.get_lon()
        levels=self.grid.get_z()
        
        stitle=""
        if ('title' in keywords):
            stitle=keywords['title']

        if ('name' in self.attr):
            add_str=self.attr['name']
        else:
            add_str=self.category
            
        add_str=add_str.strip()
        stitle=stitle+' '+add_str
        
        if ('traunit' in self.attr):
            add_str=self.attr['traunit']
        else:
                
            add_str=self.unit
        
        add_str=add_str.strip()
        
        stitle=stitle+' ('+add_str+')'
        
        keywords.update({'title':stitle})
        
                        
        if ('tau0' in self.attr):
            tau0=self.attr['tau0']
            tau0=3600.*tau0 # convert to seconds
            utc=tm.tai85_to_utc(tau0)
            stitle=stitle+' '+utc
                
        level=0
        if ('level' in keywords):
            level=keywords['level']
            del keywords['level']
        
        factor=1.0
        
        if ('factor' in keywords):
            factor=keywords['factor']
            del keywords['factor']
        
        if ('show_map' in keywords):
            show_map=keywords['show_map']
            del keywords['show_map']  
        
        vals=self.data[:,:,level]
        
        vals=array(vals)
        vals=factor*vals

        gpl.plot_map(vals, rlon, rlat, **keywords)
        
        gpl.show_plot()
        
        
        return 0
    
	
        
class tracer_info:
    """ class for tracer info file 
    see otool_gcfile_io.py for more information

    """
    def __init__(self, flnm):
        
        """
        Initialization 
        Inputs:
        ------------------------------
        1. flnm:<str>: name of the tracerinfo file 

        """
        self.name=list()
        self.fullname=list()
        self.tracer=list()
        self.scale=list()
        self.unit=list()
        self.c=list()
        self.molew=list()
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                print 'Tracer info file '+flnm+' not found'
                file_found=False
		
        if (file_found):
            lines=fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0]<>'#'):
                    sname, sfullname, smw, sc, stra, sscal, sunit= \
                        line[0:8], line[9:39], line[39:49], line[49:53], line[53:62], line[62:72], line[73:] 
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.name.append(sname.strip())
                    self.fullname.append(sfullname.strip())
                    self.molew.append(float(smw))
                    self.c.append(float(sc))
                    self.scale.append(float(sc))
                    self.tracer.append(int(stra))
                    self.unit.append(sunit.strip())
		
            # loop line end 
                    
            self.tracer=array(self.tracer)
		

    def get_tracer_info(self, tracer_in):
        """
        Get information for tracer
        
        Inputs:
        ----------------------------------------------
        1. tracer_in:<int>: tracer ID
        
        Returns:
        ------------------------------------------------
        1. name:<str>: tracer name
        2. fullname:<str>: tracer full name
        3. molew: <float>: moleculr mass
        4. unit: <str>: unit of the tracer
            
        """
        if (len(self.tracer)==0):
            # if (no tracer found)
            return bpch2_fill_str_val, bpch2_fill_str_val, \
                bpch2_fill_flt_val, bpch2_fill_flt_val, bpch2_fill_str_val
        
		
        idx=where(tracer_in==self.tracer)
        # idx=compress(idx)
        if (size(idx)>=1):
            name=self.name[idx[0]]
            fullname=self.fullname[idx[0]]
            scale=self.scale[idx[0]]
            unit=self.unit[idx[0]]
            c=self.c[idx[0]]
            molew=self.molew[idx[0]]
            return name, fullname, molew, scale, unit 
        else:
            return bpch2_fill_str_val, bpch2_fill_str_val, \
                bpch2_fill_flt_val, bpch2_fill_flt_val, bpch2_fill_str_val
        
    def load_tracer_info(self, flnm):
        """
        read tracer infor from tracerinfo file 
        
        Inputs
        --------------------------------------
        1. flnm:<str>: file name 
        
        """
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                print 'No tracer info file: ', flnm.strip() 
                file_found=False
                    
        if (file_found):
            lines=fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0]<>'#'):
                    
                    sname, sfullname, smw, sc, stra, sscal, sunit= \
                        line[0:8], line[9:39], line[39:49], line[49:53], line[53:62], line[62:72], line[73:] 
                        # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.name.append(sname.strip())
                    self.fullname.append(sfullname.strip())
                    self.molew.append(float(smw))
                    self.c.append(float(sc))
                    self.scale.append(float(sc))
                    self.tracer.append(int(stra))
                    self.unit.append(sunit.strip())
                    
            # loop line end 
                    
            self.tracer=array(self.tracer)
    
		

				
class diag_info:
    """class for diag info file 
    
    see otool_gcfile_io.py for more information

    """
    def __init__(self, flnm):
        
        """
        Initialization  
        
        Inputs:
        ----------------------------------------
        1. flnm:<str>: name of the file 
        
        2. 
        
        """
        self.category=list()
        self.comment=list()
        self.offset=list()
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                print 'Diag info file '+flnm+' not found'
                file_found=False
        if (file_found):
            lines=fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0]<>'#'):
                    soffset, scategory, scomment= line[0:8], line[9:49], line[49:] 
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.offset.append(int(soffset))
                    self.category.append(scategory.strip())
                    self.comment.append(scomment.strip())
            # loop line end 
                        

    def get_offset(self, category):
        """ get offset for given  category
        
        Inputs:
        ------------------------------------
        1. category:<str>: tracer category 
        
        Returns:
        --------------------------------------
        1. offset:<int>: offset of the category 
        
        """
        
        nid=len(self.category)
        found=False
        scate=category.strip()
        for id in range(nid):
            if (scate in self.category[id]):
                found=True
                offset=self.offset[id]
        
        
        if (found):
            return offset
        else:
            return bpch2_fill_flt_val
        
    def load_diag_info(self, flnm):
        """
        Read diag info from file 
        
        """
        file_found=False
        if (flnm.strip()<>""):
            try:
                fl=open(flnm.strip(), "r")
                file_found=True
            except IOError:
                print 'no diag info file: ',  flnm.strip()
                file_found=False

        if (file_found):
            lines=fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0]<>'#'):
                    soffset, scategory, scomment= line[0:8], line[9:49], line[49:] 
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.offset.append(int(soffset))
                    self.category.append(scategory.strip())
                    self.comment.append(scomment.strip())
                    
                    
        
            # loop line end 

                    			
class bpch2_file_rw:
    
    """ the class for functions read and write BPCH2 data 
    
    Members:
    ---------------------------------------------------
    1.flnm:<str>: name of bpch2 flnm
    2. cur_tracer:<int>: index for current tracer
    3. maxtracer:<int>: maximal tracer number 
    4. self.title:<str>: title of the bpch2 file 
    5. modelname:<str>: name of GEOS-Chem simulation 
    6. data:<list>: list of bpch2 data set read from file 
    7. stat:<int>: IO status
    8. ntracers:<int>: number of tracers read in 
    9. read_mode:<int>: modes (-1, 0,1) for file access 
    10. tracerinfo:<tracer_info>: class for information on tracers 
    11. diaginfo:<diag_info>: class for diaginfo 
    

    Functions:
    -------------------------------------------------------
    1. __init__: initialize
    2. open_bpch2_w: open bpch2 file for outputs
    3. close_bpch2_file: close bpch2 file
    4. write_bpch2_data: write selected data to  bpch2 file
    5. add_bpch2_data: add a bpdata into self.data list    
    6. get_data: fetch bpch2 data according specifications
    7. read_bpch2_data: read in single pch2 data record
    8. print_datainfo: print out info on selected bpch2 data
    9. write_bpch2_head: write head to one bpch2
    10. open_bpch2_for_write: one wrapper for open_bpch2_w
    11. read_mode_0: Read tracer information and head from bpch2 file only
    12. read_mode_1: Read all tracer data from bpch2 file
    13. read_mode_2: read only selected data from bpch2 file 
    
    """	
    def __init__(self, flnm, mode, maxtracer=300, \
                     tmpflnm=None,\
                     do_read=0, \
                     ftracerinfo=def_tracerinfo_file, \
		     fdiaginfo=def_diaginfo_file,\
		     sel_categorys=None, sel_tracers=None, sel_taus=None ):
 		
        """	initialize
        Inputs:
        -------------------------------------------------------
        1. flnm:<str>: bpch file name
        2. mode:<str>: read and write bpch2 data
        3. maxtracer:<int>: the max number of tracers
        4. tmpflnm: <string>: file name for  temperary saving 
        5. do_read: <int>: if do_read=1, data itself will be read in instead of just its information 
        6. ftracerinfo:<str>: name of tracer info file 
        7. fdiagfo:<str>: name of diag info file 
        8. sel_categorys:<list, t:str>: selected category 
        9. sel_tracers: <list, t:int>: selected tracer
        10. sel_taus: <list, t:float>: selected tau (time)
        
        """

        # S1: initialize class member 
        
        self.flnm=flnm
        self.cur_tracer=-1
        self.maxtracer=maxtracer
        self.title=None
        self.modelname=None
        
        self.data=list()
        self.stat=-1
        self.ntracers=0
        self.read_mode=-1
        
# grid 

        # S2: initialize diag and tracer info 
        self.tracerinfo=tracer_info(ftracerinfo)
        self.diaginfo=diag_info(fdiaginfo)
		
        
        if (tmpflnm==None):        	
            flnm_tmp=flnm.strip()+'_info_py'
        else:
            flnm_tmp=tmpflnm.strip()
        
        # S3: read bpch2 file  
        if (mode=='r'):
            if (do_read==0):
                # #c: read head 
                self.read_mode_0(flnm, flnm_tmp, maxtracer)
            elif (do_read==1):
                self.read_mode_1(flnm)
            elif (do_read==2):
                self.read_mode_2(flnm, sel_categorys, sel_tracers, sel_taus)
        
				
    def open_bpch2_w(self, funit, title=""):
        """
        open bpch2 file for outputs
        
        Inputs:
        ---------------------------------------
        1. funit:<int>: ID for file IO 
        2. titel:<str>: title of the data set 
        
        Returns:
        ----------------------------------------
        1. status:<int>: status of the file  
        
        """
        
        self.funit=funit
        print 'funit', funit
        stat=bpch2_rw_py.open_bpch2_for_write(funit,self.flnm, title)
        self.stat=stat
        return stat
	
    def close_bpch2_file(self):
        """
        close bpch2 file
        """

        if (self.funit<>None):
            stat=bpch2_rw_py.close_bpch2_file(self.funit)
            self.stat=stat
            self.funit=None
            
    def write_bpch2_data(self, data_id):
        """
        write selected data to  a bpch2 file
        
        Inputs:
        -------------------------------------
        1. data_id:<int>: index for data to be written
        
       
        """

        stat=-1
        if (data_id<self.ntracers and self.funit<>None):
            bpdata=self.data[data_id]
            stat=bpdata.write_to_bpch2_file(self.funit)
        
        return stat
	
    def add_bpch2_data(self,  bpdata, data_id=None):
        
        """
        add a bpdata into self.data list 
        
        Inputs:
        -------------------------------------
        1. bpdata:<bpch2_data>: bpch2 data class to be added to the list 
        2. data_id:<int>: data ID 
        
        
                
        """

        if (data_id==None):
            
            self.data.append(bpdata)
            self.ntracers=self.ntracers+1
        else:
            self.data[data_id]=bpdata
        
			
    def get_data(self, categorys=None,tracers=None, \
                     taus=None, tranames=None):
        
        """ fetch bpch2 data according specifications
        
        Inputs:
        --------------------------------------------------------------
        1. categorys:<list, t:str>: list of category 
        2. taus:<list, t:float>: list of time 
        2. tracers:<list, t:int>: list of tracer ID 
        3. tranames:<list, t:str>: list of tracer names
        
        Returns:
        -------------------------------------------------------------------
        1. found_data:<list, t:bpch2_data>: list of bpch2_data meeting criteria sets (see Note 1)
        2. founded:<list, t:int>: number of bpch2 data meeting each  criteria sets (see Note 1
        
        Notes: 
        --------------------------------------------------------------------
        1. One criteria set [category, tracer,tau, traname]  is formed by one element 
        from each lof categorys, tracers, taus, and tranames lists.  
        
        """
        if (self.ntracers<=0):
            print 'tracer not found'
            return None, None
        
        # S1: check matching scores
	
        nval=0
        array_cat=None
        array_tracer=None
        array_tau=None
        array_traname=None
	
        
        if (categorys<>None):
            array_cat=array(categorys)
            nval=size(array_cat)
		
        if (tracers<>None):
            array_tracer=array(tracers)
            nval=size(array_tracer)
			
        if (taus<>None):
            array_tau=array(taus)
            nval=size(array_tau)

        if (tranames<>None):
            array_traname=array(tranames)
            nval=size(array_traname)
		
        # #c: no requirement specified. 
            
        if (nval==0):
            nval=len(self.data)
            founded=ones(nval)
            return self.data, founded
	
        
        found_data=list()
        founded=zeros(nval)
	
	# S2: check one by one for data meeting (one of) the requirement sets 
        
        for bpdata in self.data:
            
            
            score_cat=ones(nval)
            score_tracer=ones(nval)
            score_tau=ones(nval)
            score_traname=ones(nval)
            
            if (array_cat<>None):
                score_cat=where(array_cat==bpdata.category.strip(),1,0)
            #	print '1', score_cat
            if (array_tracer<>None):
                score_tracer=where(array_tracer==int(bpdata.ntracer), 1,0)
				# print '2', score_tracer
            if (array_tau<>None):
                tau0, tau1=bpdata.get_attr(['tau0', 'tau1'])
                score_tau=where(logical_and(array_tau>=tau0, array_tau<tau1), 1,0)
                
                # print '3', score_tau
                # print array_tau, tau0
				
				
            if (array_traname<>None):
                traname, tau1=bpdata.get_attr(['name', 'tau1'])
                # print 'traname', traname
                score_traname=where(array_traname==traname.strip(), 1,0)
                # print '4', score_traname
				
            final_score=score_cat*score_tracer*score_tau*score_traname
            founded=founded+final_score

            # #S3: add it to the found_data if current bpch2_data meet one criteria
            
            if (sum(final_score)>0):
               
                if (bpdata.data==None):
                    # #c: when the data has not been read in 
                    # #read in the required data 
                    ix=bpdata.grid.ix
                    jy=bpdata.grid.jx
                    lz=bpdata.grid.lx
                    
                    tracer=bpdata.ntracer
                    xtau0, xtau1=bpdata.get_attr(['tau0', 'tau1'])
                    category=bpdata.category
                    data, dunit,ios=bpch2_rw_py.read_bpch2(self.flnm, category, tracer,\
                                                               xtau0,  ix,  jy, lz)
                    print 'IOS', ios
                    
                    if (ios>0):
                        print 'error in read data', ios
                        return None, None
                    
                    elif (ios==-1):
                        print 'no data found'
                        return None, None
                    print 'after read', size(data)
                    bpdata.update_data(data)	
                # if bpdata.data end 
                    
                found_data.append(bpdata)
        # loop bpdata end 
           
        return found_data, founded
		
    def read_bpch2_data(self, data_id, save_data=False):
        """ read in single pch2 data record
        
        Inputs:
        -----------------------------
        1. data_id:<int>: data index to be read
        2. save_data:<T/F>: if true, bpchdata will be updated 

        Returns:
        ---------------------------------------
        1. data:<obj>:  if save_data==True, class bpch2_data will be return. Otherwise data (array) will be return 
        
        """
        # S1: if all data from one file have been read in 

        if (data_id > self.ntracers):
            print 'tracer not found'
            return None
        
        if (self.read_mode==1):
            bpdata=self.data[data_id]
            return bpdata
		

        # S2: if the data has not been read, read it from bpch2 file 
        
        bpdata=self.data[data_id]
        
        ix=bpdata.grid.ix
        jy=bpdata.grid.jx
        lz=bpdata.grid.lx
		
        tracer=bpdata.ntracer
        tau0, tau1=bpdata.get_attr(['tau0', 'tau1'])
        category=bpdata.category
        data, dunit, ios=bpch2_rw_py.read_bpch2(self.flnm, category, tracer,\
                                                    tau0,  ix,  jy, lz)
        
        # #c: check IO status
        if (ios>0):
            print 'error in read data', ios
            return None
        elif (ios==-1):
            print 'no data found'
            return None
        

        # S3: return data or clas  bpch2_data 
        
        if (save_data):
            bpdata.update_data(data)
            return bpdata
        else:
            return data
	
        
    def print_datainfo(self, data_id=None):
        """
        print out data info 
        
        Inputs:
        -----------------------------------------]
        1. data_id:<int>: index of data (not tracer ID) to be printed out 
        
        
        """
        # S1: print out information for selected bpch2_data

        if (data_id<>None and data_id<self.ntracers):
            
            bpdata=self.data[data_id]
            ix=bpdata.grid.ix
            jy=bpdata.grid.jx
            lz=bpdata.grid.lx
            tracer=bpdata.ntracer
            unit=bpdata.unit
            category=bpdata.category
            tau0, tau1, modelname, traname, unit, traid=bpdata.get_attr(['tau0', 'tau1', 'modelname', 'name', 'traunit', 'id'])
            unit=bpdata.unit
            print data_id, tracer, traid, traname.strip(), category.strip(), ix, jy, lz, unit.strip()
				
        # S2: or print out information fpr bpch2_data
        
        else:
            for data_id in range(self.ntracers):
                
                bpdata=self.data[data_id]
                ix=bpdata.grid.ix
                jy=bpdata.grid.jx
                lz=bpdata.grid.lx
                tracer=bpdata.ntracer
                category=bpdata.category
                
                tau0, tau1, modelname, traname, unit, traid=bpdata.get_attr(['tau0', 'tau1', 'modelname', 'name', 'traunit', 'id'])
                unit=bpdata.unit
                print data_id+1, tracer, traid, traname.strip(), category.strip(), ix, jy, lz, unit.strip()
				
    def open_bpch2_for_write(self, title, funit):
        """
        open bpch2 file for outputs
        It is a wrapper for old name open_bpch2_w

        Inputs:
        ---------------------------------------
        1. titie:<str>: title of the data set 
        2. iunit:<int>: ID for file IO 
        
        Returns:
        ----------------------------------------
        1. stat:<int>: status of the file  
        
        """

        state=bpch2_rw_py.open_bpch2_for_write(iunit,self.flnm,title)
        if (state==0):
            self.funit=iunit
            
    def write_bpch2_head(self, funit, title):
        
        """
        write head to one bpch2

        Inputs:
        ---------------------------------------
        1. funit:<int>: ID for file IO 
        2. titie:<str>: title of the data set 
        
        Returns:
        ----------------------------------------
        1. stat: status of the file 
        """
        stat=bpch2_wr_py.write_bpch2_hdr(iunit,title)
        return stat
    
    def read_mode_0(self, flnm, flnm_tmp, maxtracer):
        
        """
        Read tracer information and head from bpch2 file only
        
        Inputs:
        -------------------------------------------------------------------
        1. flnm:<str>: name of bpch2 file 
        2. flnm_tmp:<str>: file name for temporary file 
        3. maxtracer:<int>: maximal tracer 
        
        
        Returns:
        -------------------------------------------
        1. ios:<int>: IO status 
        
        """
        # #t: only the the heads will be first
        
        # S1: read bpch2 data head 

        ios,title,tracer_id, lonres, latres,\
            ix,jx,lx,\
            ifirst, jfirst, lfirst,\
            halfpolar, centre180, \
            tau0,tau1,ntracers= \
            bpch2_rw_py.read_bpch2_head(flnm,flnm_tmp, maxtracer)
                
                
        
        if (ios<>0):
            print 'error in read data', ios
            self.read_mode=-1
            self.stat=ios
            return ios
        
        # #c: read temporary file 
            
        self.title=title.strip()
        tf=open(flnm_tmp,'r')
        lines=tf.readlines()
        tf.close()
        ic=0
            
        category=list()
        modelname=list()
        unit=list()
        reserved=list()
        
        # #c: loop over temporary lines 
        
        
        if (ntracers==1):
            line=lines[0]
            terms=line.split(',')
            category.append(terms[0])
            modelname.append(terms[1])
            unit.append(terms[2])
            reserved.append(terms[3].strip())
        
        else:
            for line in lines:
                line=line.replace('\n', '')
                terms=line.split(',')
                category.append(terms[0])
                modelname.append(terms[1])
                unit.append(terms[2])
                reserved.append(terms[3])
				
        
        self.title=title
        self.stat=0
        self.data=list()
        self.ntracers=ntracers
                            
        # S3: store the data into the class bpch2_data 
                
        for isp in range(ntracers):
            # #c: create grid 
            data_grid=bpgrid(ix[isp], jx[isp], lx[isp], \
                                 ifirst[isp], jfirst[isp], lfirst[isp],\
                                 halfpolar[isp], centre180[isp], \
                                 lonres[isp], latres[isp])
            
            # #c: decode tracer ID and category ID 
            
            offset=self.diaginfo.get_offset(category[isp])
            if (offset<0):
                offset=0
                
            real_id=tracer_id[isp]+offset
            traname, trafullname, tramolew, trascale, traunit=self.tracerinfo.get_tracer_info(real_id)
            # create bpch2_data 
            
            pbdata=bpch2_data(tracer_id[isp],category[isp], unit[isp],\
                                  data_grid, None, \
                                  tau0=tau0[isp],tau1=tau1[isp], modelname=modelname[isp],\
                                  reserved=reserved[isp], offset=offset, name=traname, \
                                  fullname=trafullname, molew=tramolew, scale=trascale, \
                                  traunit=traunit, id=real_id)
            
            
            self.data.append(pbdata)
            
        
        self.stat=ios
        self.read_mode=0
        return ios
    
    def read_mode_1(self, flnm):
        
        """
        Read all tracer data from bpch2 file
        
        Inputs:
        -------------------------------------------------------------------
        1. flnm:<str>: name of bpch2 file 
        
        
        Returns:
        -------------------------------------------
        1. stat:<int>: IO status 
        """
        
        print 'read_mode_1:', self.flnm
        
        # S1: open file 
        
        funit=199
        fti,title,stat = bpch2_rw_py.open_bpch2_for_read(funit,self.flnm)
        self.title=title.strip()
        
        if (stat<>0):
            print 'error in read :',  stat
            self.stat=self.stat
            return 0

        # S2: read records till end of the file 
        
        vtracer_id,vhalfpolar,vcentre180,\
            vni,vnj,vnl,vifirst,vjfirst,vlfirst,vlonres,vlatres,\
            vtau0,vtau1,vmodelname,vcategory,vunit,\
            vreserved,vdata_array,stat = bpch2_rw_py.read_bpch2_record(funit)
        
        # loop over till 
        while (stat==0):
            # S3: create grid
            data_grid=bpgrid(vni, vnj, vnl,\
                                 vifirst, vjfirst, vlfirst,\
                                 vhalfpolar, vcentre180, \
                                 vlonres, vlatres)
            
            # S4: decode tracer and diag info
            
            offset=self.diaginfo.get_offset(vcategory)
            if (offset<0):
                offset=0
                real_id=vtracer_id+offset
            else:
                real_id=vtracer_id+offset
            
            traname, trafullname, tramolew, trascale, traunit=self.tracerinfo.get_tracer_info(real_id)
            
            # S5: fill data array 
            
            vdata=zeros([vni, vnj, vnl], float)
            # print vni, vnj,vnl
            if (vlatres<1.0):
                # print vifirst, vjfirst, vlfirst, vni, vnj, vnl
                # print shape(vdata_array), shape(vdata)
                
                
                
                # #c: nested model, we use all the data 0:vni, 0:vnj etc
                vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[0:vni, 0:vnj, (vlfirst-1):(vlfirst+vnl-1)]
            else:
                if (vni*vnj>1):
                    # print vifirst, vjfirst, vlfirst, vni, vnj, vnl
                    # #c: both vni and vnj have been set 
                    
                    vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[(vifirst-1):(vifirst+vni-1), \
                                                               (vjfirst-1):(vjfirst+vnj-1), \
                                                               (vlfirst-1):(vlfirst+vnl-1)]
                else:
                    # #c: vni or vnj not in use
                    vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[0:vni, 0:vnj, (vlfirst-1):(vlfirst+vnl-1)]
            
            # S6: create bpch2_data class
            
            pbdata=bpch2_data(vtracer_id,vcategory.strip(), vunit.strip(),\
                                  data_grid,vdata, \
                                  tau0=vtau0,tau1=vtau1, modelname=vmodelname.strip(),\
                                  reserved=vreserved, offset=offset, name=traname.strip(), \
                                  fullname=trafullname.strip(), molew=tramolew, \
                                  scale=trascale, traunit=traunit.strip(), id=real_id)
                            
            
            # S7: added to data list 
            
            self.data.append(pbdata)
            
            self.ntracers=self.ntracers+1
            
            # print 'traname', 'fullname', 'ntracer:', traname, trafullname, self.ntracers
            # S8: read in next record
            
            vtracer_id,vhalfpolar,vcentre180,\
                vni,vnj,vnl,vifirst,vjfirst,vlfirst,vlonres,vlatres,\
                vtau0,vtau1,vmodelname,vcategory,\
                vunit,vreserved,vdata_array,stat = bpch2_rw_py.read_bpch2_record(funit)
                                                
        
            
                    
        if (stat==-1 or stat==29):
            self.stat=0
            self.read_mode=1
        else: 
            self.stat=stat
            # problem occur
            self.read_mode=-1
        
        # loop while end 
         
        stat=bpch2_rw_py.close_bpch2_file(funit)
                
        return self.stat
    
    def read_mode_2(self,flnm, categorys, tracers, taus):
        """
        Read only selected tracer data from bpch2 file
        
        Inputs:
        -------------------------------------------------------------------
        1. flnm:<str>: name of bpch2 file 
        2. categorys:<list, t:str>: list of category 
        3. tracers:<list, t:int>: list of tracer ID 
        4. taus:<list, t:float>: list of time 
      
        
        Returns:
        -------------------------------------------
        1. stat:<int>: IO status 
        """
        
        # S1: open file  
        funit=199
        fti,title,stat = bpch2_rw_py.open_bpch2_for_read(funit,self.flnm)
        
        self.title=title.strip()
        
        
        if (stat<>0):
            print 'error in read :',  stat
            self.stat=self.stat
            return stat

        # S2: read in selected tracer, and stored them into data list 
        
        for iname in range(len(categorys)):
            cname=categorys[iname]
            tid=tracers[iname]
            tau0=taus[iname]
            # S3: read in record
            # print 'cnmae, tid, tau0', cname, tid, tau0
            
            vtracer_id,vhalfpolar,vcentre180,\
                vni,vnj,vnl,vifirst,vjfirst,\
                vlfirst,vlonres,vlatres,\
                vtau0,vtau1,vmodelname,vcategory,\
                vunit,vreserved,vdata_array,stat = \
                bpch2_rw_py.sel_bpch2_record(funit, cname, tid, tau0)

            
            if (stat==0):
                # S4: create grid 
                data_grid=bpgrid(vni, vnj, vnl,\
                                     vifirst, vjfirst, vlfirst,\
                                     vhalfpolar, vcentre180, \
                                     vlonres, vlatres)
              # S5: decode tracer info and diag info
                offset=self.diaginfo.get_offset(vcategory)
                if (offset<0):
                    offset=0
                                        
                real_id=vtracer_id+offset
                traname, trafullname, tramolew, \
                    trascale, traunit=self.tracerinfo.get_tracer_info(real_id)
                
                vdata=zeros([vni, vnj, vnl], float)
                # S6: fill data 
                vdata[0:vni, 0:vnj, 0:vnl]=vdata_array[(vifirst-1):(vifirst+vni-1), \
                                                           (vjfirst-1):(vjfirst+vnj-1), \
                                                           (vlfirst-1):(vlfirst+vnl-1)]
                
                # S7: creat bpch2_data 
                pbdata=bpch2_data(vtracer_id,vcategory.strip(), vunit.strip(),\
                                      data_grid,vdata, \
                                      tau0=vtau0,tau1=vtau1, modelname=vmodelname.strip(),\
                                      reserved=vreserved, offset=offset, name=traname.strip(), \
                                      fullname=trafullname.strip(), molew=tramolew, \
                                      scale=trascale, traunit=traunit.strip(), id=real_id)
                # S8: add the data to list 
                self.data.append(pbdata)
                self.ntracers=self.ntracers+1
                                    
            else:
                break
            
            
        # loop iname end
    
        stat=bpch2_rw_py.close_bpch2_file(funit)
        return stat
    

    

if (__name__=='__main__'):
    flnm='ts_satellite.ST001.EN0081-EN0146.20090225.bpch'
    path='~/local_disk/otool_data/enkf_output/2009/'
    full_flnm=path+flnm
    bpch2=bpch2_file_rw(full_flnm, 'r', do_read=1)
    print'ntrac',  bpch2.ntracers
    
    bpch2.print_datainfo()
    data=bpch2.read_bpch2_data(30, True)
    print shape(data.data)
    print data.grid.centre180
    
    data_list, founded=bpch2.get_data(None, None, 160056.0, 'PSURF')
    print 'founded', founded, len(data_list)
    for bpdata in data_list:
        traname, tau0, tau1=bpdata.get_attr(['name', 'tau0', 'tau1'])
        print bpdata.category, tau0, tau1, traname.strip()
        print bpdata.data[3,3,0]
    
    test_write=True
    if (test_write):
        bpch2_w=bpch2_file_rw("restart_new.bpch", 'w')
        ifunit=95
        bpch2_w.open_bpch2_w(ifunit,"test")
        bpch2_w.add_bpch2_data(data)
        print 'I am here'
        bpch2_w.add_bpch2_data(data)
        data=bpch2.read_bpch2_data(20, True)
        print 'the second read'
        bpch2_w.add_bpch2_data(data)
        print 'try to write'
        bpch2_w.write_bpch2_data(0)
        bpch2_w.write_bpch2_data(1)
        bpch2_w.close_bpch2_file()

    
  

  

            
            
        
