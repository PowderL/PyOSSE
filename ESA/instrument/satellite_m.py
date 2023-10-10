""" 
   Virtual class for defining satellite 
   These functions to be overrided by the user
   
   Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.02.15
    
    Classes:
    ===============================================
    1. satellite_cl: class for satellite information
    

"""


def_viewsize=2  #  km^2  
def_viewangle=0.0 # 
def_sat_alt=800.0*1000.0 #  800 km 
def_sat_speed=7.5*1000.0 #  7.5 km/s 


class satellite_cl:
    def __init__(self, name, sid):
        self.name=name 
        self.id=sid
        
    def get_groundspeed(self, sat_lon, sat_lat, **keywords):
        """
        get satellite ground speed 
        """
        return def_sat_speed
    
    def get_alt(self, sat_lon, sat_lat, tau, **keywords):
        """ get satellite altitude 
        """
        return def_sat_alt


    def get_viewsize(self, viewtype, viewmode, viewangle, **keywords):
        """
        get view size at the ground 
        """
        
        return def_viewsize
  

    def get_viewangle(self, olon, olat, sat_lon, sal_lat, sat_alt, **keywords):
        """
        get viewangle 
        """
        return def_viewangle 

    
    def get_sat_loc(self, tau, **keywords):
        
        """
        get  the location of the satellite 
        
        """
        
        return None, None
    

    def get_orb(self, tau, **keywords):
        
        """
        get  the location of the satellite 
        
        """
        lons=[]
        lats=[]
        
        return lons, lats

    
