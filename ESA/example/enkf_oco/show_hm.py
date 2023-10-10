import ESA.util.otool_ncfile_io as ncfio
import ESA.util.geo_constant as gc
from pylab import *
from numpy import *
kg_s_to_GtC_Y=gc.mc*3600.0*365*24.0/(1.0e12*gc.mco2)

flnm='./oco_inv/hm.20090125.nc'
varnames=['lon', 'lat', 'hm']
lon, lat, hm=ncfio.ncf_read(flnm, varnames)
print lon[1000], lat[1000]
print 1.0e6*hm[1000, 80:120]

pcolor(hm)
show()



