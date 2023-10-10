import ESA.util.otool_ncfile_io as ncfio
import ESA.util.geo_constant as gc

from numpy import *
kg_s_to_GtC_Y=gc.mc*3600.0*365*24.0/(1.0e12*gc.mco2)

flnm='gc_co2.20090101.nc'
varnames=['BB', 'OC', 'FF', 'BAL', 'BF']
BB, OC, FF, BAL, BF=ncfio.ncf_read(flnm, varnames)
print 'BB:',  kg_s_to_GtC_Y*sum(BB)
print 'OC:',  kg_s_to_GtC_Y*sum(OC)
print 'FF:',  kg_s_to_GtC_Y*sum(FF)
print 'BAL:',  kg_s_to_GtC_Y*sum(BAL)
print 'BF:',  kg_s_to_GtC_Y*sum(BF)



