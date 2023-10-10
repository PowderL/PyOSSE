import numpy as npy
import sigma_pres_mod as pm
nlvl, ap, bp=pm.get_gc_ap_bp(5, 1)
print nlvl
ap=ap[0:nlvl]
bp=bp[0:nlvl]
ps=npy.zeros(2, float)
ps[0]=1000.0
ps[1]=850.0
pres, pres_edge=pm.get_model_pres(ps, ap, bp)
print npy.shape(pres), npy.shape(pres_edge)
pres=pres[:,0,:]
pres_edge=pres_edge[:,0,:]

print 'pres at level 0', pres[:,0]
print 'pres edge at level 0', pres_edge[:,0]

print 'pres at level -1', pres[:,-1]
print 'pres edge at level -1', pres_edge[:,-1]



ps=1030.0
print npy.shape(ps)

pres, pres_edge=pm.get_model_pres(ps, ap, bp)


print npy.shape(pres), npy.shape(pres_edge)
pres=pres[0,0,:]
pres_edge=pres_edge[0:,0,:]
