import interpolate_f as isch
from numpy import *
a=arange(300)
b=-2
pr=isch.search_upper_boundary(a, b)
print pr

# print a[pr], b
a[11]=-999.0
b=30
pl, pr=isch.getpos_mask_outsider(a, b, mask_val=-999.0)
print pr
print pl


b=[30.5, 20.5, -999.0, 11.5, 310, -2]
pl, pr=isch.getpos_mask_outsider(a, b, mask_val=-999.0)
print pr
print pl


pl, pr,wgt=isch.getwgt(a, b, mask_val=-999.0)
print pr
print pl
print wgt 

# print a[pr]
# print b

# print a[pl]

zx=isch.get_interpl(a, pl,pr,wgt, mask_val=-999.0)
print zx

zx=isch.get_interpl_1d(a, pl,pr,wgt, mask_val=-999.0)
print shape(zx)

print zx

pl, pr,wgt=isch.getwgt_mask_outsider(a, b, mask_val=-999.0)
print pr
print pl
print wgt 

zx=isch.get_interpl_1d(a, pl,pr,wgt, mask_val=-999.0)
print shape(zx)
print zx



