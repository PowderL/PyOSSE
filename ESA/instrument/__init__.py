"""
directory for codes to model instrument, which include

1. orbit_m.py  # read orbit 
2. ak_m.py    # read averaging kernel tables, and generate scene-dependent average kernels 

3. err_m.py   # read averaging kernel tables, and generate scene-dependent average kernels 
 
# scene  screening 
4. cloud_m.py  # use cloud map to screen contaminated scenes. 
5. aerosol_m.py  # get aerosol aod from climatology PDFs. 

"""
import cloud_m
import cloud_file_m
import aod_m
import aod_file_m
import ak_m
import ak_file_m
import orbit_m
import orbit_file_m

