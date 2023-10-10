"""
Functions to divide the regions into smaller areas 

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2013.01.12
    
    

    Functions:
    ==============================================================
    1. divide_region_to_subregion: divide invidual region to sub-regions
    2. divide_region_to_box: divide invidual region to boxes of given sizes.
    3. divide_map: divide selected regions into required sub-region or boxes. 
    
    Classes:
    ==============================================================
    1. reg_op_cl:class for defining actions  to split regions. 
    
    """

import ESA.util.otool_obj as oob
import ESA.util.message_m as msm
import ESA.util.geo_constant as gc
import ESA.util.otool_var_io as ovar
import ESA.util.otool_ncfile_io as ncfio
import ESA.atmosphere.ctm_grid_2d as grd_2d
import ESA.util.flux_regrid as fgrd
import ESA.util.gen_plots as gpl
import  numpy as npy
import pylab as plb



class reg_op_cl:
    """
    class for paramter and action to split one region 
    
    Members:
    ------------------------------------------------------
    1. reg_id:<int>: region ID for parent region
    2. nx: <int>: 
    ---If split_op=divide_region_to_subregion, nx is sub-region number 
    ---along longitude.
    ---If split_op=divide_region_to_box, nx is 
    ---the number of cells along longitude for each box. 
    
    3. ny:<int>:
    ---If split_op=divide_region_to_subregion, ny is sub-region number 
    ---along latitude.
    ---If split_op=divide_region_to_box, ny is 
    ---the number of cells along latitude for each box. 
    
    4. split_op:<func>: actions to split the region
    

    Functions:
    -----------------------------------------
    1. __init__: initialization
    2. do_split:  split regions as required
    

    """
    def __init__(self, reg_id, nx, ny, \
                     split_op):
        
        """    
        Initialization 
        
        Inputs:
        ---------------------------------------
        1. reg_id:<int>: region ID for parent region
        2. nx: <int>: 
        ---If split_op=divide_region_to_subregion, nx is sub-region number 
        ---along longitude.
        ---If split_op=divide_region_to_box, nx is 
        ---the number of cells along longitude for each box. 
        
        3. ny:<int>:
        ---If split_op=divide_region_to_subregion, ny is sub-region number 
        ---along latitude.
        ---If split_op=divide_region_to_box, ny is 
        ---the number of cells along latitude for each box. 
        4. split_op:<func>: actions to split the region
        
        """
        
        self.reg_id=reg_id
        self.nx=nx
        self.ny=ny
        self.split_op=split_op
    
    def do_split(self, inmap, lon, lat,  **keywords):
        """
        split region into sub-regions or boxes

        Inputs:
        ----------------------------------------------------
        1. inmap:<array, (nlon, lat)>: region map 
        2. lon:<array, (nlon)>: longitude 
        3. lat:<array, (nlat)>: latitude 
        4. keywords:<dict>: dictionary for extra inputs
        ---reserved keywords
        ---nx:<int>: replacement for self.nx
        ---ny:<int>: replacement for self.ny
        
        
        Returns:
        -------------------------------------------------
        1. combine_map:<array, (nlon, nlat)>: map for sub-regions 
        2. map_lst:<list, t:array>: layer map for sub-regions

        
        Notes:
        ------------------------------------------------------
        1. combine_map masks sub-regions from 1 to its total
        
        """

        if ('nx' in keywords):
            self.nx=nx
        
        if ('ny' in keywords):
            self.ny=ny
            
        combine_map, map_lst=self.split_op(inmap, \
                                               lon, lat, \
                                               self.nx, self.ny, self.reg_id)
        return combine_map, map_lst
    
    
        
def divide_region_to_subregion(type_map_in, \
                                   lon, lat, nsb_lon, \
                                   nsb_lat, reg_id):
    
    
    """
    
    divide a region into sub-regions defined by (nsb_lon, nsb_lat)
    
    Inputs:
    --------------------------------------
    1. type_map_in:<array, (nlon, nlat)>: the map defining the regions 
    2. lon:<array, (nlon,)>: longitude 
    3. lat:<array, (nlat,)>: latitude
    4. nsb_lon:<integer>: number of sub regions along lon 
    5. nsb_lat:<integer>: number of sub regions along lat 
    6. reg_id:<integer>: region id
    
    Returns:
    --------------------------------------------
    1. combine_map:<array, (nlon, nlat)>: new 
    ---regions masked with sub-region number 
    2. new_map_lst:<list, t:array>: the layer map 
    ---with each layer for each sub-regions
    

    Notes:
    -----------------------------------
    1.Algorithm for dividing regions 
    ---> separate regions into (nsub_lat) latitude bands according to 'area'
    ---> divede each latitude band into (nsub_lon) longitude blocks according to 'area'

    2. combine_map masks sub-regions from 1 to its total
    
    """
    
    
    fpi=npy.pi/180.0

    # S1 get parent region 
    
    type_map=npy.where(type_map_in==reg_id, 1, 0)
    
    
    nx, ny=npy.shape(type_map)
    
    nlayer=nsb_lat*nsb_lon
    
    new_map=npy.zeros([nx, ny, nlayer+8]) # +8 to accomodate points not allocated into any sub regions
 
    cap_map=npy.zeros([nx, ny, nlayer])  # for reconstructing map to asorb left point
    
    if (nsb_lat*nsb_lon==1):
        # ##c if no sparation needed 
        new_map_lst=[]
        new_map_lst.append(type_map)
        return  type_map, new_map_lst
    
    
    ## T1 get grid boxes (cells) in the region   
    
    idx=npy.where(type_map==1)
    xcell_id=idx[0]
    ycell_id=idx[1]
    
    xcell=lon[xcell_id] # lons of  region cells
    ycell=lat[ycell_id] # lats of  region cells
    
    ## T2 get rectangle covering the whole region   
    
    
    xmin, xmax=npy.min(xcell), npy.max(xcell) 
    ymin, ymax=npy.min(ycell), npy.max(ycell)
    
    # ##c xrange 
    
    xid_rect=npy.where((lon>=xmin) & (lon<=xmax))
    
    xid_rect=npy.squeeze(xid_rect)  ## ids for lon range 
    
    if ((xmax-xmin)>300): 
        # ##c  if rectangle is too large along longitude, 
        # ##c  we need to shift points with lon<0 by 360 degree
        
        nx=npy.size(lon)
        # ##c find the regional cells needed to be shifted
        # ##c shift cell longitude by 360 degree (or by nx)
        
        neg_id=npy.where(xcell<0)
        xcell_id[neg_id]=xcell_id[neg_id]+nx
        
        
        # ##c shift longitude range by 360 degree (or by nx)
        
        x_sq=lon[xid_rect]
        neg_id=npy.where(x_sq<0)
        xid_rect[neg_id]=xid_rect[neg_id]+nx
        
        # ##c sort xid_rect in an assending order
        
        xid_rect=npy.sort(xid_rect) 

    
    # ##c get latitude range  

    yid_rect=npy.where((lat>=ymin) & (lat<=ymax))
    yid_rect=npy.squeeze(yid_rect)
    

    
    
    
    
    
    # S3 divide latitude bands

    ## T1 collect area and cell idx along y (latitude )
    
    sel_cell_id_list=list()  # list for cell id with with same latitude
    areas=list() # list for areas with same latitude 
    
    
    for iy in yid_rect:
        cur_idx=npy.where(ycell_id==iy)
        # id in xcell_id, ycell_id
        
        cur_idx=npy.squeeze(cur_idx)
        if (npy.size(cur_idx)>=1):
            cix, ciy=xcell_id[cur_idx], ycell_id[cur_idx]
            ###C area_y=areas of cells with latitude==lat(iy)
            
            area_y=npy.sum(abs(npy.cos(fpi*lat[ciy]))) 
            
            areas.append(area_y)
            sel_cell_id_list.append(npy.array(cur_idx))
        
        else:
            pass
    
    ny_usd=len(areas) # number of latitude points falling into the regions
    
    areas=npy.array(areas)
    
    

    area_div_y=npy.sum(areas)/nsb_lat  
    
    # S4  divide into sub-regions (see Note 1)
    
    # ##C Variables
    # ##C ---1:  area_rect_x:<array>:  areas of the sub-latitude bands for each x grid 
    # ##  ---2:  sub_cur_idx:<list, t:integer>: each of its elements is the id for 
    # ## ------cells in the sub-latitude bands. 
    
    
    area_rect_x=npy.zeros(npy.size(xid_rect), float)
    sub_cur_idx=list()
    
    
    
    ## T1  loop over every used y value. 
    
    # ##C region number 
    
    region_count=0
    area_sum=0.0
    
    maxlvl=-1 # levels (i.e, sub region) already generated . 
    
    
    for iy in range(ny_usd):
        
        area_sum=area_sum+areas[iy]  # aggregated areas of latitude band 
        
        ## T2 allocate area of the cells of the given latitude
        ## --to area_rect_x according to their x values
        
        
        cur_cell_idx=sel_cell_id_list[iy]


        if (npy.size(cur_cell_idx)>=1):
            
            tmp_xorg_id=xcell_id[cur_cell_idx]
            # ##c find the poistion of cells in rectangle xid_rect
            
            xpos_in_rect=npy.searchsorted(xid_rect, tmp_xorg_id) 
            xpos_in_rect=npy.squeeze(xpos_in_rect)
            
            if (npy.size(cur_cell_idx)==1):
                ciy=cur_cell_idx
            else:
                ciy=cur_cell_idx[0]
                
            # ##c calculated cell 'area'
                
            cell_area=abs(npy.cos(fpi*lat[ycell_id[ciy]]))
            
            # ##c 1. accumulate cell areas into area_rect_x
            # ##c 2. insert cell idx into sub_cur_idx 
            
            for ii in range(npy.size(cur_cell_idx)):
                if (npy.size(xpos_in_rect)==1):
                    sqix=xpos_in_rect
                    cix=cur_cell_idx
                else:
                    sqix=xpos_in_rect[ii]
                    cix=cur_cell_idx[ii]
                
                area_rect_x[sqix]=area_rect_x[sqix]+cell_area
                sub_cur_idx.append(cix) 
                
            ## T3: divide latitude bands into sub-regions 
            if (area_sum>=area_div_y or iy==ny_usd-1):
                
                # ##  if the sum of areas of the latitude latitudes 
                # ##c  reach the threshold
 
                sub_cur_idx=npy.array(sub_cur_idx)
                
                # ##c divide cells of sub_cur_idx along the x (longitude)
                
                area_div_x=npy.sum(area_rect_x)/nsb_lon
                area_rect_x=npy.add.accumulate(area_rect_x)
                
                x_div_ar=npy.arange(0, nsb_lon,1)
                x_div_ar=area_div_x*x_div_ar
                
                # ##c xpos_in_lvl is the layer number for each x values
                             
                
                
                xpos_in_lvl=npy.searchsorted(x_div_ar, area_rect_x)
                xpos_in_lvl=xpos_in_lvl-1

                # ##C select id in map (cix, ciy)
                
                cix_lst=xcell_id[sub_cur_idx]
                ciy_lst=ycell_id[sub_cur_idx]

                ## T4:  project cell into x rectangle range
                
                
                tmp_x_rect_pos=npy.searchsorted(xid_rect, cix_lst)
                ## T5: find sub region (level number) for each cell according to clvl 
                ## (i.e, x values)
                
                for sel_ii in range(npy.size(sub_cur_idx)):
                    if (npy.size(cix_lst)==1):
                        cix=cix_lst
                        ciy=ciy_lst
                        clvl=tmp_x_rect_pos
                    else:
                        cix=cix_lst[sel_ii]
                        ciy=ciy_lst[sel_ii]
                        clvl=tmp_x_rect_pos[sel_ii]
                    
                    ilvl=region_count+xpos_in_lvl[clvl]
                    if (ilvl>maxlvl):
                        maxlvl=ilvl
                    
                    # ##c project back to normal index if necessary

                    if (cix>=nx):
                        cix=cix-nx

                    
                    
                    new_map[cix, ciy, ilvl]=type_map[cix, ciy]

                region_count=maxlvl+1
                
                ## T6: re-initialize for next latitude band
                
                area_sum=0.0
                
                area_rect_x=npy.zeros(npy.size(xid_rect), float)
                sub_cur_idx=list()
                
        # print iy, area_sum, ny_sq
    
        # add the last 2 maps into the nlayer-1 ones
        
    ## T7 asorb mis-allocated points to the top layer 
                
    add_map=npy.sum(new_map[:,:,nlayer:], axis=2)
    cap_map[:,:,0:nlayer]=new_map[:,:,0:nlayer]
    cap_map[:,:,nlayer-1]=cap_map[:,:,nlayer-1]+add_map[:,:]
    
    # S5 construct combined map, and remove layers without any cells. 
    
    combine_map=npy.zeros(npy.shape(cap_map[:,:,0]))
    
    new_map_lst=[]
    
    icount=0
    for ilvl in range(nlayer):
        tx=cap_map[:,:,ilvl]
        
        if (npy.max(tx.flat)>0):
            new_map_lst.append(tx)
            icount=icount+1
            tx=npy.where(tx==1, icount,0)
            combine_map=combine_map+tx
    
    
    
                        
    return combine_map, new_map_lst




def divide_region_to_box(type_map_in, lon, lat, nx_box, ny_box, reg_id):
    
    """
    
    divide a region into sub-regions defined by (nsb_lon, nsb_lat)
    
    Inputs:
    --------------------------------------
    1. type_map_in:<array, (nlon, nlat)>: ratio  map defining one region. 
    
    2. lon:<array, (nlon,)>: longitude 
    3. lat:<array, (nlat,)>: latitude
    4. nx_box:<integer>: number of longitude grids for one box.   
    5. ny_box:<integer>: number of longitude grids for one box.
    6. reg_id:<integer>: region id
    
    Returns:
    --------------------------------------------
    1. combine_map:<array, (nlon, nlat)>: new regions masked with sub-region number 
    2. new_map_lst:<list, t:array>: the layer map with each layer for each sub-regions
    
    
    
    """
    
    
    
    # S1 get parent region 
    
    type_map=npy.where(type_map_in==reg_id, 1, 0)
    
    nx, ny=npy.shape(type_map)
    new_map=npy.zeros([nx, ny])
    
    
    idx=npy.where(type_map==1)
    xcell_id=idx[0]
    ycell_id=idx[1]
    
    
    xcell=lon[xcell_id] # lons of  region cells
    ycell=lat[ycell_id] # lats of  region cells
    
    ## T2 get rectangle covering the whole region   
    
    # ##c longitude and latitude range
    
    xmin, xmax=npy.min(xcell), npy.max(xcell) 
    ymin, ymax=npy.min(ycell), npy.max(ycell)
    
    # ##c longitude index  range
        
    xid_rect=npy.where((lon>=xmin) & (lon<=xmax))
    xid_rect=npy.squeeze(xid_rect)  ## ids for lon range 
    
    len_xid=npy.size(xid_rect)

    # ##c latitude index range 
    
    
    yid_rect=npy.where((lat>=ymin) & (lat<=ymax))
    yid_rect=npy.squeeze(yid_rect)


    yid_rect=npy.squeeze(yid_rect)

    len_yid=npy.size(yid_rect)
    
    
    # check each cell within the rectangle, and assign valid cell to each boxes
    
    ireg=1  # region number starting with 1
    
    
    ycnt=0
    max_reg=0
    for idy in yid_rect[::-1]:  # over y starting from top (North)
        

        new_x_map=type_map[:, idy]  # map along latitude 
        new_x_id=npy.where(new_x_map>0)
        new_x_id=npy.squeeze(new_x_id)
        
        xcnt=0
        ireg_x=0
        
        
        for idx in new_x_id:  # over x
            new_map[idx, idy]=ireg+ireg_x
            type_map[idx, idy]=-1
            xcnt=xcnt+1
            
            if ((ireg+ireg_x)>max_reg):
                max_reg=ireg+ireg_x
            
            if (xcnt>=nx_box):
                ireg_x=ireg_x+1
                xcnt=0
        
        # ##c idx
        ycnt=ycnt+1
        
        if (ycnt>=ny_box):
            ireg=max_reg+1
            ycnt=0
        
    
    # #for idy
    
                        
    new_map_lst=list()
    nlayer=int(npy.max(new_map.flat))
    

    
    for ilvl in range(nlayer):
        tx=npy.where(new_map==ilvl+1, 1, 0)
        new_map_lst.append(tx)
        
    
    combine_map=new_map
    return combine_map, new_map_lst

    





def divide_map(type_map, lon, lat, reg_op_lst):
    
    """
    divide regional map into more sub-regions. 

    Inputs:
    -------------------------------------------------------
    1. type_map_in:<array, (nlon, nlat)>: the map defining the regions 
    2. lon:<array, (nlon,)>: longitude 
    3. lat:<array, (nlat,)>: latitude
    4. reg_op_lst:<list, t:reg_op_cl>: list for class reg_op_cl on 
    ---how to divide each regions
    
    
    Returns:
    ---------------------------------------------------------
    1. new_map:<array, (nlon, nlat)>: map for more regions
    2. pid_lst:<list, t:int>: list  for parent ids of the sub regions
    
    
    """
    icount=0
    all_region=[]
    new_map=None
    pid_lst=[]
    
    nlvl=0

    
    for reg_op in reg_op_lst:
        
        combine_map, new_map_lst=reg_op.do_split(type_map, lon, lat)
        
            
        if (new_map==None):
            new_map=combine_map
            # print npy.shape(new_map)
            total_reg=int(npy.max(new_map.flat))
        
        else:
            # shift the combine map by nlvl 
            
            combine_map=npy.where(combine_map==0, 0, combine_map+nlvl)
            
            # print npy.shape(combine_map), npy.shape(new_map)
            new_map=new_map+combine_map
        
        nlvl=int(npy.max(new_map.flat))
        
        nlayer=len(new_map_lst)
        
        reg_pid=nlayer*[reg_op.reg_id]
        pid_lst=pid_lst+reg_pid
        

    
    
    
    
    return new_map, pid_lst




 

    
# <<<<  TEST  >>>>

if (__name__=='__main__'):
    
    varnames=['longitude', 'latitude', 'layer', 'map', 'flux', 'area', 'pid']
    ncflnm='t3_reg_flux_05x05.nc'
    lon, lat, layer, reg_map, flux, area= ncfio.ncf_read(ncflnm, varnames)
    
    reg_op_lst=[]
    nsnow=1
    nland=11
    nocean=11
    reg_id=1
    
    reg_op=reg_op_cl(reg_id,1,1, divide_region_to_subregion)
    reg_op_lst.append(reg_op)
    
    
    for i in range(nland):
        reg_id=reg_id+1
        reg_op=reg_op_cl(reg_id, 3, 3, divide_region_to_subregion)
        if (reg_op.reg_id==4):
            print 'Divide region 4 into boxes ' 
            reg_op=reg_op_cl(reg_id, 20, 20, divide_region_to_box)
        
        reg_op_lst.append(reg_op)
    
    for i in range(nocean):
        reg_id=reg_id+1
        reg_op=reg_op_cl(reg_id, 2, 2, divide_region_to_subregion)
        reg_op_lst.append(reg_op)


    
    print npy.shape(reg_map)
    
    new_map, pid_lst= divide_map(reg_map, lon, lat, reg_op_lst)
    print npy.max(new_map.flat)
    
    
    print pid_lst
    
    gpl.plot_map(new_map, lon, lat, use_pcolor=1)
    
    plb.show()
    
    
    

    

                              
