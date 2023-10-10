""" generic class for grid used by numerical models.  
   
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.08.23
    this module defines grid class for high-dimension fields


Classes:
===================================================
1. gp_grid_cl: grid 


Functions:
=====================================================
1.  coordinate_sample_index: Generate (0,1) (i.e, low-high) 
----list for the corners around a  given point in a high-dimension space. 
    
2. build_grid_from_matrix: convert a matrix into a gridded data set. 


"""
import numpy as npy
import gp_axis_m as axis_m
import otool_obj as oob
import message_m as msm
import copy as obj_copy

def coordinate_sample_index(n):

    """ Generate (0,1) (i.e, low-high) list for the corners around a 
    given point in a high-dimension space. These indexes can be 
    used to sample the (left-right) boundary pairs of the point along each directions. 
    
    Inputs:
    ----------------------------------------------
    1. n:<int>: dimension number of the space (==the number of left-right boundary pairs)
    
    Returns:
    ----------------------------------------------
    1. cor_lst:<list>: list of the (0,1)  for the surrounding corners of 
    ---a point in a N-dimensional space. 
    For example,  in a 2D space, cor_lst=[[0,0], [0,1], [1,0], [1,1]] for four corners. 
    
    """
    # for dimension size n, we have 2^n corners encircling one given point. 
    
    
    cor_id_lst=list()
    zx=n*[0]
    
    cor_id_lst.append(npy.array(zx))
    
    ibx=range(n)
    
    for i in range(1, 2**n):
        for ib in ibx:
            zx[ib]=zx[ib]+1
            if (zx[ib]<2):
                break
            else:
                zx[ib]=0
        
        cor_id_lst.append(npy.array(zx))
    
    
    
    return cor_id_lst


class gp_grid_cl: 
    
    """ grid class
    
    
    Members:
    ------------------------------------------------------
    1. ndim:<int>: number of the dimensions. 
    2. dims:<list>: shape of the grid 
    3. __axis_dict:<dictionary>: dictionary for id of gp_axis
    4. axis_lst:<list>: list of the axises
    5. attr_dict:<dictionary>: dictionary for the attributes
    
    
    Functions:
    ------------------------------------------------------
    1. __init__: initialize the class
    2. __getitem__: overide index function to return axis. 
    3. set_attr:  set attribute
    4. get_attr:  get attribute
    5. allocate_point: find grid corners surrounding a given point.  
    6. allocate_point_subspace: find grid corners surrounding a point given in a subplace
    
    7. update_axis: replace or add axis to the grid 

    8. get_axis_ids: get IDs for given axis names
    
    9. get_value_at_point: interpolate data field to a point defined at the grid (N-d) space. 
    
    10. get_profile: interploate profiles of data field to a point defined at (N-1)-d space 
    
    
    """
    
    
    def __init__(self, gp_axis_lst, **keywords):
        
        """initialize grid 
        Inputs:
        1. gp_axis_lst:<list>:list of gp_axis for definition of the grid.
        2. keywords:<dict>: attributes

        Notes:
        ----------------------------------------------------
        each axis in gp_axis_lst  is expected have different names 
        
        """
        
        # shape of the grid 
        self.ndim=len(gp_axis_lst)
        self.dims=[]
        
        for gax in gp_axis_lst:
            self.dims.append(gax.ax_np)
        
        # shape we expected
        
        # private member
            
        self.__axis_dict={}
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_grid})
        
        self.__axis_lst=list(gp_axis_lst)
        
        ax_id=0
        # construct axis_dict for name --> axis_id
        
        for gax in self.__axis_lst:
            self.__axis_dict.update({gax.ax_name:ax_id})
            ax_id=ax_id+1
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
            
    
    def __getitem__(self, ax_id):  
        """ get axis  
        
        Inputs:
        --------------------------------
        1. ax_id:<integer/str>: Id for the axis to be found
        
        Returns:
        ----------------------------------
        1, gax:<gp_axis_cl>: axis if the ax_id is valid
        
        """
        
        if (type(ax_id)==type('name')):
            if (ax_id in self.__axis_dict):
                ls_id=self.__axis_dict[ax_id]
                gax=self.__axis_lst[ls_id]
               
                
                return gax
    
            else:
                
                msg='axis '+ax_id+' not exist'
                msm.show_err_msg(msg)
                return None
            
        elif (ax_id<self.ndim):
            gax=self.__axis_lst[ax_id]
            return gax
        else:
            msg='Axis Id number is larger than the grid dimensions'
            msm.show_err_msg(msg)
            return None

    
    def get_axis(self, ax_id):

        """wrapper for __getitem__
        
        Inputs:
        --------------------------------
        1. ax_id:<integer/str>: Id for the axis to be found
        
        Returns:
        ----------------------------------
        1, gax:<gp_axis_cl>: axis if the ax_id is valid
        
        
        """

        if (type(ax_id)==type('name')):
            if (ax_id in self.__axis_dict):
                ls_id=self.__axis_dict[ax_id]
                gax=self.__axis_lst[ls_id]
                return gax
        
            else:
                
                msg='axis '+ax_idx+' not exist'
                msm.show_err_msg(msg)
                return None
            
        elif (ax_id<self.ndim):
            gax=self.__axis_lst[ax_id]
            return gax
        else:
            msg='Axis Id number is larger than the grid dimensions'
            msm.show_err_msg(msg)
            return None
    
    
    def update_axis(self, gp_axis):
    
        """ add or replace axis 
        
        Inputs:
        --------------------------------------------------
        1. gp_axis:<gp_axis_cl>: class to be added or replaced
                
        Output:
        -------------------------------------------------
        2. gax:<gp_axis>: axis if the Id is correct
        
        """
        
        ax_name=gp_axis.ax_name
        
        if (ax_name in self.__axis_dict):
            # replace existing one 
            ax_id=self.__axis_dict[ax_name]
            self.__axis_lst[ax_id]=gp_axis
            # update dims
            self.dims[ax_id]=gp_axis.ax_np
            
            
        else: # append it
            
            self.ndim=self.ndim+1
            self.__axis_lst=self.__axis_lst+[gp_axis]
            self.__axis_dict.update({ax_name:self.ndim-1})
            
            if (oob.get_ot_type(self.dims)==oob.ot_list):
                self.dims.append(gp_axis.ax_np)
            elif (oob.get_ot_type(self.dims)==oob.ot_array):
                self.dims=concatenate((self.dims, [gp_axis.ax_np]))
            
            else:
                self.dims=self.dims+[gp_axis.ax_np]
            
                                     
    
    def get_attr(self, attr_name):
        
        """ find the attribute
        
        Inputs:
        ----------------------------------
        1. attr_name:<str>: name of the attribute
        
        Returns:
        ------------------------------------------
        1. attr: <obj>: if the attribute exists. 
        """
        
        if (attr_name in self.attr_dict):
            return self.attr_dict[attr_name]
        else:
            msg=attr_name
            msm.show_err_msg(msg, msm.msm_no_attr)
            return None
    
    def set_attr(self, attr_name, attr_value):
        
        """ set attribute
        Inputs:
        ----------------------------
        1. attr_name:<str>: name of the attribute
        2. attr_value:<obj>: the value of the attribute 
        """
        
        self.attr_dict.update({attr_name:attr_value})
        
    
    def allocate_point(self, p):
       
        """ find grid points around a given point, and  
        calculate their weights according to the distances from this point.  
        
        Inputs:
        ----------------------------------------------------------------
        1. p:<array, (nax,)>: the coordiantes of the given point.  
        
        Returns:
        -------------------------------------------------------------------
        1. grd_pt_lst: <list, t:tuple>: list of the index of the surrounding grid points.
        ---If the dimension size nax is n, the length of grd_pt_lst is 2^n, 
        with each element being an array of size n. 
                
        2. grd_wgt_lst: <list, t:float>: Weights of the the surrounding points,         
        ---which are calculated from their distances from the point.
        
        """
        
        
        p1_lst=list() # list for left-boundary
        p2_lst=list() # list for right-boundary
        w1_lst=list() # list for weighting 

        
        #S1 get position and weights for boundary points at each axis
        
        iax=0
        for val in p:
            gax=self.__axis_lst[iax]
            p1, p2, wgt=gax.getwgt(val)
            p1_lst.append(p1)
            p2_lst.append(p2)
            w1_lst.append(wgt)
            iax=iax+1
            
        
        # S2 coordiante indexes for the corners (boundaries)
        # cor_lst is a list of left-rigt (0 or 1) coordiantes at ndim-shace, 
        # its length is 2^ndim 
             

        
        cor_lst= coordinate_sample_index(self.ndim)
        

        # S3 sample the left-right boundary pairs to 
        # decide the coordinate indexes of all the corners (boundaries)
                
        px=[p1_lst,p2_lst] 
        
        
        px=npy.array(px)
        px=npy.transpose(px)
        
        # px is a matrix of [ndim, 2]
        # for example, if we have 3  grid of (x,y,z)
        # px will be in the shape of 
        # X:   | x1,  x2|
        # Y:   | y1,  y2|
        # Z:   | z1,  z2|
        
        
        w1_lst=npy.array(w1_lst)
        w2_lst=1.0-w1_lst
        wx=[w1_lst,w2_lst]
        
        wx=npy.array(wx)
        wx=npy.transpose(wx)
        
        # wx is weighting matrix in the shape of px 
        
        
        idx=range(iax)
        

        
        grd_pts_lst=list()
        wgt_lst=list()

        # S4 Generate a  list of 2^ndim corners, 
        #  each element is for a coordinate in ndim-space
        
        
                     
        for c_cor in cor_lst:
            
            tx=(idx, c_cor[idx])
            
            # idx--> which row  (axis id)
            # c_cor[idx]: 0 (left) or 1 (right) for column 
             
            # T coordiante and weight for one corner
            pt=px[tx]
            wgt=wx[tx]
            pt=tuple(pt)
            
            total_wgt=1.0
            
            for sgl_wgt in wgt:
                total_wgt=total_wgt*sgl_wgt
            
            grd_pts_lst.append(pt)
            wgt_lst.append(total_wgt)
            
            
        return grd_pts_lst, wgt_lst
    
    def get_axis_ids(self, axis_name_lst):
        
        """ Find id for given names 
        
        Inputs:
        ----------------------------------------------------------------
        1. axis_name_lst:<list, t:str>: the names
        
        Returns:
        ------------------------------------------
        1. axis_id_lst:<list,t:integer>: the index
        """
        axis_id_lst=[]
        for name in axis_name:
            idx=self.__axis_dict[name]
            axis_id_lst.append(idx)
    
        return axis_id_lst
    
    
    def allocate_point_subspace(self, p, axis_id_lst):
       
        """ find grid points around a given point, and  
        calculate their weights according to the distances from this point.  
        
        Inputs:
        ----------------------------------------------------------------
        1. p:<array, (nuse,)>: the coordiantes of the given point.  
        2. axis_id_lst:<list, t:integer>: the use the index
        
        Returns:
        -------------------------------------------------------------------
        1. grd_pt_lst: <list, t:tuple>: list of the index of the surrounding grid points.
        ---If the dimension size nax is n, the length of grd_pt_lst is 2^n, 
        with each element being an array of size n with coordinate [:] for axis=axis_id
        
                
        2. grd_wgt_lst: <list, t:float>: Weights of the the surrounding points,         
        ---which are calculated from their distances from the point.
        
        """
        
        
        p1_lst=list() # list for left-boundary
        p2_lst=list() # list for right-boundary
        w1_lst=list() # list for weighting 

        vtype=oob.get_ot_type(axis_id_lst)

        if (vtype==oob.ot_list):
            ndim_ax=len(axis_id_lst)
            
        elif (vtype==oob.ot_array):
            ndim_ax=npy.size(axis_id_lst)
            
        else:
            axis_id_lst=[axis_id_lst]
            ndim_ax=1
            
        
        vtype=oob.get_ot_type(p)
        
        if (vtype==oob.ot_list):
            
            ndim_p=len(p)
            if (ndim_p<>ndim_ax):
                msg='number of axis is different from the size of  given coordiante'
                msm.show_err_msg(msg)
                return None, None
        
        elif (vtype==oob.ot_array):
            ndim_p=npy.size(p)
            if (ndim_p<>ndim_ax):
                msg='number of axis is different from the size of  given coordiante'
                msm.show_err_msg(msg)
                return None, None

        else:
            p=[p]
            ndim_p=1
            if (ndim_p<>ndim_ax):
                msg='number of axis is different from the size of  given coordiante'
                msm.show_err_msg(msg)
                return None, None

        
        

        # S1 get position and weights for boundary points along each axis
        
        iax=0
        iuse=0
        
        
        for iax in range(self.ndim):
            
            if (iax in axis_id_lst):
                val=p[iuse]
                gax=self.__axis_lst[iax]
                p1, p2, wgt=gax.getwgt(val)
                p1_lst.append(p1)
                p2_lst.append(p2)
                w1_lst.append(wgt)
                iuse=iuse+1
        
                
        
        # S2 coordiante indexes for the corners (boundaries)
             

        
        
        # S3 sample the left-right boundary pairs to 
        # decide the coordinate indexes of all the corners (boundaries)
                
        px=[p1_lst,p2_lst] 
        
        
        px=npy.array(px)
        px=npy.transpose(px)
        
        # px is a matrix of [ndim, 2]
        # for example, if we have 3  grid of (x,y,z)
        # px will be in the shape of 
        # X:   | x1,  x2|
        # Y:   | y1,  y2|
        # Z:   | z1,  z2|
        
        
        w1_lst=npy.array(w1_lst)
        w2_lst=1.0-w1_lst
        wx=[w1_lst,w2_lst]
        
        wx=npy.array(wx)
        wx=npy.transpose(wx)
        
        # wx is weighting matrix in the shape of px 
        
        
        

        
        grd_pts_lst=list()
        wgt_lst=list()
        # S4 Generate a  list of 2^ndim corners, 
        #  each element is for a coordinate in ndim-space
        
        idx=range(iuse)
        cor_lst= coordinate_sample_index(iuse)
        # cor_lst is a list of left-rigt (0 or 1) coordiantes at ndim-shace, 
        # its length is 2^ndim 
        
        all_sel=slice(None, None, None) # [:]
        
        for c_cor in cor_lst:
            
            tx=(idx, c_cor[idx])
            
            # idx--> which row  (axis id)
            # c_cor[idx]: 0 (left) or 1 (right) for column 
             
            # T coordiante and weight for one corner
            pt=px[tx]
            
            
            wgt=wx[tx]
            
            total_wgt=1.0
            
            for sgl_wgt in wgt:
                total_wgt=total_wgt*sgl_wgt
            
            # T project pt to full space 

            fpt=[] # full coordinates
            
            iz=0 # index inside pt
            
            # print 'pt', pt
            
            for ix in range(self.ndim):
                
                if (ix in axis_id_lst):
                    
                    ipz=pt[iz]
                    fpt.append(ipz)
                    iz=iz+1
                
                else:
                    
                    fpt.append(all_sel)
                    
                    
            
            fpt=tuple(fpt)
            grd_pts_lst.append(fpt)
            wgt_lst.append(total_wgt)
        
        
        
        # re-organize the tuple with coordiante=[:] for axis=axis_id
        
            
        return grd_pts_lst, wgt_lst



    
    def get_value_at_point(self, p, field):
        
        """ calculate value at a given p by interpolating the 
        gridded data (field) at surrounding grid points
        
        Inputs:
        ----------------------------------------------------------------
        1. p:<array, (nax,)>: the coordiantes of the given point.  
        2. field:<array, (nx, ny, nz, ...)>:gridded data
        
        
        Returns:
        ----------------------------------------------------------------
        val:<float>: values at point p
        """
        
        if (oob.get_ot_type(p)==oob.ot_list):
            if (len(p)<>self.ndim):
                msm.show_err_msg('Coordinate of the point has a size different from the grid')
                return None
        
        else:
            if (npy.size(p)<>self.ndim):
                msm.show_err_msg('Coordinate of the point has a size different from the grid')
                return None

        
        grd_pts_lst, wgt_lst=self.allocate_point(p)
        
        
        icp=0 #
        val=0.0
    
        # loop over each corner (boundary) point
        
        for wgt in wgt_lst:
            xp=grd_pts_lst[icp]
            val=val+wgt_lst[icp]*field[xp]
            icp=icp+1
        return val
    
    
    def get_profile(self, p, field, axis_id_lst):
        
        
        """ extract profile at a  given 'horizontal location'
        
        Inputs:
        ----------------------------------------------------------------
        1. p:<array, (nuse,)>: the sub coordiantes of the given points.  
        2. field:<array, dims>: gridded data 
        3. axis_id_lst:<integer>: index of the axis, where the axises to be used
        
        Returns:
        ---------------------------------------------------------------
        
        1. profile:<array, ...>: the profile at the locations
        
        """
        
        

        
        
        #S1 get coordiante indexes and weight for the corners around the point 
        

        pt_lst, wgt_lst=self.allocate_point_subspace(p, axis_id_lst)
        
        #S2 interpolate the data to it 

        for ip in range(len(pt_lst)):
            cpt=pt_lst[ip]  # full coordinate for the corners
            twgt=0.0        # total weigt
            if (ip==0):
                profile=wgt_lst[ip]*field[cpt]
                twgt=wgt_lst[ip]
            else:
                profile=profile+wgt_lst[ip]*field[cpt]
                twgt=twgt+wgt_lst[ip]
        
        return profile
    
        
    def copy_attr_dict(self):
        """ 
        get one copy of the attr_dict
        """
        return dict(self.attr_dict)
    
    def shape(self):
        """ Return the shape of the grid 
        """
        return self.dims
    
    def size(self,ax_id=None):
        """
        Return size along one direction or total size of the grid 

        Inputs:
        --------------------------------
        1. ax_id:<integer/str>: Id for the axis to be found
        
        Returns:
        1. nsize:<integer>: the size of axis if axis_id is provided, or the total size if no ax_id is specified. 
        
        
        """
        
        nsize=1
        
        if (ax_id==None):
            
            for idim in self.dims:
                nsize=nsize*idim
                
            return nsize
        

        
        if (type(ax_id)==type('name')):
            if (ax_id in self.__axis_dict):
                ls_id=self.__axis_dict[ax_id]
                nsize=self.dims[ls_id]
                return nsize
            
        
            else:
                
                msg='axis '+ax_idx+' not exist'
                msm.show_err_msg(msg)
                return None
            
        elif (ax_id<self.ndim):
            nsize=self.dims[ax_id]
            return nsize
        else:
            msg='Axis Id number is larger than the grid dimensions'
            msm.show_err_msg(msg)
            return None
        
            
        return nsize

        
    
    
    def copy(self):
        
        """ 
        make a copy of itself
        """
        
        attrs=self.copy_attr_dict()
        new_grid=gp_grid_cl(self.__axis_lst, **attrs)
        
        
            
        return new_grid
    
        
## >S2 functions 

def build_grid_from_matrix(ax_name_lst, \
                               ax_colno_lst, \
                               data, zname=None, zval=None,\
                               mask_val=oob.fill_val):
    
    """convert a matrix into a gridded data set
    
    
    Inputs:
    -----------------------------------------------------
    1. ax_name_lst:<list, t:str>: axis names
    2. ax_colno_lst:<list, t:integer>: columns of the axis in the table 
    3. data:<array, (nrow, ncol>: data
    --- if (zname and zval) are given,   data is supposed to be in the column of 
    --- [ax_1, ax_2, ax_3,.., ax_n, z1, z2, z3, ...]
    
    ---if (zname and zval) are not given, data is supposed to be in the form of 
    --- [ax_1, ax_2, ax_3, ...,ax_n, val]
    
    4. zname:<str>: name of the z-axis ('vertical' axis)
    5. zval:<array>:  z ('vertical') data: 
    


    Returns:
    --------------------------------------------
    1. grd:<gp_grd_cl>: the grid defined by  [ax_name_lst, yname]
    2. data:<array, (nrow, nax)>: gridded data
    
    """
    
    
    # S1  split axis column from the data 
    
    nrow, ncol=npy.shape(data)
    val_lst=list()
    dcol=0
    
    
    for icol in range(ncol):
        if (icol in ax_colno_lst):
            pass
        else:
            val_lst.append(data[:, icol])
            dcol=dcol+1
    
    if (dcol==1):
        # simple data set 
        val_lst=npy.array(val_lst)
        val_lst=npy.squeeze(val_lst)
        
    else:
        
        # transpose
        val_lst=npy.array(val_lst)
        val_lst=npy.transpose(val_lst)
        
    
    # S2 setup  axis to add
    
    
    
    nax=len(ax_colno_lst)
    
    dims=list()
    all_axis=list()
    
    # print ax_name_lst
    # print ax_colno_lst
    
    for iax in range(nax):
        
        ax_nm=ax_name_lst[iax]
        axcol=ax_colno_lst[iax]
        # print axcol
        ax_val=data[:, axcol]
        
        
        ## removing repeated values from axis
        
        
        ax_val=set(ax_val)
        ax_val=list(ax_val)
        ax_val=npy.sort(ax_val)
        
        axis_new=axis_m.gp_axis_cl(ax_nm, ax_val)
        dims.append(npy.size(ax_val))
        
        all_axis.append(axis_new)
    
    #c#c for 
        
        
    if (zname<>None):
        ## the row data is given for 'y' axis
        
        axis_z=axis_m.gp_axis_cl(zname, zval)
        all_axis.append(axis_z)
        
        nz=npy.size(zval)
        dims.append(nz)
        
        if (nz<>ncol-nax):
            msg=r'size is not much: %i vs %i' % (ncol-nax, nz)
            msm.show_err_msg(msg)
            grd=gp_grid_cl(all_axis)
            return grd, None
        
        
        # outdata is in shape of [nax1, nax2, .., naxn, ny]
        
        outdata=npy.zeros(dims, float)
        
        

        all_sel=slice(None, None, None)
        
        for irow in range(nrow):
            xcor=[]
            for iax in range(nax):
                ax_sel=all_axis[iax]
                axcol=ax_colno_lst[iax]
                val=data[irow, axcol]
                
                # get the position 
                pos=ax_sel.get_closest_point(val, mask_val)
                
                xcor=xcor+[pos]
            
            
            
            xcor.append(all_sel)
            xcor=tuple(xcor)
            outdata[xcor]=val_lst[irow,:]
        
            
            
            
        
    else:
        
        # data_lst is single column
        outdata=npy.zeros(dims, float)
        
        for irow in range(nrow):
            xcor=tuple()
            for iax in range(nax):
                ax_sel=all_axis[iax]
                axcol=ax_colno_lst[iax]
                val=data[irow, axcol]
                pos=ax_sel.get_closest_point(val, mask_val)
                xcor=xcor+tuple([pos])
    
            
            outdata[xcor]=val_lst[irow]
        
    grd=gp_grid_cl(all_axis)
    
    return grd, outdata

            
if (__name__=='__main__'):
     
    print '> test 1: form a grid 10(lon)x10(lat) with 47 levels <<<'
    
    rlon=npy.arange(-180, 180.0, 10)
    axis_lon=axis_m.gp_axis_cl('lon', rlon, ax_unit='Deg')
    
    
    rlat=npy.arange(-90, 90.0, 10.0)
    axis_lat=axis_m.gp_axis_cl('lat', rlat, ax_unit='Deg')
    
    
    z=npy.arange(47)
    axis_z=axis_m.gp_axis_cl('z', z, ax_unit='None')


    axis_list=[axis_lon, axis_lat, axis_z]

    mod_grid=gp_grid_cl(axis_list)
    print 'dims:', mod_grid.dims
    print 'ndim:', mod_grid.ndim
    
    
    print ' '
        
    
    print '>> test 2: find the axis  of lon  <<<'
    gax=mod_grid['lon']
    
    print 'axis name:', gax.ax_name
    print 'axis value:', gax[:]
    print ' '
    
    
    print '>>> test 3: find axis of No. 1  <<<'
    gax=mod_grid[1]
    
    print 'axis name:', gax.ax_name
    print 'axis value:', gax[:]

    
    print ''
     

    print '>>>> test 4: set attribute  <<<<'
    mod_grid.set_attr('z-type', 'model level')
    
    tz=mod_grid.get_attr('z-type')
    
    print 'attribute value for z-type:', tz
    print ' '
    

    print '>>>> test 5: find the grid points surrounding point [35.0, 35.0, 0] <<<<'  
    
    px=[35.0, 35., 0]
    data=npy.zeros(mod_grid.dims, float)
    
    for ix in range(mod_grid.dims[0]):
        for jy in range(mod_grid.dims[1]):
            lon=rlon[ix]
            lat=rlat[jy]
            
            data[ix, jy,:]=lon+lat
    
   
    
    
    pt_lst, wgt_lst=mod_grid.allocate_point(px)
    
    print 'px', px
    print 'neighbouring grid points:', pt_lst
    print 'weights:', wgt_lst
    val=0.0
    
    for ip in range(len(pt_lst)):
        xt=pt_lst[ip]
        wt=wgt_lst[ip]
        vx=data[xt]
        
        val=val+wgt_lst[ip]*vx
        
    print 'Interpolated value:', val
    print  'True value:',px[0]+px[1]
    val2=mod_grid.get_value_at_point(px, data)
    print 'get_value_at_point', val2
    print ''
    
    # get profile 
    
    print '>>>>>> test 6: allocate_point_subspace:'
    pt=[32, 0]
    print pt

    pt_lst, wgt_lst=mod_grid.allocate_point_subspace(pt, [0, 1])
    print 'found corners:', pt_lst
    print 'weights:', wgt_lst
    
    accurate_val=pt[0]+pt[1]
    twgt=0.0
    for ip in range(len(pt_lst)):
        xt=pt_lst[ip]
        twgt=0.0
        
        if (ip==0):
            profile=wgt_lst[ip]*data[xt]
            twgt=wgt_lst[ip]
        else:
            profile=profile+wgt_lst[ip]*data[xt]
            twgt=twgt+wgt_lst[ip]

    
    print 'interpolated profile', profile
    print 'True value', accurate_val
    
    
            
    print '>>>>>>>> test 7: get sections at:'
    
    pt=32.0
    
    pt_lst, wgt_lst=mod_grid.allocate_point_subspace(pt, 0)
    
    print 'corners', pt_lst
    print 'weights', wgt_lst
    
    accurate_val=pt
    twgt=0.0
    for ip in range(len(pt_lst)):
        cpt=pt_lst[ip]
        twgt=0.0
        
        if (ip==0):
            profile=wgt_lst[ip]*data[cpt]
            twgt=wgt_lst[ip]
        else:
            profile=profile+wgt_lst[ip]*data[cpt]
            twgt=twgt+wgt_lst[ip]
            
    
    print 'interpolated field:' 
    print  'shape:', npy.shape(profile)
    print  'value:', profile[:,0]
    print 'True value:', accurate_val+rlat
    
    
    #>> test get profile
    print '>>>>>>>>>> test8:   get_profile '

    pt=45
    profile=mod_grid.get_profile(pt, data, 1)
    print  'shape:', npy.shape(profile)
    print  'value:', profile[:,0]
    print 'True value:', pt+rlon
    
    


    print '>>>>  build grid from a xy-z table'
    

    a=npy.arange(100.0)
    

    a=npy.reshape(a, [10, -1])
    
    ax1=[0]*5+[0.1]*5
    ax2=[0, 1.0, 2.0, 3.0, 4.0]*2
    
    
    
    ax1=npy.array(ax1)
    
    ax2=npy.array(ax2)
    
    a=npy.column_stack((ax2[:, npy.newaxis], a))
    a=npy.column_stack((ax1[:, npy.newaxis], a))
    zname='z'
    zval=npy.arange(10.0)
    
    
    print 'data shape', npy.shape(a)
    print 'ax1', a[:,0]
    print 'ax2', a[:,1]
    print 'zval', zval[:]
    
    grd, data=build_grid_from_matrix(['x', 'y'], \
                                         [0,1], \
                                         a, zname=zname, zval=zval)
    print 'after build grid' 
    print 'shape of data:' , npy.shape(data)
    print 'data at [1,0,:]:',data[1, 0, :]
    print  'ax1, ax2 at:',   ax1[0], ax2[1]
    print  'original data',  a[5, :]
    

    
    print '>>>> build grid from a xy table'

    a=npy.arange(10.0)
    
    ax1=[0]*5+[0.1]*5
    ax2=[0, 1.0, 2.0, 3.0, 4.0]*2
    
    
    
    ax1=npy.array(ax1)
    print 'ax1:', ax1[:]
    
    ax2=npy.array(ax2)
    print 'ax2:', ax2[:]
    
    a=npy.column_stack((ax2[:, npy.newaxis], a[:, npy.newaxis]))
    a=npy.column_stack((ax1[:, npy.newaxis], a))
    
    print 'shape of the table:', npy.shape(a)
    
    zname=None
    zval=None
    
    
    grd, data=build_grid_from_matrix(['x', 'y'], \
                                         [0,1], \
                                         a, zname=zname, zval=zval)
    
    print 'after build grid'
    print 'shape of data:', npy.shape(data)
    
    print 'data at [1,:]', data[1, :]
    print 'data at [:,1]', data[:, 1]
    print  'original data', a 
    
   
    
                
    
