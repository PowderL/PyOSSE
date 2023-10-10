""" (class ) for two or three or four dimenisonal (lon-lat-press-time etc) data set

Authors: L. Feng, Edinburgh University
History: v0.5, 2012.06.28
History: v0.95, 2012.10.28
this class  defines a grid to represent fields, and 
routines for interplotating fields. 


""" 

import numpy as npy
# from scipy import *

import geo_constant as gc
import gp_axis_m as gax
import message_m as msm
import gp_grid_m as gp_grid
     
class gp_data_cl:
    """ general container for gridded geo-physical variables
    
    Members:
    ------------------------------------------------------------------------------
    
    1. name: <str>: name of variable 
    2. category:<str>:category of variable 
    3. id: <numberic>:id of variable
    4. group:<str>:  group name of the variable 
    4. grid:<gp_grid_cl>:grid
    5. __attr_dict:<dict>: dictionary for name and ndim
    6. data
    
    Functions:
    ==================================================================
    # Interpolation and index
    1. gp_interp:interpolation data set to a given location
    2. gp_get_prof:get profile at given locations
    3. gp_interp_2d: a fast version for interpolate 2-D data set to given locations
    4. gp_axis_transform_2d: bilinear interoplation
    5. interp_along_axis:get (n-1)-Dimensional fields at locations along one chosen axis. 
    6. mylice: get indexes to access chosen elements of the data.  
    7. get_slice: 
    
    
    8. save_to_netcdf:save data to a given netcdf file
    9. set_global_attr: set one global attribute
    10. get_global_attr: get one global attribute
    11. remove_global_attr:remove one global attribute
    
    
    
    # disk IO
    11. save_gp_to_netcdf: save selected data set to one netcdf 
    12. load_gp_from_file: virtual function to be overrided
    
    
    
    """ 
    
    def __init__(self, data, \
                     name="", \
                     axis_set=[],  \
                     id=0, \
                     category="", \
                     unit="",\
                     data_attr={}):
        
        """ initialization
        
        Inputs:
        ----------------------------------------------------------
        1. data:<array>: array of the data set
        2. axis_set:<list>: list of (gp_)axis for grid
        3. name:<str>: name of the data set
        3. category:<str>: category of the data set
        4. unit:<str>:unit of the data set
        5. data_id:<obj>:id of the data set
        7. attr:<dict>: dictionary of the attributes
        
        """ 
        
        
        self.data=data
        self.name=name
        self.dims=npy.shape(self.data)
        
        
        self.axis_set=axis_set
        self.naxis=len(axis_set)
        if (self.naxis==npy.size(dims)):
            self.set_grid(grid_set)
            
        else:
            self.grid=None
            self.is_grid_ready=False
        
            
            
        self.__attr_dict=dict()
        self.id=id
        
        for attr_name in attr:
            val=attr[attr_name]
            self.__attr_dict.update({attr_name:val})
            if (attr_name=='name'):
                self.name=val
            elif (attr_name=='unit'):
                self.unit=val
            elif (attr_name=='category'):
                self.category=category
            elif (attr_name=='id'):
                self.id=val

    def get_attr(self, attr_name):
        """
        Inputs:
        ------------------------------------------------------------------------
        1. attr_name:<str>: name of attribute to be retrieved. 
        
        
        Returns:
        -------------------------------------------------------------------------
        1. val: <obj>: value of the attribute
        
        """
        if (attr_name in self.__attr_dict):
            return self.__attr_dict[attr_name]
        else:
            
            error_msg='No attribute name of '+attr_name
            msm.show_err_msg(error_msg)
            return None
    
    def set_attr(self, attr_name, val):
        """
        Inputs:
        ------------------------------------------------------------------------
        1. attr_name:<str>: name of attribute to be retrieved. 
        2. val: <obj>: value of the attribute
        
        """
        
        self.__attr_dict.update({attr_name:val})
        if (attr_name=='name'):
            self.name=val
        elif (attr_name=='unit'):
            self.unit=val
        elif (attr_name=='category'):
            self.category=category
        elif (attr_name=='id'):
            self.id=val
        elif (attr_name=='grid'):
            self.grid=val
        
        
    def __getitem__(self, idx):
        """
        Inputs:
        ------------------------------------------------------------------------
        1. idx:<tuple>: coordinate 
        
        Returns:
        -------------------------------------------------------------------------
        1. val: <obj>: value of the attribute
        
        """
        return self.data[idx]
    
    

    def set_grid(self, axis_list):
        self.grid=gp_grid_cl(axis_list)
        self.axis_list=axis_list
        if (len(self.axis_list)==npy.size(self.dims)):
            self.is_grid_ready=True
        else:
            self.is_grid_ready=False
    
    

    
        
    def get_gp(self, xcor, fill_val=None, threshold=0.0):  # pos as single point
        """ Interpolate data to a given location 
        
        Inputs:
        1.xcor: <array>: coordinates of the location
        2. fill_val:<float>:Optional. If data[idx] ==fill_val, 
        --- the weight of value at the axis point will be equal to zero
        3. threshold:<float>: if the total weight for xocs is smaller <threshold, its value 
        will be set to fill value  
        
        Returns:
        1. gp: variable value at the given locations
        2. wgt: weight information as (lp, rp, wgt) for every axis
        """
        
        
        axis_set=self.axis_set
        
        naxis=len(axis_set)

        if (len(rpos)<>naxis):
            return gax.gp_fill_value
        
        #S1 check the location and weights
        
        
        for ix in range(naxis):
            val=xcor[ix]
            ax=axis_set[ix]
            p1, p2,wgt=ax.getwgt(val)
            pos_info.append([p1, p2, wgt])
        
        out_val=0.0
        #S2 there are  2^naxis corners (or boundaries) around a given point
        
        npt=2**naxis
        sum_wgt=0.0
        
        #S3 Interpolate the values of all npt to the location
        #Note: 
        #1. if the value at one corner is invalid (i.e., ==fill_value), 
        #  its weigh will be  assigned to zeros
        
      
        for ii in arange(nval):
            pos=zeros(naxis)
            wgt=1.0
            it=ii
            ip=0
            ## construct corner coordinate and weights
            
            for ix in arange(naxis):
                pt=pos_info[ix]
                isel=mod(it,2)
                it=it/2
                pos[ix]=pt[isel]
                wgt=wgt*abs(isel-pt[2])
            
            ipos=tuple(pos)
            # print ii, ipos
            # print wgt
            
            if (gp[ipos]<>fill_val):
                # print gp[ipos], wgt
                sum_wgt=sum_wgt+wgt
                out_val=out_val+gp[ipos]*wgt
        if (sum_wgt>threshold):
            out_val=out_val/sum_wgt
        else:
            out_val=fill_val
        
        # if (naxis==1):
        #    print out_val, pos_info
        # del axis_set
        return out_val, pos_info

    
    def interp_2d(self, cor_lst): 
        """ Interpolate variable at 2D-grid values to a given locations in a two dimensional grid 
        
        Inputs:
        1. cor_lst: <list>: coordinate of the set of locations 
        (i.e., [[P0], [P1],...]
        
        Returns:
        1. gpv:<array>: variable values at the locations
        2. loc_info:<list>: list of weighs and locations of the corners information :
        [[lpx, rpx, wgx,lpy, rpy, wgy], ...]
        
        
        """
        
        #S1 check the number of points
        
        
        rset=npy.array(cor_lst)
        dims=npy.shape(rset)
        
        if (len(dims)==2): # the co-ordinate is given as [x, y]
            npt=1
        else:
            npt=dims[0]
        
            
        pos_info=list()
        gpv=zeros(npt, float)
        gp=self.data
        
        if (npt==1):
            rpos=rset[0]
            print 'rpos', rpos, shape(rpos)
            
            xp1, xp2, wgx=self.axis_set[0].getwgt(rpos[0])
            yp1, yp2, wgy=self.axis_set[1].getwgt(rpos[1])
            gpv=wgx*(wgy*gp[xp1,yp1]+(1.0-wgy)*gp[xp1, yp2])+ \
                     (1.0-wgx)*(wgy*gp[xp2,yp1]+(1.0-wgy)*gp[xp2, yp2]) 
            pos_info.append([xp1, xp2, wgx, yp1, yp2, wgy])
            return gpv, pos_info
        
        for ip in range(npt):
            rpos=cor_lst[ip]
            xp1, xp2, wgx=self.axis_set[0].getwgt(rpos[0])
            yp1, yp2, wgy=self.axis_set[1].getwgt(rpos[1])
            gpv[ip]=wgx*(wgy*gp[xp1,yp1]+(1.0-wgy)*gp[xp1, yp2])+ \
                     (1.0-wgx)*(wgy*gp[xp2,yp1]+(1.0-wgy)*gp[xp2, yp2])
            pos_info.append([xp1, xp2, wgx, yp1, yp2, wgy])           
            
        return gpv, pos_info
    
     
    def get_prof(self, xpos, xaxis=0, ):
        
        """ interpolate the data to a given locations [xvalx, yvals]
        Inputs:
        ------------------------------------------------------------------
        1. xpos: <array>: coordinate of the point. Its size is equal to ndim-1, where ndim         is the dimension size of the data set.  
        2. xaxis:<integer>: id of the axis for the profile to be retrieved 
        
        Returns:
        ------------------------------------------------------------------
        1. ref_x:<array>: grid along the chosen x
        2. gp_prof:<array>: profiles at the given locations 
        
        Notes:
        -------------------------------------------------------------------
        1. the stored data set  is assumed to have a dimension size of 
        ndim equal to len(xpos)+1
        
        """
        
        
        ip=0
        #S1 check out the grid for the given variable # 
        
        axis_set=self.axis_set
        ax=axis_set[xaxis]
        naxis=npy.size(axis_set)
        ref_x=ax[:]
        gp_prof=zeros(npy.size(ref_x), float)
        sum_wgt=array(gp_prof)
        pos_info=list()
        
        #S2 check left/right boundaries and interpolation weight 
        
        for ix in range(0, naxis):
            ## only for ix<>xaxis as the profile will be at direction of xaxis
            
            if (ix<>xaxis):
                val=xpos[ip]
                ax=axis_set[ix]
                p1, p2,wgt=ax.getwgt(val)
                pos_info.append([p1, p2, wgt])
                ip=ip+1
        
        #S3 get data set
        gp=self.data
        

        #S4 interpolation along different directions
        
        # for n-dimension, we have 2^n-1 grid point
        
        nval=2**(naxis-1)
        
        nx=size(ref_x)
        
        idx=arange(nx)
        
        for ii in arange(nval):
            pos_xl=()
            pos_xr=()
            wgt=1.0
            it=ii
            ip=0
            
            for ix in arange(naxis):
                if (ix<>xaxis):
                    pt=pos_info[ip]
                    isel=mod(it,2)
                    it=it/2
                    
                    if (ix<xaxis):
                        ## form the coordiante with 1 less than the data set dimension
                        
                        pos_xl=pos_xl+tuple([pt[isel]])
                        
                    if (ix>xaxis):
                        pos_xr=pos_xr+tuple([pt[isel]])
                        
                    wgt=wgt*abs(isel-pt[2])
                    ip=ip+1
            
            tuple_pos=self.myslice_sgl(pos_xl, None, pos_xr)
            
            # print tuple_pos
            # print shape(gp)
        
            gp_tmp=gp[tuple_pos]
            
            if (size(gp_tmp)>1):
                gp_tmp=npy.squeeze(gp_tmp)
                usd_idx=npy.where(gp_tmp<>fill_val)
                
            
                sum_wgt[usd_idx]=sum_wgt[usd_idx]+wgt
                gp_prof[usd_idx]=gp_prof[usd_idx]+wgt*gp_tmp[usd_idx]
                
            else:
                if (gp_tmp<>fill_val):
                    sum_wgt=sum_wgt+wgt
                    gp_prof=gp_prof+wgt*gp[tuple_pos]
                    
        
        #S3 only for wgt > threshold the point is valid
        
        if (fill_val<>None):
            gp_prof=where(sum_wgt>threshold, gp_prof/sum_wgt, fill_val)
        else:
            gp_prof=where(sum_wgt>threshold, gp_prof/sum_wgt, gax.gp_fill_val)
        #
        return ref_x, gp_prof 



    def save_gp_to_netcdf(self, gpfile):
        """ dump chosen data to a  netcdf file 
        Inputs:
        1. gpfile:<str>:file name
        2. gpname:<str>:variable name to be save

        Returns:
        
        """
        import netCDF_gen as nf
        
        dimTypes=list()
        dimVars=list()    
        dimNames=list()
        
        for ax in self.axis_set:
            dimTypes.append('f')
            axvar=array(ax[:])
            axname=ax.get_attr('name')
            dimNames.append(axname)
            dimVars.append(axvar)
            
        nf.netCDF_def_dims(gpfile,dimNames,dimTypes, dimVars)
            
        varName=self.name
        varType='f'
        varData=self.data
        
        nf.netCDF_var_write(gpfile,dimNames,varName, varType, varData)
    
    

    def myslice(self, pos_xl, idx, pos_xr):
        """ return a set of index to access elements of n-d array
        
        Inputs:
        -------------------------------------------------
        1. pos_xl, pos_xr: <tuple>: Tuple of index at the left and right sides of the chosen axis
        2. idx:<list>:  list of indexs along the the axis
        
        Returns:
        -------------------------------------------------
        ps: <list>: a list of tuple of the n-d index  
        
        """
        f=lambda ix:self.myslice_sgl(pos_xl, ix, pos_xr) 
        ps=map(f,idx)
        # ps=array(ps)
        del f
        return  ps
    
    def myslice_sgl(self, pos_xl, id, pos_xr):
        
        """ Form a myslice to access element of n-dimension array 
        Inputs:
        --------------------------------------------------------------------
        1. pos_xl, pos_xr:<tuple>: the tuple of index at the left and right 
        ---side of the chosen axis     2. id: <integer/None>: the location along the chosen axis
        
        Returns:
        -----------------------------------------------------
        1. sl: <tuple>: the tuple of the n-Dimension index. n=len(pos_xl)+len(pos_xr)+1
        If id==None, sl=tuple([pos_xl, :, pos_xr]). If not, sl=tuple([pos_x, id, pos_xr])
        
              
        """
        
        if (id==None):
            slp=tuple([slice(None, None, None)]) # slp=[:]
        else:
            slp=tuple([id])
        sl=pos_xl+slp
        sl=sl+pos_xr
        return sl
        
        # print sl
    
    
    def remove_attr(self, attr_name):
        """ remove attribute 
        1. Inputs
        attr_name:<str>: the name of the attribute to be removed 
        """
        
        if (attr_name in self._attr_dict):
            del self.__attr_dict[attr_name]
    
