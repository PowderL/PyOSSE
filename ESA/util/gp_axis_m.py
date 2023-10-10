""" 

Generic class for coordiante axis

    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.28
    History: v0.95, 2012.09.12
    
    this module define axis classes
    
    Classes:
    ======================================================================
    1. gp_axis_cl: class for generic axis
    


"""

import numpy as npy
import otool_obj as oob
import interpolate_f as ipf
import message_m as msm


         
class gp_axis_cl: 

    """ class axis 
    
    A. Members:
    ------------------------------------------
    1. ax_name:<str>: Name of the axis
    2. ax_grid:<array>: centers of grid boxes. 
    3. ax_edge:<array>: edges of the grid boxes
    4. unit:<str>: units of the axis
    5. ax_min, ax_max:<numeric>: ranges of axis
    6. attr_dict:<dictionary>: attributes of the axis
    7. ax_np:<int>: number of the points
    
    B. Functions:
    -------------------------------------------
    1.  __init__: initialize the class
    2.  __getitem__:get values at given index
    3.  set_attr:  assign one attribute to axis
    4.  get_attr:  get one attribute for axis
    5.  set_axis: set/reset axis intervals
    6.  get_axis: get axis intervals
    7.  set_unit: set axis units 
    8.  get_unit: get axis units
    9.  set_edge: set axis edge
    10. get_edge: get axis edges
    11. getpos:get left/right boundaries for given value(s) 
    12. getwgt:get boundaries and weighting factors for given value(s) 
    13. get_closest_point: find the closest grid point for an array or a point 
        
    

    """
       
    
    def __init__(self, ax_name, ax_grid=None, ax_edge=None, \
                     unit=None, **keywords):
        """ initilialization 
        
        Inputs:
        --------------------------------------------------------
        1. ax_name:<str>: String of the axis
        2. ax_grid:<array, (ax_np,)>: Centers of the (ax_np) grid boxes, 
        --- assumed to be in ascending order.
        3. ax_edge:<array, (ax_np+1,)>: edges of the grid box. 
        4. unit: <str>: unit of the grid values
        5. keywords:<dict>: attributes
        
        Returns:
        ----------------------------------------------
        """
        
        self.ax_name=ax_name
        
        # use dictionary as a container to store additional information
        
        self.attr_dict={'name':ax_name}
        self.unit=unit
        self.attr_dict.update({"unit":unit})
        self.attr_dict.update({"ot_type":oob.ot_axis})
        
        
        if (ax_grid==None):
            self.ax_grid=npy.arange(10)
        else:
            self.ax_grid=npy.array(ax_grid)
            
        
        self.ax_edge=ax_edge
        
        self.ax_min=npy.min(self.ax_grid)
        self.ax_max=npy.max(self.ax_grid)
        self.ax_np=npy.size(self.ax_grid)
        self.dims=npy.shape(self.ax_grid)
        self.ndim=npy.size(self.dims)
        
        for keyname in keywords:
            keyval=keywords[keyname]
            self.attr_dict.update({keyname:keyval})
        
    def __getitem__(self, idx):  
        
        """ override [],  etc called as axis[id]
        
        Inputs:
        --------------------------------------------------------------
        1. idx:<int/array>: index(es)
        
        Returns:
        -------------------------------------------------------------
        1. vals:<numeric/array>: val(s) for the given index(es)
        """ 
        
        vals=self.ax_grid[idx]
        return vals

    def shape(self):
        """
        get the shape of the axis)
        """
        return self.dims
    
    def size(self):
        """
        get the size of the axis)
        """
        return self.ax_np
    
    def getpos(self, val, mask_val=oob.fill_val, mask_outsider=False):
        
        """ get left-right node pairs (ie., boundaries) for given value(s) 
        
        Inputs: 
        ----------------------------------------
        1. val:<array, nval>: values to be allocated. 
        2. mask_outsider:<T/F>: 
        if switched on, locations for values outside the axis range will 
        be filled with oob.fill_val_int. Otherwise:
        ---if (val>ax_max): the left-right node pairs are set to the right-most
        ---if (val<ax_min): the left-right node pairs is set to the left-most
        
        Returns:
        ---------------------------------------
        1. lp, rp:<array, (nval,)>: the left and right axis nodes (i.e., boundaries)
        around the given value(s)
        
        
        """

        #S1. locate  upper boundaries
        
        # for values samller than ax_min, rp will be filled with 0
        # for values larger than ax_max,  rp will be filled with ax_np
        
        if (mask_outsider):
            lp, rp=ipf.getpos_mask_outsider(self.ax_grid, val, \
                                                mask_val=mask_val)
        else:
            lp, rp=ipf.getpos(self.ax_grid, val, \
                                  mask_val=mask_val)
        
        
        #S2. check the number of the array dimensions
            
        if (oob.get_ot_type(val)==oob.ot_array):
            return lp, rp
        else:
            # if val is secular 
            return lp[0], rp[0]
    
        
    def set_edge(self, vfirst=None, vlast=None):
        
        """ set edges of grid boxes 
        
        Inputs: 
        ----------------------------------------
        1. vfirst: <float>: if given, it will be used to replace ax_edge[0]
        2. vlast:  <float>: if given, it will be used to replace ax_edge[self.ax_np]
        
        Returns:
        ----------------------------------------
        1. self.ax_edge, <array, (ax_np+1,)>: edges of the grid boxes
        
        
        
        """
        # S1 set edges values 
        self.ax_edge=npy.zeros(self.ax_np+1, float)
        ax_right=self.ax_grid[1:]
        ax_left=self.ax_grid[0:self.ax_np-1]
        self.ax_edge[1:self.ax_np]=0.5*(ax_right+ax_left)
        #S2 correct the first and last one 
        
        if (vfirst==None): # no first value is given
            grd_hw=self.ax_edge[1]-self.ax_grid[0]
            self.ax_edge[0]=self.ax_grid[0]-grd_hw
        else:
            self.ax_edge[0]=vfirst
           
           
        if (vlast==None): # no first value is given
            grd_hw=self.ax_grid[-1]-self.ax_edge[-2]
            self.ax_edge[-1]=self.ax_grid[-1]+grd_hw
        else:
            self.ax_edge[-1]=vlast
            
        return npy.array(self.ax_edge)
    
    

            
    def get_edge(self):
        """ get edges of grid boxes 
        
        """

        return npy.array(self.ax_edge)
    


            
    def copy_attr_dict(self):
        """ get one copy of the attr_dict
        """
        return dict(self.attr_dict)
    
    
    def getwgt(self, val, mask_val=oob.fill_val, mask_outsider=False):
        
        """ get boundaries and weighting factors for given value(s) 
        
        Inputs: 
        ----------------------------------------
        1. val:<array, (nval,)>: values to be allocated. 
        2. mask_outsider:<T/F>: 
        if switched on, locations for values outside the axis range will 
        be filled with mask value. Otherwise:
        ---if (val>ax_max): the left-right pairs  are set to the right-most
        ---if (val<ax_min): the left-right pairs is set to the left-most
        
        Returns:
        ---------------------------------------
        1. lp, rp:<array, (nval,)>: the left and right axis nodes (i.e., boundaries)
        around the given value(s)
        
        2. wgt:<array, (nval)>: the weighting factors for the left boundary, being calculated by
        wgt=(val-ax_grid[lp])/(ax_grid[rp]-ax_grid[lp]
        
        
        """
         
    

        #S1. locate  upper boundaries
        
        # for values samller than ax_min, rp will be filled with 0
        # for values larger than ax_max,  rp will be filled with ax_np
        
        if (mask_outsider):
            lp, rp, wgt=ipf.getwgt_mask_outsider(self.ax_grid, val, mask_val=mask_val)
        else:
            # print val
            # print self.ax_grid
            
            lp, rp, wgt=ipf.getwgt(self.ax_grid, val, mask_val=mask_val)
            # print lp, rp, wgt
            
        
        
        
        #S2. check the number of the array dimensions
       
        if (oob.get_ot_type(val)==oob.ot_array):
            return lp, rp, wgt
        else:
            return lp[0], rp[0], wgt[0]
            
    
    
    
    def set_attr(self, attr_name, value):
        
        """  assign one attribute to axis
        
        Inputs:
        -------------------------------------------------------
        1. attr_name: <str>: attribute name
        2. value: <obj>:attribute value
        """
        
        self.attr_dict.update({attr_name:value})
        
        if (attr_name=='unit'):
            self.unit=value
        elif (attr_name=='name'):
            self.ax_name=name
        
    def get_attr(self, attr_name):
        """  get one attribute for axis
        Inputs:
        ------------------------------------------------
            1. name:<str>: attribute name
        
        Returns:
        -----------------------------------------------
            1. val: <obj>: attribute value
         """
        
        if (attr_name in self.attr_dict):
            val=self.attr_dict[attr_name]
        else:
            val=oob.fill_val
        
        return val
    
    def set_axis(self, ax_grid):
        
        """  set/reset axis grid 
        
        Inputs:
        -----------------------------------------------
        1. ax_grid:<array>:axis intervals 
        """

        self.dims=npy.shape(ax_grid)
        self.ndim=npy.size(self.dims)
       

        self.ax_grid=npy.array(ax_grid)

        self.ax_np=npy.size(ax_grid)
        self.ax_max=npy.max(ax_grid)
        self.ax_min=npy.min(ax_grid)
        
        
    def get_axis(self):
        
        """  retrieve axis grid 
        
        Returns:
        --------------------------------------------
        1. ax_grid:<array>:axis intervals 
        """
        

        return self.ax_grid

        
    def set_unit(self, unit):
        """  set axis unit 
        
        Inputs:
        ---------------------------------------
        1. unit:<str>:unit
        """
        self.unit=unit
    
    def get_unit(self):
        """  get axis unit 
        
        Returns:
        --------------------------------------
        1. unit:<str>:unit
        """
        return self.unit
    
                  
    def interpolate(self, f, x, mask_val=oob.fill_val, mask_outsider=False):
        
        """ interpolate f(x0) given at axis nodes to f(x) 
        
        Inputs: 
        ----------------------------------------
        1. f:<array, (ax_np, ...,) >: values at aixs nodes
        2. x:<array, nval>: values to be allocated. 
        3. mask_val:<float>: filling value for bad or missing values
        4. mask_outsider:<T/F>: 

        if switched on, locations for values outside the axis range will 
        ---be filled with mask value. Otherwise:
        ---if (val>ax_max): the left-right node pairs  are set to the right-most
        ---if (val<ax_min): the left-right node pairs is set to the left-most
        
        Returns:
        ---------------------------------------
        1. new_f:<array, (nx, ...,):values at x
        
        """
        # s1 check the size
        dims=npy.shape(f)
        
        ndim=npy.size(dims)
        nx=npy.size(x)
        new_dims=tuple([nx])
        ny=1
        
        if (not oob.check_ot_shape(f, self.dims)):
            msg=dims
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
        
        for idx in dims[1:]:
            new_dims=new_dims+tuple([idx])
            ny=ny*idx
        
        # S2 get location and weight of x at the axis 
        
        lp, rp, wgt=self.getwgt(x, mask_val=mask_val, mask_outsider=mask_outsider)
        # S3 interpolation
        
        if (ndim>1):
            # reshape f to [ax_np,-1]

            f=npy.reshape(f, [self.ax_np, -1])
            
            
            new_f=ipf.get_interpl_1d(f, lp, rp, wgt, mask_val=mask_val)
            # reshape back 
            
            new_f=npy.reshape(new_f, new_dims)
            if (npy.size(x)==1):
                if (not oob.get_ot_type(x)==oob.ot_array):
                    new_f=new_f[0]
        else:
            new_f=ipf.get_interpl(f, lp, rp, wgt, mask_val=mask_val)
            if (not oob.get_ot_type(x)==oob.ot_array):
                new_f=new_f[0]
        
        return new_f
    
    
    def get_closest_point(self, xval, mask_val, mask_outsider=False):
        
        """
        get the closet point to given vals
        
        Inputs:
        ----------------------------
        
        1. xval:<array, (nx,)>: values to be allocated
        2. masl_val:<float/int>: filling for bad or missing values
        3. mask_outsider:<T/F>: if True, All the outsider will be filled with mask_val. 
       
        Returns:
        ----------------------------------
        1. pt:<array, (nx,)>: the closest grid point 
        
        """
        
        p1, p2, wgt=self.getwgt(xval, mask_val, mask_outsider=mask_outsider)
        
        if (oob.get_ot_type(xval)==oob.ot_array):
            # if wgt >0.5 wgt=1.0
            chose_idx=(p1<>int(mask_val)) & (wgt<>mask_val) & (wgt>=0.5)
            
            wgt=npy.where(chose_idx, 1.0, wgt)
            
            # if (wgt<0.5, wgt=0.0
            chose_idx=(p1<>int(mask_val)) & (wgt<>mask_val) & (wgt<0.5)
            wgt=npy.where(chose_idx, 0.0, wgt)
        
            # only value with valid index will be used
            
            chose_idx=npy.where((p1<>int(mask_val)) & (wgt<>mask_val))
            
            pt=1.0*p1
            
            pt[chose_idx]=wgt*p1[chose_idx]+(1.0-wgt)*p2[chose_idx]
            pt=pt.astype(int)
            
        else:
            
            if (p1==int(mask_val)):
                pt=mask_val
            else:  
                if (wgt>=0.5):
                    pt=p1
                else:
                    pt=p2
        
        return pt
    

    def copy(self):
        
        """
        make a copy of itself without shared instance
        """
        

        attrs=self.copy_attr_dict()
        if ('name' in attrs):
            del attrs['name']
        
        if ('unit' in attrs):
            del attrs['unit']
        
        
        new_axis=gp_axis_cl(self.ax_name, self.ax_grid, \
                                ax_edge=self.ax_edge, \
                                unit=self.unit,\
                                **attrs)
        return new_axis
    




             
        
    
 #<<< TESTS >>>
                    
        
if (__name__=='__main__'):
    """ tests of  gp_axis """
    import numpy.random as rnd
    print ">>>t1: build a axis for (0, 6.0,1.0) with name='a', unit='ppm' <<<"
    
    a=npy.arange(6.0)
    test_gx=gp_axis_cl('a', a, unit='ppm')
    print 'set edge'
    ax_edge=test_gx.set_edge()
    print 'center:', test_gx[:]
    print 'edge:', ax_edge
    
    print ax_edge
    print 'axis name:', test_gx.ax_name
    print 'axis unit:', test_gx.unit
    print 'axis intervals:', test_gx[:]
    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    print '>>>t2: interpolation-test (array)<<<' 
    
    c=a+rnd.rand(6)
    
    print 'original array:', c
    
    p1, p2=test_gx.getpos(c)
    
    print 'allocate the array into axis'
    print 'left bd:', p1
    print 'right bd:', p2
    
    p1, p2,wgt=test_gx.getwgt(c)
    print 'weight (left):', wgt
    
    c2=wgt*test_gx[p1]+(1-wgt)*test_gx[p2]
    print 're-construct array'
    print 'reconstructed array:', c2
    
    
    print 'mask outsiders'
    print 'original array:', c
    p1, p2,wgt=test_gx.getwgt(c, mask_outsider=True)
    p1, p2=test_gx.getpos(c, mask_outsider=True)
    print 'left bd:', p1
    print 'right bd:', p2
    print 'weight (left):', wgt
    

    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    print '>>>t4: interpolation-test (single value)<<<' 
    
    c=3.5
    
    print 'original val:', c
    p1, p2,wgt=test_gx.getwgt(c)

    p1, p2=test_gx.getpos(c)
    
    print 'left bd:', p1
    print 'right bd:', p2
    print 'weight (left):', wgt
    
    
    c2=wgt*test_gx[p1]+(1-wgt)*test_gx[p2]
    print 're-constructed val:', c2
    
    print 'mask outsider'
    c=-2.0
    print 'original val:', c
    p1, p2,wgt=test_gx.getwgt(c, mask_outsider=True)
    p1, p2=test_gx.getpos(c, mask_outsider=True)

    print 'left bd:', p1
    print 'right bd:', p2
    print 'weight (left):', wgt
    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    
   
    print '>>>t5: reset the axis to (0., 7.0, 1.0)<<<'
    b=npy.arange(7.0)
    # reset the axis
    test_gx.set_axis(b)
    
    print 'new points:', test_gx.ax_np
    print 'new min value:', test_gx.ax_min
    print 'new max value:', test_gx.ax_max
    
    print 'get axis:'
    t_axis=test_gx.get_axis()
    print t_axis
    
    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    

    print '>>>t6: reset unit to ppb'

    test_gx.set_unit('ppb')

    print 'get unit'
    ux=test_gx.get_unit()
    print ux

    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    
    
    
    print ">>>t7: set attribute {'model':GEOS5}"
    
    test_gx.set_attr('model', 'GEOS5')
    
    print "get attribute"
    
    
    attr=test_gx.get_attr('model')
    
    print attr

    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    
    
    
    print ">>>t8: test interpolation"

    print '******* 1D data'
    print 'axis', test_gx.dims
    
    f=npy.arange(test_gx.ax_np)
    print 'f', f
    
    new_x=npy.arange(4.0)
    new_x=new_x+0.2
    new_f=test_gx.interpolate(f, new_x)
    print 'new x', new_x
    print 'new f', new_f
    
    
    new_x=0.2
    new_f=test_gx.interpolate(f, new_x)
    
    print 'new x', new_x
    print 'new f', new_f

    
    print '****** 2D data'
    
    print 'axis', test_gx.dims
    
    f=npy.arange(3*test_gx.ax_np)
    f=npy.reshape(f,[test_gx.ax_np, -1])
    print 'f', f
    
    new_x=npy.arange(4.0)
    new_x=new_x+0.2
    new_f=test_gx.interpolate(f, new_x)
    print 'new x', new_x
    print 'new f', new_f
    
    
    new_x=0.2
    new_f=test_gx.interpolate(f, new_x)
    
    print 'new x', new_x
    print 'new f', new_f
 
    
    
    print '<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>'
    print '' 
    
    new_x=npy.array([0.2,2.2])
    
    pt=test_gx.get_closest_point(new_x, -999.0)
    print pt
    print new_x
    
    
    
    new_x=npy.array([2.7])
        
    pt=test_gx.get_closest_point(new_x, -999.0)
    print pt
    print new_x
    
    
    new_x=npy.array([2.7, -30.0])
    
    pt=test_gx.get_closest_point(new_x, -999.0)
    print pt
    print new_x
    
    
             

      
