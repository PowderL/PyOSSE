""" Class for vertical interpolation of model profiles.  
   
    Authors: L. Feng, Edinburgh University
    History: v0.9, 2012.06.30
    History: v0.95, 2012.12.17
    
    1. Classes:
    ==============================================
    1. vintpl_cl: class for vertical interpolation
    

"""
import numpy as npy
import gp_axis_m as axis_m
import vertical_profile as vpf
import vertical_column as vcol

import sample_model_field as grd_pf
import message_m as msm
import otool_obj as oob



class vintpl_cl:
    
    """class for vertical interpolation and integration 
    
    Members:  
    -------------------------------------------------------------------------------
    
    # model grid 
    1.  grd_pres:<array, ([nx], [ny], ..., nz)>:  vertical pressure levels at horizontal locations. 
    ---It can be 1D, 2D, 3D or more, but the last dimension must be vertical levels. 
    
    2.  lg_grd_pres:<array, ([nx], [ny], ..., nz)>: log10(grd_pres)
    
    3.  grd_dims:<tuple/array>: shape of lg_grd_pres
    4.  grd_ndim:<integer>: length of grd_dims
    
    5.  do_reverse:<T/F>: True if the original grid pressure is in descending order. 
    --- class member grd_pres will  always be stored in ascending order.
    
    # interpolation 
    
    6.  grd_vpl:<array, ([nx], [ny], nz_ob)>: lower boundary for interpolations from  
    ob (target) grid to model grid
    
    7.  grd_vpr:<array, ([nx], [ny], nz_ob)>: upper boundary for interpolations from  
    ob (target) grid to model grid
    
    8.  grd_vwgt:<array, ([nx], [ny], nz_ob)>: weighting factor for lower boundary 
    
    # column
    9.  grd_colwgt<array, , ([nx], [ny], nz)>: mass weight at levels 
    
    
    # Observation (target) grid 
    
    10. ob_pres:<array, ([nx], [ny], ..., nz_ob)>:  ob pressue levels 
    
    11. lg_ob_pres:<array, ([nx], [ny], ..., nz_ob)>: log10(lg_ob_pres)
    
    12. ob_dims:<tuple>: shape of lg_grd_pres
    
    13. ob_ndim:<integer>: length of ob_dims
    
    
    14. ob_vpl:<array, ([nx], [ny], ..., nz)>: lower boundary for interpolations from  model grid to ob grid
    15. ob_vpr:<array, ([nx], [ny], ..., nz)>: upper boundary
    16. ob_vwgt:<array, ([nx], [ny], ..., nz)>: weighting factor for lower boundary
    17. do_ob_reverse:<T/F>: Ture when  the original ob_pressure is in descending order
    
    # others 
    
    18. mask_val: <float>: mask for missing or bad data
    19. do_reshape:<T/F>:  True when the dimension of the pressure is larger than 3
    --- when it is true, the lg_grd_pres will  be reshaped as (reshape_nx, nz) 
    --- before calling profile or column fortran functions.
    --- the result will be reshaped back to high-dimension data
    
    20. reshape_nx:<integer>: the size of the leading column
    
    
    
    
    
    Functions:
    -------------------------------------------------------------------------------
    1.  __init__: initialization
    2. init_interp: set parameters for interoplation 
    3. interpolate_ob_prof:  Interpolate ob profiles to model pressure grid
    4. interpolate_mod_prof:  interpolate field from model grid to ob grid 
    5. get_mod_column: calculate column values for profiles at model vertical grid
    6. get_ob_column: calculate column values for profile  at ob vertical grid
    7. get_ob_ak_column: calculate column values for profiles weighted by averging kernel at ob vertical grid. 
    
    
    
    
    """
    
    def __init__(self, pres, is_in_log, do_reverse, mask_val=oob.fill_val):
        
        """  initialize horizontal interpolation class 
        Inputs:
        1. pres:<array, ([nx], [ny], nz>: 
        ---Pressure grid for profiles to be interpolated to or from.     
        ---The last column of the array is assumed to be vertical levels 
                
        2. is_in_log: <T/F>: Ture the pressure is given in log10 
        3. do_reverse:<T/F>: True, if the pressure is given in decending order/ 
        4. mask_val: <float>: masks for bad or missing value 
        
               
        
        """
        # masks for missing or bad data
        
        self.attr_dict={}
        self.attr_dict.update({'ot_type': oob.ot_vintpl})
                
        self.mask_val=-999.0
        self.set_mask(mask_val)
        
        
        dims=npy.shape(pres)
        ndim=npy.size(dims)
        nz=dims[-1]
        
        self.grd_dims=dims
        self.grd_ndim=ndim
        self.reshape_nx=1
        
        self.do_reverse=do_reverse
        self.do_reshape=False
        
        if (do_reverse):

            if (ndim==1):
                pres=pres[::-1]
            elif (ndim==2):
                pres=pres[:,::-1]
            elif (ndim==3):
                pres=pres[:,:, ::-1]
            else:
                
                self.do_reshape=True
            
        
        
        if (not is_in_log):
            # if the value is not given in log10 
            
            self.grd_pres=pres
            
            # set log10(pres)

            usd_idx=npy.where(pres<>self.mask_val)
            lg_grd_pres=npy.array(pres)
            lg_grd_pres[usd_idx]=npy.log10(self.grd_pres[usd_idx])
            self.lg_grd_pres=lg_grd_pres
            
        else:
            
            # if pres is in log10 
            
            self.lg_grd_pres=npy.array(pres)
            
            # set pressure
            usd_idx=npy.where(pres<>self.mask_val)
            grd_pres=npy.array(pres)
            grd_pres[usd_idx]=10.0**(pres[usd_idx])
            self.grd_pres=grd_pres
            
        
        
        
        
        # coefficients for profile from model grid to obs grid 
        
            
        self.grd_vpl=None
        self.grd_vpr=None
        self.grd_vwgt=None
        
        # coefficients for column values
        
        if (ndim==1):
            self.grd_colwgt=vcol.get_col_wgt_0d(self.lg_grd_pres,  \
                                                    mask_val=self.mask_val)
        elif (ndim==2):
            self.grd_colwgt=vcol.get_col_wgt_1d(self.lg_grd_pres,  \
                                                    mask_val=self.mask_val)
        elif (ndim==3):
            self.grd_colwgt=vcol.get_col_wgt_2d(self.lg_grd_pres,  \
                                                    mask_val=self.mask_val)
        else:
            # high-dimension
            
          
            # reshape 
            nx=1
            for nd in self.grd_dims[0:self.grd_ndim-1]:
                nx=nx*nd
                
            self.reshape_nx=nx
            
            tmp_pres=npy.reshape(self.lg_grd_pres, [nx, self.grd_dims[-1]])
            if (self.do_reverse):
                tmp_pres=tmp_pres[:, ::-1]
                
            self.grd_colwgt=vcol.get_col_wgt_1d(tmp_pres,  \
                                                    mask_val=self.mask_val)
            
        
        
        # coefficients for profile from obs grid to model grid 
        
        self.opres=None
        self.lg_opres=None
        
        self.ob_vpl=None
        self.ob_vpr=None
        self.ob_vwgt=None
        
        self.do_ob_reverse=None
        self.ob_dims=None
        self.ob_ndim=0
        
        self.is_ready=False
        
    def set_mask(self, mask):
        """
        set values for data mask
        
        Inputs:
        ------------------------------------------
        1. mask:<float>: mask for bad or missing values
        
        """
        
        self.mask_val=mask
        self.attr_dict.update({'mask_val':mask})
    

    def set_attr(self, name, val):
        
        """
        add or replace attributes
        
        Inputs:
        ---------------------------------
        1. name:<str>: attribute name
        2. val: <obj>: value of the attribute
        
        
        """
        if (name=='mask_val'):
            self.set_mask(val)
        else:
            self.attr_dict.update({name:val})
    
    def get_attr(self, name):
        
        """
        check attribute
        Inputs:
        ---------------------------------------------
        1. name:<str>: attribute name
        
        Returns:
        ----------------------------------------------
        1. val:<obj>: value of the attribute
        
        """
        
        if (name in self.attr_dict):
            return self.attr_dict[name]
        else:
            msm.show_err_msg(name, msm.msm_no_attr)
            return None

        
    def init_interp(self, opres, is_in_log=False, do_ob_reverse=False):
        """  initialize interpolation coefficients 
        Inputs:
        
        1. opres:<array>: ob pressure grid to be projected. 
        2. is_in_log:<T/F>: optional. Ture if opres is given in log10 space
        3. do_ob_reverse:<T/F>: optional. True of opres is given at descending order 
        
        """
        
        # S1 size check 
        
        dims=npy.shape(opres)
        ndim=npy.size(dims)
        if (ndim<>self.grd_ndim): 
            print 'ob dims:', dims
            print 'grd_dims:', self.grd_dims
            
            msg='Number of  dimensions of the observation grid is different from model grid'
            msm.show_err_msg(msg,  msg_title=msm.msm_wrong_dim)
            return False
        
        if (not oob.check_ot_shape(opres, self.grd_dims[0:-1])):
            
            print 'ob dims:', dims
            print 'grd_dims:', self.grd_dims
            
            msg='Shape of opres is not expected'
            msm.show_err_msg(msg,  msg_title=msm.msm_wrong_dim)
            return False
        
        
        
        
        self.do_ob_reverse=do_ob_reverse
        self.ob_dims=dims
        ob_ndim=npy.size(dims)
        
        self.ob_ndim=ob_ndim
        
        
        # S2 reverse and calculate log10 (pressure) when  necessary 
        
        if (self.do_ob_reverse):
            if (self.ob_ndim==1):
                opres=opres[::-1]
            elif (self.ob_ndim==2):
                opres=opres[:, ::-1]
            elif (self.ob_ndim==3):
                opres=opres[:, :,::-1]
            else:
                # handle later when fortran functions are to be called. 
                pass

        
        
        
                
        if (not is_in_log):
            
            self.opres=npy.array(opres)
            usd_idx=npy.where(opres<>self.mask_val)
            
            # to log10(pressure)
            
            self.lg_opres=npy.array(opres)
            self.lg_opres[usd_idx]=npy.log10(opres[usd_idx])
        else:
            
            self.lg_opres=npy.array(opres)
            usd_idx=npy.where(opres<>self.mask_val)
            
            # to pressure 
            
            self.opres=npy.array(opres)
            self.opres[usd_idx]=10**(opres[usd_idx])
            
        # S3 calculate interpolation coefficients
        
        
        if (self.ob_ndim==1):  # opres in shape of (nz)
            
            # coefficents: obs profile ---> model profile 
            self.grd_vpl, self.grd_vpr, self.grd_vwgt=\
                vpf.get_vertical_wgt_0d(self.lg_opres, self.lg_grd_pres, mask_val=self.mask_val)
            
            # coefficents: model profile ---> obs profile
            self.ob_vpl, self.ob_vpr, self.ob_vwgt=\
                vpf.get_vertical_wgt_0d(self.lg_grd_pres, self.lg_opres,  mask_val=self.mask_val)
        
            # columns
            self.ob_colwgt=vcol.get_col_wgt_0d(self.lg_opres, mask_val=self.mask_val)
            

        elif (self.ob_ndim==2): # shape of (nx, nz)
            
            # coefficents from obs space to model space
            self.grd_vpl, self.grd_vpr, self.grd_vwgt=\
                vpf.get_vertical_wgt_1d(self.lg_opres, self.lg_grd_pres, mask_val=self.mask_val)
            
            # coefficents from model space to obs space
            
            self.ob_vpl, self.ob_vpr, self.ob_vwgt=\
                vpf.get_vertical_wgt_1d(self.lg_grd_pres, self.lg_opres, mask_val=self.mask_val)
            
            # columns
            self.ob_colwgt=vcol.get_col_wgt_1d(self.lg_opres, mask_val=self.mask_val)
            
        elif (self.ob_ndim==3):  # shape of (nx, ny,nz)
            
            
            # coefficents: model profile --> obs profile 
            self.grd_vpl, self.grd_vpr, self.grd_vwgt=\
                vpf.get_vertical_wgt_2d(self.lg_opres, self.lg_grd_pres, mask_val=self.mask_val)
            
            # coefficents: obs profiles ---> model profiles 
            self.ob_vpl, self.ob_vpr, self.ob_vwgt=\
                vpf.get_vertical_wgt_2d(self.lg_grd_pres, self.lg_opres, mask_val=self.mask_val)
            
            # columns
            
            self.ob_colwgt=vcol.get_col_wgt_2d(self.lg_opres, mask_val=mask_val)
        
        else:  # high dimension
            
            tmp_lg_opres=npy.reshape(self.lg_opres, [self.reshape_nx, self.ob_dims[-1]])
            if (self.do_ob_reverse):
                tmp_lg_opres=tmp_lg_opres[:,::-1]
                
            
            
            tmp_lg_pres=npy.reshape(self.lg_grd_pres, [self.reshape_nx, self.grd_dims[-1]])
            
            if (self.do_reverse):
                tmp_lg_pres=tmp_lg_pres[:,::-1]
            
            self.grd_vpl, self.grd_vpr, self.grd_vwgt=\
                vpf.get_vertical_wgt_1d(tmp_lg_opres, \
                                            tmp_lg_pres, mask_val=self.mask_val)
            
            self.ob_vpl, self.ob_vpr, self.ob_vwgt=\
                vpf.get_vertical_wgt_1d(tmp_lg_pres, tmp_lg_opres, mask_val=self.mask_val)

            del tmp_lg_opres
            del  tmp_lg_pres
            
            

        
        self.is_ready=True
    
    def interpolate_mod_prof(self, fld):
        
        """
        Interpolate model profiles to the ob (target) pressure grid defined at init_interp
        
        Inputs:
        --------------------------------------------------------------------------------
        1. fld:<array, ([nx],[ny], nz, [ne]>: 1D, 2D or 3D gridded data at model grid.  Its shape is \
        assumed to be (nz, [ne]) 1D;   (nx, nz, [ne]), (2D); or (nx, ny, nz, [ne]),  (3D).
        
        ---[nx] * [ny]: Number of horizontal locations; 
        ---nz:  Size of the model vertical grid. 
        ---[ne]: Tracer (time step) number 
        
        
        Outputs:
        -----------------------------------------------------------------------------------
        1. prof_at_ob:<array,([nx],[ny], nz_ob, [ne]) >: gp_fld at the new (ob) vertical grid 
        of nz_ob levels

        
        """

        
        # S1:  check the data size
        
        dims=npy.shape(fld)
        ndim=npy.size(dims)
        
        
        # S2: check status and shape
        
        
        if (not self.is_ready):
            msg='interpolation has not been initiailized'
            msm.show_err_msg(msg)
            return None
        
        if (not oob.check_ot_shape(fld, self.grd_dims)):
            
            msg='Data have a shape different from obs pressure'
            msm.show_err_msg(msg)
            return None
        
        
        # S3: reverse if necessary 
        
        
        
        if (self.do_reverse):
            if (ndim==self.grd_ndim): # single tracer
                # single tracer
                
                if (ndim==1):
                    data=fld[::-1]
                
                elif (ndim==2):
                    data=fld[:, ::-1]
                
                elif (ndim==3):
                
                    data=fld[:, :, ::-1]
                else:
                    # reshape and reverse
                    data=npy.reshape(fld, [self.reshape_nx, dims[-1]])
                    data=data[:,::-1]
            
                    
                    
                
                    
            
            elif (ndim==self.grd_ndim+1):  
                # multiple tracers (or time)
                
                if (ndim==2):
                    data=fld[::-1,:]
                elif (ndim==3):
                    data=fld[:, ::-1,:]
                elif (ndim==4):
                    data=fld[:, :, ::-1,:]
                else:
                    # reshape and reverse
                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                    data=data[:,::-1,:]
                    
                    
                    
                
        else:
            
            data=fld
            
            if (self.do_reshape):
                if (ndim==self.grd_ndim):
                    data=npy.reshape(fld, [self.reshape_nx, dims[-1]])
                elif (ndim==self.grd_ndim+1):
                    # reshape as [nx, nz, ne]
                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                
        
        # S4: Interpoplation
        
            
        if (ndim==self.grd_ndim): 
            ## B1 single tracer 
            if (ndim==1):
                prof_at_ob=vpf.prof_vertical_intpl_0d(self.ob_vpl, \
                                                          self.ob_vpr, self.ob_vwgt, data,mask_val=self.mask_val)
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[::-1]
                    
            elif (ndim==2):
                prof_at_ob=vpf.prof_vertical_intpl_1d(self.ob_vpl, \
                                                          self.ob_vpr, self.ob_vwgt, data,mask_val=self.mask_val)
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:, ::-1]
            
            elif (ndim==3):
                prof_at_ob=vpf.prof_vertical_intpl_2d(self.ob_vpl, \
                                                          self.ob_vpr, self.ob_vwgt, data, mask_val=self.mask_val)

                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:,:, ::-1]
            else:
                prof_at_ob=vpf.prof_vertical_intpl_1d(self.ob_vpl, \
                                                          self.ob_vpr, self.ob_vwgt, data,mask_val=self.mask_val)
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:, ::-1]
                
                # reshape
                new_dims=self.ob_dims
                # reshape to same as shape of  obs pressure
                prof_at_ob=npy.reshape(prof_at_obs, new_dims)
                
                
                
            
        elif (ndim==self.grd_ndim+1):  
            ## B2:  multiple tracer (or time)
            
            ## T1:
            
            if (ndim==2):
                prof_at_ob=vpf.prof_vertical_intpl_0d_em(self.ob_vpl, \
                                                             self.ob_vpr, self.ob_vwgt, data, mask_val=self.mask_val)
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[::-1,:]
                    
                
            elif (ndim==3):
                prof_at_ob=vpf.prof_vertical_intpl_1d_em(self.ob_vpl, \
                                                             self.ob_vpr, self.ob_vwgt, data, mask_val=self.mask_val)
            
                
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:, ::-1,:]

            elif (ndim_grd==3):
                prof_at_ob=vpf.prof_vertical_intpl_2d_em(self.ob_vpl, \
                                                             self.ob_vpr, self.ob_vwgt, data, mask_val=self.mask_val)
                
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:, :, ::-1,:]

            else:
                prof_at_ob=vpf.prof_vertical_intpl_1d_em(self.ob_vpl, \
                                                             self.ob_vpr, self.ob_vwgt, data, mask_val=self.mask_val)
                
                if (self.do_ob_reverse):
                    prof_at_ob=prof_at_ob[:, ::-1,:]

                # reshape 
                new_dims=self.ob_dims
                new_dims=new_dims+tuple([-1])
                prof_at_ob=npy.reshape(prof_at_obs, new_dims)

        else:
            # too many dimensions 
            
            msg='Data have too many more dimensions than pressure field'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
            
                
        return prof_at_ob


    
    def interpolate_ob_prof(self, ob_fld, mask_val=oob.fill_val):
        

        """
        
        Interpolate ob profiles to model pressure grid
        
        Inputs:
        --------------------------------------------------------------------------------
        1. ob_fld:<array, ([nx],[ny], nz0>: 2D or 3D gridded data at ob grid.  Its reshape is 
        assumed to be (nz) 1D;   (nx, nz), (2D); or (nx, ny, nz),  (3D).
        
        ---[nx] * [ny]: Number of horizontal locations; 
        ---nz0:  Size of the vertical grid. 
        
        2. mask_val:<float>: missing or bad value

        Outputs:
        -----------------------------------------------------------------------------------
        1. prof_at_ob:<array,([nx],[ny], nz>: gp_fld at the new (model) vertical grid 
        
        """
        
        
        
        
        # S1: check the data size
        dims=npy.shape(ob_fld)
        ndim=npy.size(dims)
        
        ndim_ob=self.ob_ndim
        ndim_grd=self.grd_ndim
        
        
        # S2:  check 
        ## T1 check if it is ready 
        if (not self.is_ready):
            
            msg='interpolation has not been initiailized'
            msm.show_err_msg(msg)
            return None
        
        ## T2 check the dimension 
        
        if (not oob.check_ot_shape(ob_fld, self.ob_dims)):
            msg='Data have a shape different from the stored ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
            
         
        # S3: reverse input profiles if necessary
        
                
        if (self.do_ob_reverse):
            if (ndim==1):
                data=ob_fld[::-1]
            elif (ndim==2):
                data=ob_fld[:, ::-1]
            elif (ndim==3):
                data=ob_fld[:, :, ::-1]
            else:
                data=reshape(ob_fld, [self.reshape_nx, dims[-1]])
                data=data[:,::-1]
                             
             
        else:
            # no reverse
            
            data=ob_fld
            if (self.do_reshape):
                data=reshape(ob_fld, [self.reshape_nx, dims[-1]])
                
        
        # S4: Interpoplation
                
        
        if (ndim==1):
            prof_at_grd=\
                vpf.prof_vertical_intpl_0d(self.grd_vpl, self.grd_vpr, self.grd_vwgt, \
                                               data, mask_val=self.mask_val)
            if (self.do_reverse):
                prof_at_grd=prof_at_grd[::-1]
       
        elif (ndim==2):
                prof_at_grd=\
                    vpf.prof_vertical_intpl_1d(self.grd_vpl, self.grd_vpr, self.grd_vwgt, \
                                                   data, mask_val=self.mask_val)
                if (self.do_reverse):
                    prof_at_grd=prof_at_grd[:, ::-1]
        
        elif (dim==3):
            prof_at_grd=\
                vpf.prof_vertical_intpl_2d(self.grd_vpl, self.grd_vpr, self.grd_vwgt, \
                                               data, mask_val=self.mask_val)
            if (self.do_reverse):
                prof_at_grd=prof_at_grd[:, :, ::-1]
        
        else:
            prof_at_grd=\
                vpf.prof_vertical_intpl_1d(self.grd_vpl, self.grd_vpr, self.grd_vwgt, \
                                               data, mask_val=self.mask_val)
        
            if (self.do_reverse):
                prof_at_grd=prof_at_grd[:, :, ::-1]
                
            prof_at_grd=npy.reshape(prof_at_grd, self.grd_dims)
        
            
        
        return prof_at_grd
    
    
    def get_mod_column(self, gp_fld):
        
        """
        
        get column values for profiles at model pressure grid
        
        Inputs:
        --------------------------------------------------------------------------
        1. gp_fld:<array>: 1D, 2D or 3D gridded data at model (or slice) grid.  Its shape is 
        assumed to be ([nx], [ny], nz, [ne])
        
        ---[nx] * [ny]: Number of horizontal locations; 
        ---nz:  Size of the model vertical grid. 
        ---[ne]: Tracer (time step) number 
        
        2. mask_val:<float>: missing or bad value

        Outputs:
        
        1. col_at_grd <array,([nx], [ny], [ne]) >: column values 
        
        """
        # S1: check dimensions
        
        dims=npy.shape(gp_fld)
        ndim=npy.size(dims)
        
        # grd (model\slot) pressure dimension 
        
        ndim_grd=self.grd_ndim
        
        if (not oob.check_ot_shape(gp_fld, self.grd_dims)):
            msg='Data have a shape different from the stored ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
        
        # S2: reverse and reshape if necessary 
        
        if (self.do_reverse):
            
            if (ndim==self.grd_ndim):
                
                # single tracer
                
                if (ndim==1):
                    data=gp_fld[::-1]
                elif (ndim==2):
                    data=gp_fld[:, ::-1]
                elif (ndim==3):
                    data=gp_fld[:, :, ::-1]
                else:
                    data=npy.reshape(gp_fld, [self.nx,dims[-1]])
                    data=data[:,::-1]
                    
            else:
                # multiple tracers
                
                if (ndim==2):
                    data=gp_fld[::-1,:]
                elif (ndim==3):
                    data=gp_fld[:, ::-1,:]
                elif (ndim==4):
                    data=gp_fld[:, :, ::-1,:]
                else:
                    data=npy.reshape(gp_fld, [self.nx,dims[-2], dims[-1]])
                    data=data[:,::-1,:]
                    
        else:
            
            data=gp_fld
            
            if (self.do_reshape):
                if (ndim==self.grd_ndim):
                    data=npy.reshape(gp_fld, [self.nx, -1])
                elif (ndim==self.grd_ndim+1):
                    data=npy.reshape(gp_fld, [self.nx, dims[-2], dims[-1]])
                else:
                    msg='Data have too many more dimensions than pressure field'
                    msm.show_err_msg(msg, msm.msm_wrong_dim)
                    return None
        
        # S3 integration
                
        
        if (ndim==ndim_grd):         
            
            if (ndim_grd==1):
                col_gp=vcol.col_int_0d(data, self.grd_colwgt, mask_val=self.mask_val)
                
            elif (ndim_grd==2):
                col_gp=vcol.col_int_1d(data, self.grd_colwgt, mask_val=self.mask_val)
            
            elif (ndim_grd==3):
                col_gp=vcol.col_int_2d(data, self.grd_colwgt, mask_val=self.mask_val)
            
            else:
                col_gp=vcol.col_int_1d(data, self.grd_colwgt, mask_val=self.mask_val)
                
            # reshape to grid 

            if (self.do_reshape):
                col_gp=npy.reshape(col_gp, self.grd_dims[0:-1])
                

        elif (ndim==ndim_grd+1):
            
            if (ndim==2):
                col_gp=vcol.col_int_0d_em(data, self.grd_colwgt,mask_val=self.mask_val)
            elif (ndim==3):
                col_gp=vcol.col_int_1d_em(data, self.grd_colwgt, mask_val=self.mask_val)
            elif (ndim==4):
                col_gp=vcol.col_int_2d_em(data, self.grd_colwgt, mask_val=self.mask_val)
            else:
                col_gp=vcol.col_int_1d_em(data, self.grd_colwgt, mask_val=self.mask_val)

            # reshape if necessary
            
            if (self.do_reshape):
                new_dims=self.grd_dims[0:-1]+tuple([-1])
                col_gp=npy.reshape(col_gp, new_dims)
            
                
        
                

        return col_gp

    
    def get_ob_column(self, fld):
        
        """
        
        get column values at the observation pressure grid defined at init_interp

        Inputs:
        1. fld:<array, ([nx], [ny], nz, [ne])>:1D, 2D or 3D data at ob grid.
        
        
        
        Outputs:
        1. col_ob <array>: column values 
        
        """
        
        dims=npy.shape(fld)
        ndim=npy.size(dims)
        
        ndim_ob=self.ob_ndim
        
        # S1 check initialisation of obs grid 
        
        if (not self.is_ready):
            msg='interpolation has not been initiailized'
            msm.show_err_msg(msg)
            return None
        
        # S2 check Size

        if (not oob.check_ot_shape(fld, self.ob_dims)):
            msg='Data have a shape different from the stored ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
        # S3 reverse and reshape 
        
        if (self.do_ob_reverse):
            if (ndim==self.ob_ndim):
                
                if (ndim==1):
                    data=fld[::-1]
                elif (ndim==2):
                    data=fld[:, ::-1]
                elif (ndim==3):
                    data=fld[:, :, ::-1]
                else:
                    data=reshape(fld, [self.reshape_nx, -1])
                    data=data[:,::-1]
            
            elif (ndim==self.ob_ndim+1):
                # multiple tracers

                if (ndim==2):
                    data=fld[::-1,:]
                elif (ndim==3):
                    data=fld[:, ::-1,:]
                elif (ndim==4):
                    data=fld[:, :, ::-1,:]
                else:
                    # high-demnsion data

                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                    data=data[:,::-1,:]
            else:
                msg='Data have too many more dims than ob pressure'
                msm.show_err_msg(msg, msm.msm_wrong_dim)
                return None
            
        else:

            data=fld
            if (self.do_reshape):
                if (ndim==self.ob_ndim):
                    data=reshape(fld, [self.reshape_nx, -1])
                
                elif (ndim==self.ob_ndim+1):
                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                
                else:
                     msg='Data have too many more dims than ob pressure'
                     msm.show_err_msg(msg, msm.msm_wrong_dim)
                     return None

        
        # S3: col integration 
        
        if (ndim==self.ob_ndim):
        
            if (ndim==1):       
                col_ob=vcol.col_int_0d(data, self.ob_colwgt, mask_val=self.mask_val)
            
            elif (ndim==2):
                col_ob=vcol.col_int_1d(data, self.ob_colwgt, mask_val=self.mask_val)
            
            elif (ndim==3):
                col_ob=vcol.col_int_2d(data, self.ob_colwgt, mask_val=self.mask_val)

            else:
                
                # high-dimension data
                col_ob=vcol.col_int_1d(data, self.ob_colwgt, mask_val=self.mask_val)
                
                # reshape 
                col_ob=npy.reshape(col_ob, self.ob_dims[0:-1])

                
        elif (ndim==self.ob_ndim+1):
            # mult-tracer 
            
            if (ndim==2):       
                col_ob=vcol.col_int_0d_em(data, self.ob_colwgt, mask_val=self.mask_val)

            elif (ndim==3):
                col_ob=vcol.col_int_1d_em(data, self.ob_colwgt, mask_val=self.mask_val)
            
            elif (ndim==4):
                col_ob=vcol.col_int_2d_em(data, self.ob_colwgt, mask_val=self.mask_val)
            
            else:
                col_ob=vcol.col_int_1d_em(data, self.ob_colwgt, mask_val=self.mask_val)
                
                # reshape
                new_dims=self.ob_dims[0:-1]+tuple([-1])
                col_ob=npy.reshape(col_ob, new_dims)
        else:
            msg='Data have too many more dims than ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
            
        return col_ob


    
    def get_ak_ob_column(self, fld, ak):
        
        """
        
        Get column values weigted by averging kernels at the observation pressure grid defined at init_interp
        
        Inputs:
        ===============================================================
        1. fld:<array, ([nx], [ny], nz,[ne])>:1D, 2D or 3D data at ob grid.
        2. ak:<array, ([nx], [ny], nz)>: averaging kernels 
        
        
        Outputs:
        ==================================================================
        1. col_ob <array>: column values 
        
        """
        
        dims=npy.shape(fld)
        ndim=npy.size(dims)
        
        ndim_ob=self.ob_ndim
        
        # S1: check initialisation of obs grid 
        
        if (not self.is_ready):
            msg='interpolation has not been initiailized'
            msm.show_err_msg(msg)
            return None
        
        # S2: check Size
        
        if (not oob.check_ot_shape(fld, self.ob_dims)):
            msg='Data have a shape different from the stored ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
         
        
        if (self.do_ob_reverse):
            # do reverse 
            
            if (ndim==self.ob_ndim):
                
                if (ndim==1):
                    data=fld[::-1]
                    cak=ak[::-1]
                
                elif (ndim==2):
                    data=fld[:, ::-1]
                    cak=ak[:, ::-1]
                
                elif (ndim==3):
                    data=fld[:, :, ::-1]
                    cak=ak[:, :, ::-1]
                
                else:
                    data=npy.reshape(fld, [self.reshape_nx, -1])
                    data=data[:,::-1]
                    
                    cak=npy.reshape(ak, [self.reshape_nx, -1])
                    cak=cak[:,::-1]
                    
            
            elif (ndim==self.ob_ndim+1): # 
                # ensemble 
                print 'do ensemble -ndim-ob_ndim:',  ndim, self.ob_ndim
                
                
                
                if (ndim==2):
                    data=fld[::-1,:]
                    cak=ak[::-1]

                elif (ndim==3):
                    data=fld[:, ::-1,:]
                    cak=ak[:, ::-1]

                elif (ndim==4):
                    data=fld[:, :, ::-1,:]
                    cak=ak[:, :, ::-1]
                    
                else:
                    # high-dimension data 
                    
                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                    data=data[:,::-1,:]
                    
                    cak=npy.reshape(ak, [self.reshape_nx, -1])
                    cak=cak[:,::-1]
            else:
                
                msg='Data have too many more dims than ob pressure'
                msm.show_err_msg(msg, msm.msm_wrong_dim)
                return None
            
        else: # no reverse
            data=fld
            cak=ak
            
            if (self.do_reshape):
                if (ndim==self.ob_ndim):
                    data=reshape(fld, [self.reshape_nx, -1])
                    cak=npy.reshape(ak, [self.reshape_nx, -1])
                    
                elif (ndim==self.ob_ndim+1):
                    print 'ensemble reshape:', [self.reshape_nx, dims[-2], dims[-1]]
                    print 'original shape:',   npy.shape(fld)
                    
                    data=npy.reshape(fld, [self.reshape_nx, dims[-2], dims[-1]])
                    cak=npy.reshape(ak, [self.reshape_nx, -1])
                else:
                    msg='Data have too many more dims than ob pressure'
                    msm.show_err_msg(msg, msm.msm_wrong_dim)
                    return None

        
        # S3 col integration 
        
        if (ndim==self.ob_ndim):
            if (ndim==1):       
                col_ob=vcol.ak_col_int_0d(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            elif (ndim==2):
                col_ob=vcol.ak_col_int_1d(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            elif (ndim==3):
                col_ob=vcol.ak_col_int_2d(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            else:
                col_ob=vcol.ak_col_int_1d(data, self.ob_colwgt, cak, mask_val=self.mask_val)
                # reshape 
                col_ob=npy.reshape(col_ob, self.ob_dims[0:-1])
        
        elif (ndim==self.ob_ndim+1):
            
            if (ndim==2):       
                col_ob=vcol.ak_col_int_0d_em(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            
            elif (ndim==3):
                print 'ak_col_int data:', data[0,:,0]
                
                col_ob=vcol.ak_col_int_1d_em(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            
            elif (ndim==4):
                col_ob=vcol.ak_col_int_2d_em(data, self.ob_colwgt, cak, mask_val=self.mask_val)
            
            else:
                # high-demsion data
                
                col_ob=vcol.col_int_1d_em(data, self.ob_colwgt, cak, mask_val=self.mask_val)
                # reshape
                new_dims=self.ob_dims[0:-1]+tuple([-1])
                col_ob=npy.reshape(col_ob, new_dims)
        else:
            msg='Data have too many more dims than ob pressure'
            msm.show_err_msg(msg, msm.msm_wrong_dim)
            return None
        
        return col_ob


    

        
#<<< TEST >>> 

if (__name__=='__main__'):
    
    print 'for test, see ESA/observation/satellite_operator.py' 
    
