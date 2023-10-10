subroutine get_closest_point(nx0,x0, nx1, x1, pos)
  
  ! find the position in axis x0 closest to points of axis x1 
  !
  ! Input:
  ! 
  ! integer::nx0  # size of original array x0  
  ! integer::nx0  # size of target array x1
  ! real*8::x0(nx0)  # original array (in ascending order)
  ! real*8::x1(nx1)  # target array 
  !
  ! Outputs
  ! integer::pos(nx1)  # closet points for x1(nx1)
  ! 
  
  
  implicit none
  ! inputs
  integer, intent(in)::nx0, nx1
  real*8, intent(in)::x0(nx0), x1(nx1)
  ! outputs
  integer, intent(out)::pos(nx1)
  
  ! locals
  integer::i, j, pst
  real*8::val
  
  

  do i=1, nx1 ! over each point in x1

     val=x1(i)
     pos(i)=-999
     
     do j=1, nx0  ! x0
        
        if (x0(j)>=val) then
           if (j==1) then
              pos(i)=j
           else
              pos(i)=j
              if ((x0(j)-val)>(val-x0(j-1))) pos(i)=j-1
           end if
           exit
           
        end if  ! if
        
        
     end do ! x0
     
     if (pos(i)<0) then 
        pos(i)=nx0
     end if
     
  end do ! x1

  
end subroutine get_closest_point


subroutine get_closest_point_cycle(nx0,x0, nx1, x1, period, pos)
  
  ! find the position in axis x0 closest to point of axis x1 
  ! here x0 and x1 are assumed to be in periodic cycle
  !
  ! Inputs:
  ! ------------------------------------------------------
  ! integer::nx0  # size of original array x0  
  ! integer::nx0  # size of target array x1
  ! real*8::x0(nx0)  # original array (in ascending order)
  ! real*8::x1(nx1)  # target array 
  ! real*8::period   # periodic value
  
  !
  ! Outputs
  !---------------------------------------------------------
  ! integer::pos(nx1)  # closet points for x1(nx1)
  ! 

  
  
  implicit none
  ! inputs
  integer, intent(in)::nx0, nx1
  real*8, intent(in)::x0(nx0), x1(nx1)
  real*8, intent(in)::period  ! period such 360 for longitude 
  

  ! outputs
  integer, intent(out)::pos(nx1)
  
  
  ! locals
  integer::i, j, pst
  real*8::val, shift_val
  
  

  do i=1, nx1 ! over each point in x1

     val=x1(i)
     pos(i)=-999
     
     
     do j=1, nx0  ! x0
        
        if (x0(j)>=val) then
           if (j==1) then
              pos(i)=j
              shift_val=x0(nx0)-period
              if (abs(val-shift_val)<(x0(1)-val)) pos(i)=nx0
           else
              pos(i)=j
              if ((x0(j)-val)>(val-x0(j-1))) pos(i)=j-1
           end if
           exit
           
        end if  ! if
        
     
     end do ! x0
     
     if (pos(i)<0) then 
        pos(i)=nx0
        shift_val=x0(1)+period
        ! 
        if (abs(val-shift_val)<(val-x0(nx0))) pos(i)=1
        
     end if
     
  end do ! x1

  

end subroutine get_closest_point_cycle



subroutine get_segment_wgt(x, nx, y, ny, pos, wgt)
  
  ! get the postion and weighting factor for 
  ! small ny-1 segments defined by y(ny) in 
  !  a larger nx-1 segments defined by x(nx)
  
  ! Inputs
  !----------------------------------------------------------------
  ! integer::nx  # size of x 
  ! integer::ny  # size of y
  ! real*8 ::x(nx) # original nodes 
  ! real*8 ::y(ny) # new nodes
  
  ! outputs
  !-------------------------------------------------------------
  ! integer, intent(out)::pos(ny) ! postion of segment y in segment x
  ! integer, intent(out)::wgt(ny) ! portion (weight) of segment y in segement x
  
  !  
  ! 
  implicit none
  
  integer, intent(in)::nx
  integer, intent(in)::ny
  real*8, intent(in)::x(nx)
  real*8, intent(in)::y(ny)
  
  integer, intent(out)::pos(ny) ! postion of segments defined by  y in segments defined by x
  real*8, intent(out)::wgt(ny) ! portion of segment defined by y located in segement defined by x
  
  ! location variable 
  integer::ix, iy, pst
  
  real*8::rx
  
    
  !c ------ix-------ix+1---------ix+2
  !c  iy-1---iy---y+1
  ! 
  !c pos(iy) denote the left of the first one in x larger than 
  !c y(iy). In the above case pos(iy)=ix
  !  

 ! assign pos and wgt value for each values in y 
  pst=1
  wgt=0.0
  
  do iy=1, ny-1
     pos(iy)=1
     wgt(iy)=0.0
     
     do ix=pst, nx-1
        rx=x(ix+1) ! upper x boundary for segement ix 
      
        
        if (rx>=y(iy+1)) then ! y segment iy at left side of x segment ix+1
           pos(iy)=ix
           if (y(iy+1)<=y(iy)) then
              wgt(iy)=0.0
           else
              wgt(iy)=(y(iy+1)-x(ix))/(y(iy+1)-y(iy))
           end if
           
           if (wgt(iy)>1) then 
              wgt(iy)=1 ! totally inside x segment ix 
           else if (wgt(iy)<0) then
              wgt(iy)=0.0
              
           end if
           
           pst=1
           exit
        end if
     end do ! ix
     
  end do !  iy
  
  ! pos(ny) is in no use 
  pos(ny)=nx-1
  wgt(ny)=1.0
  
end subroutine get_segment_wgt


subroutine regrid_flux_bs(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, flux1, flux2)
  
  !  regrid the flux from larger grid boxes to smaller ones 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  Inputs:
  ! ----------------------------------------------
  !  real*8:: lon1 & lat1 # array of nlon or nlat   edges of original grid box 
  !  real*8:: lon2 & lat2 #  array  of nlon2, nlat2, edges of target grid box
  !  real*8:: flux1: array of (nlon-1, nlat-1)  original surface flux,which is 
  !  assumed to be given at the centers of the grid boxes
  !  
  !  Outputs: 
  !  ----------------------------------------------
  !  real*8:: flux2 #  outarray of (nlon2, nlat2) surface flux at target  grid box   
  
  implicit none
  ! inputs
  ! ! original box edge and values
  integer, intent(in)::nlon, nlat
  real*8, intent(in)::lon1(nlon), lat1(nlat)
  real*8, intent(in)::flux1(nlon-1, nlat-1)
  
  ! target box edges 
  
  integer, intent(in)::nlon2, nlat2
  real*8, intent(in)::lon2(nlon2), lat2(nlat2)

  ! outputs 
  
  real*8, intent(out)::flux2(nlon2-1, nlat2-1)

  ! local variables

  integer, allocatable:: plon(:), plat(:)
  
  real*8 :: flux_11, flux_12, flux_21, flux_22
  
  real*8::fll, ful, flu, fuu, hflon, hflat



  real*8, allocatable::wlon(:)
  real*8, allocatable::wlat(:)
  
  
  integer::llon, ulon, llat, ulat, ix, iy
  
  ! the left ledge 

  allocate(plon(nlon2))
  allocate(wlon(nlon2))
  
  allocate(plat(nlat2))
  allocate(wlat(nlat2))
  

  
  flux2=0.0
  
  ! get segment weight for longitude
  call get_segment_wgt(lon1, nlon, lon2, nlon2, plon, wlon)

  ! get segment weight for latitude
  call get_segment_wgt(lat1, nlat, lat2, nlat2, plat, wlat)
  
  do ix=1, nlon2-1
     ulon=plon(ix)  ! right segment   
     llon=ulon-1    ! left segment 
     if (llon<1) llon=1 
     
     do iy=1, nlat2-1
        ulat=plat(iy) ! up segment 
        llat=ulat-1   ! down segment 
        
        if (llat<1) llat=1
        
        fll=(1.0-wlon(ix))*(1.0-wlat(iy))*flux1(llon, llat)
        ful=wlon(ix)*(1.0-wlat(iy))*flux1(ulon, llat)
        flu=(1.0-wlon(ix))*wlat(iy)*flux1(llon, ulat)
        fuu=wlon(ix)*wlat(iy)*flux1(ulon, ulat)
        
        flux2(ix, iy)=fll+ful+flu+fuu
        ! weigts (for testing)
        fll=(1.0-wlon(ix))*(1.0-wlat(iy))
        ful=wlon(ix)*(1.0-wlat(iy))
        flu=(1.0-wlon(ix))*wlat(iy)
        fuu=wlon(ix)*wlat(iy)
        
        
   end do ! iy
     
  end do ! ix 
  
  
  deallocate(plon)
  deallocate(wlon)
  
  deallocate(plat)
  deallocate(wlat)
  
end subroutine regrid_flux_bs


subroutine regrid_flux_bs_out(nlon, nlat, flux1, plon, wlon, plat, wlat, &
     nlon2, nlat2, flux2)
  
  !  regrid the flux from larger grid boxes to smaller ones 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  Inputs:
  ! ----------------------------------------------
  !  integer:: nlon, nlat  # dimension of grid boundary for flux1
  !  real*8:: flux1   # array of (nlon-1, nlat-1)  original surface flux,which is
  !     assumed to be given at the centers of the grid boxes
  !  integer::plon(nlon2) #  position of segments defined by lon2 in segments defined by lon1
  !  integer::wlon(nlon2) #  portion  of lon2 segment in segment lon1
  !  integer::plat(nlat2) #  position of segments defined by lat2 in segments defined by lat1
  !  integer::wlat(nlat2) #  portion  of lat2 segment in segment lat1 
  ! 
  !  Outputs: 
  !  ----------------------------------------------
  !  real*8:: flux2: out array of (nlon2-1, nlat2-1) surface flux at target  grid box   
  

  implicit none
  ! inputs
  ! ! fluxes and grid
  integer, intent(in)::nlon, nlat
  real*8, intent(in)::flux1(nlon-1, nlat-1)
  
  ! target segment and weight 
  
  integer, intent(in)::nlon2, nlat2
  integer, intent(in):: plon(nlon2), plat(nlat2)
  real*8,  intent(in):: wlon(nlon2), wlat(nlat2)
  
  
  
  ! outputs 
  
  real*8, intent(out)::flux2(nlon2-1, nlat2-1)

  ! local variables

  
  
  real*8::fll, ful, flu, fuu



  
  integer::llon, ulon, llat, ulat, ix, iy
  
  
  flux2=0.0
  
  
  
  do ix=1, nlon2-1
     ulon=plon(ix)  ! left  
     llon=ulon-1    ! right 
     if (llon<1) llon=1 
     
     do iy=1, nlat2-1
        ulat=plat(iy) ! up
        llat=ulat-1   ! dowm
        if (llat<1) llat=1
        
        fll=(1.0-wlon(ix))*(1.0-wlat(iy))*flux1(llon, llat)
        ful=wlon(ix)*(1.0-wlat(iy))*flux1(ulon, llat)
        flu=(1.0-wlon(ix))*wlat(iy)*flux1(llon, ulat)
        fuu=wlon(ix)*wlat(iy)*flux1(ulon, ulat)
        
        flux2(ix, iy)=fll+ful+flu+fuu
        
        
        if ((flux2(ix,iy)>flux1(llon, llat)).and.(flux2(ix,iy)>flux1(ulon, llat))&  
             .and.(flux2(ix,iy)>flux1(llon, ulat))& 
             .and.(flux2(ix,iy)>flux1(ulon, ulat))) then
           print *, 'warning regrid_flux_bs_out'

           print *, 'max wrong', ix, iy
           print *, 'total', flux2(ix,iy)
           print *, flux1(llon, llat), flux1(ulon, llat), flux1(llon, ulat) ,flux1(ulon, ulat)
           print *,fll, ful, flu, fuu
           
        end if
        
        
        
     end do ! iy
     
  end do ! ix 
  
  
  
end subroutine regrid_flux_bs_out


subroutine regrid_flux_smp(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, influx, outflux)
  
  !  regrid the flux from a smaller grid boxes to larger grid boxes 
  !  
  !  Inputs:    
  !  ----------------------------------------------------------------------
  !  integer ::nlon # size of longitude
  !  integer ::nlat # size of latitude
  !  
  !  real*8 ::lon1(nlon)  # longitude  lat1(nlat)
  !  real*8 ::lat1(nlat)  # latitude
  !  real*8 ::influx(nlon, nlat)  # fluxes values per grid 
  !  integer ::nlon2 # size of target longitude 
  !  integer ::nlat2 # size of target latitude 
  !  real*8  :: lon2 # target longitude 
  !  real*8  :: lat2 # target latitude 
  ! 
  !  outputs:
  !  -------------------------------------------------------
  !  real*8::outflux  # flux in  target grid 


  !  Notes: 
  !  1. we need to check the map to make sure the surface type is consistent with each other 
  
  !  2. all longitude, latitude  and flux  are given at  centre of  grid boxes.  
  
  implicit none
  integer, intent(in)::nlon, nlat
  real*8, intent(in)::lon1(nlon), lat1(nlat)
  real*8, intent(in)::influx(nlon, nlat)


  integer, intent(in)::nlon2, nlat2
  real*8, intent(in)::lon2(nlon2), lat2(nlat2)
  real*8, intent(out)::outflux(nlon2, nlat2)
  
  ! local variables
  integer, allocatable:: plon(:), plat(:)
  
  integer::ilon2, ilat2, ilon, ilat
  
  
  allocate(plon(nlon))
  allocate(plat(nlat))

  


  ! search the lon1
 
  
  outflux=0.0
  
  


  ! search for closest lon2 location for origin lon
  call get_closest_point(nlon2, lon2, nlon, lon1, plon)
  
  ! search for closest lat2 location for origin lat
  
  call get_closest_point(nlat2, lat2,nlat, lat1, plat)
! add the flux from small boxes to larger (output) boxes 

  do ilon=1, nlon  
     ilon2=plon(ilon)
     do ilat=1, nlat
        ilat2=plat(ilat)
        outflux(ilon2, ilat2)=outflux(ilon2, ilat2)+influx(ilon, ilat)
     end do ! ilat
  end do   ! ilon
  

  deallocate(plon)
  deallocate(plat)


end subroutine regrid_flux_smp


subroutine extract_map_to_multi_layer(rlon1, rlat1, nlon, nlat, &
     rlon2, rlat2, nlon2, nlat2, & 
     inmap, ml, outmap)
  
  !  extract map to multiple-layer maps 
  ! Inputs:
  ! -----------------------------------------
  ! integer::nlon, nlat # original longitude and latitude size
  ! real*8::rlon1(nlon), rlat1(nlat)  # original longitude and latitude grid
  ! integer::inmap(nlon, nlat) #  regional map
  ! integer::ml  # ml=max(inmap) is the layer number of the map 
  ! integer::nlon2, nlat2  # output longitude and latitude size 
  ! real*8::rlon2(nlon2), rlat2(nlat2) # output longitude and latitude grid
  
  ! outputs
  ! ----------------------------------------------------
  ! real*8 ::outmap(nlon2, nlat2, ml)
  
  
  integer, intent(in)::nlon, nlat ! longitude and latitude size
  real*8, intent(in)::rlon1(nlon), rlat1(nlat)  ! original longitude and latitude grid
  integer, intent(in)::inmap(nlon, nlat) !  regional map
  integer, intent(in)::ml  ! ml=max(inmap) is the layer number of the map 
  integer, intent(in)::nlon2, nlat2
  real*8, intent(in)::rlon2(nlon2), rlat2(nlat2)
  
  ! outs
  real*8, intent(out)::outmap(nlon2, nlat2, ml)

  
  ! locals
  real*8, allocatable::reg_cnt(:,:,:)
  integer, allocatable:: plon(:), plat(:)
  integer::ilat2, ilon2
  integer::ilat, ilon, iz
  real*8 :: wgt
  
  
  allocate(plon(nlon))
  allocate(plat(nlat))
  allocate(reg_cnt(nlon2, nlat2,ml))
  
  reg_cnt=0.0
  outmap=0.0
   ! search for closest lon2 location for origin lon
  call get_closest_point(nlon2, rlon2, nlon, rlon1, plon)
  
  ! search for closest lat2 location for origin lat
  
  call get_closest_point(nlat2, rlat2, nlat, rlat1, plat)

  
  ! count inmap to layer map 
  
  do ilon=1, nlon
     ilon2=plon(ilon)
     do ilat=1, nlat
        ilat2=plat(ilat)
        iz=inmap(ilon, ilat)
        reg_cnt(ilon2, ilat2, iz)=reg_cnt(ilon2, ilat2, iz)+1
     end do ! ilat
  end do   ! ilon
  
  ! convert counts to ratio at (lon2, lat2)
  
  
  do ilon=1, nlon2
     do ilat=1, nlat2
        wgt=sum(reg_cnt(ilon, ilat, :))
        if (wgt>0) outmap(ilon, ilat,:)=reg_cnt(ilon, ilat, :)/wgt
     end do
end do

deallocate(reg_cnt)
deallocate(plon)
deallocate(plat)

end subroutine extract_map_to_multi_layer


subroutine extract_flux_to_multi_layer(rlon1, rlat1, nlon, nlat, &
     rlon2, rlat2, nlon2, nlat2, & 
     inmap, ml, influx, outmap, outflux)
  
  !  extract map and flux to multiple-layer basis functions 
  !  Inputs:
  ! -----------------------------------------
  ! integer::nlon, nlat # original longitude and latitude size
  ! real*8::rlon1(nlon), rlat1(nlat)  # original longitude and latitude grid
  ! integer::inmap(nlon, nlat) #  regional map
  ! integer::influx(nlon, nlat) #  emissions per grid box (not per area)
  ! integer::ml  # ml=max(inmap) is the layer number of the map 
  ! integer::nlon2, nlat2  # output longitude and latitude size 
  ! real*8::rlon2(nlon2), rlat2(nlat2) # output longitude and latitude grid
  
  ! outputs
  ! ----------------------------------------------------
  ! real*8, intent(out)::outmap(nlon2, nlat2, ml)
  ! real*8, intent(out)::outflux(nlon2, nlat2, ml)
  
    ! --------------------------------------
  
  ! inputs
  integer, intent(in)::nlon, nlat ! longitude and latitude size
  real*8, intent(in)::rlon1(nlon), rlat1(nlat)  ! original longitude and latitude grid
  integer, intent(in)::inmap(nlon, nlat) !  regional map
  integer, intent(in)::ml  ! ml=max(inmap) is the layer number of the map 
  real*8, intent(in)::influx(nlon, nlat) !  emissions from each grid box (not per area)
  
  
  integer, intent(in)::nlon2, nlat2
  real*8, intent(in)::rlon2(nlon2), rlat2(nlat2)
  
  ! outputs
  real*8, intent(out)::outmap(nlon2, nlat2, ml)
  real*8, intent(out)::outflux(nlon2, nlat2, ml)
  
  
  ! locals
  real*8, allocatable::reg_cnt(:,:,:)   ! counts of regions at each target grid   
  integer, allocatable:: plon(:), plat(:) ! longitude and latitude lower boundary
  

  integer::ilat2, ilon2
  integer::ilat, ilon, iz
  real*8 :: wgt
  
  
  allocate(plon(nlon))
  allocate(plat(nlat))
  allocate(reg_cnt(nlon2, nlat2,ml))
  
  reg_cnt=0.0
  outmap=0.0
  outflux=0.0
  
  ! search for closest lon2 location for origin lon
  call get_closest_point(nlon2, rlon2, nlon, rlon1, plon)
  
  ! search for closest lat2 location for origin lat
  
  call get_closest_point(nlat2, rlat2,nlat, rlat1, plat)
  
  
  ! count inmap to layer map 
  ! accumulate flux to outflux layer
  
  do ilon=1, nlon
     ilon2=plon(ilon)
     do ilat=1, nlat
        ilat2=plat(ilat)
        iz=inmap(ilon, ilat)
        outflux(ilon2, ilat2, iz)=outflux(ilon2, ilat2, iz)+influx(ilon, ilat)
        reg_cnt(ilon2, ilat2, iz)=reg_cnt(ilon2, ilat2, iz)+1
     end do ! ilat
  end do   ! ilon
  
  ! convert counts to ratio at (lon2, lat2)
  
  
  do ilon=1, nlon2
     do ilat=1, nlat2
        wgt=sum(reg_cnt(ilon, ilat, :))
        if (wgt>0) outmap(ilon, ilat,:)=reg_cnt(ilon, ilat, :)/wgt
     end do
end do

deallocate(reg_cnt)
deallocate(plon)
deallocate(plat)

end subroutine extract_flux_to_multi_layer



subroutine get_upper_index(nx0,x0, nx1, x1, pos)

  ! find the position of the upper boundary in axis \
  ! x0 for each point of axis x1 
  
  ! Input:
  !--------------------------------------------------- 
  ! integer::nx0  # size of original array x0  
  ! integer::nx0  # size of target array x1
  ! real*8::x0(nx0)  # original array (in ascending order)
  ! real*8::x1(nx1)  # target array 
  !
  ! Outputs
  ! --------------------------------------------------------------
  ! integer::pos(nx1)  # upper boundary for x1(nx1)
  ! 
  
  
  implicit none
  ! inputs
  integer, intent(in)::nx0, nx1
  real*8, intent(in)::x0(nx0), x1(nx1)
  ! outputs
  integer, intent(out)::pos(nx1)
  
  ! locals
  integer::i, j, pst
  real*8::val
  
  do i=1, nx1 !x1
     val=x1(i)
     pos(i)=-999
     
     do j=1, nx0  ! x0
        if (x0(j)>=val) then
           pos(i)=j 
           exit
        end if  ! if 
     end do ! x0 
     
     if (pos(i)<0) then 
        pos(i)=nx0
     end if
     
  end do ! x1



end subroutine get_upper_index





subroutine refine_map(lon1, lat1, nlon, nlat, & 
     lon2, lat2, nlon2, nlat2, inmap, outmap)
  
  !  regrid map for coarse grid to a finer grid 
  !  
  !  Inputs:
  !-----------------------------------------------------------------
  !  integer:: nlon, nlat # size of original longitude and latitude 
  !  real*8:: lon1, lat1  # original longitude and latitude 
  !  integer:: nlon2, nlat2 # size of target longitude and latitude 
  !  real*8,  lon2 & lat2   # target longitudes and latitudes 
  !  real*8:: inmap(nlon, nlat) #  original region map
  !  flux1    in  array of (nlon, nlat)  origin surface flux

  !  Outputs
  !--------------------------------------------------------------------
  !  :: outmap(nlon2, nlat2) # surface map 

  
implicit none
integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
integer, intent(in)::inmap(nlon, nlat)

integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
integer, intent(out)::outmap(nlon2, nlat2)

integer, allocatable:: plon(:), plat(:)
integer ::plon1, plat1, plon2, plat2


! locals

integer::pst, sel_ix, sel_iy, ix, iy



integer::latp1, latp2, lonp1, lonp2

real*8::rlon, rlat



allocate(plon(nlon2))
allocate(plat(nlat2))

! search the closest point 

call get_closest_point(nlon, lon1, nlon2, lon2, plon)
call get_closest_point(nlat, lat1, nlat2, lat2, plat)



do ix=1, nlon2
   sel_ix=plon(ix)
   do iy=1, nlat2
      sel_iy=plat(iy)
      outmap(ix, iy)=inmap(sel_ix, sel_iy)
   end do
end do

deallocate(plon)
deallocate(plat)

end subroutine refine_map



subroutine read_std_map(lon, lat, type_map, nx, ny)
  ! read type mape from binary file smoothmap.fix.2.bin
  
  !  Inputs:
  !-----------------------------------------------------------------
  !  integer:: nx, ny # size of original longitude and latitude 

  !  outputs
  !------------------------------------------------------------
  !  real*4:: lon, lat  # original longitude and latitude 
  !  real*4:: type_map  # map 
  
implicit none
real*4, intent(out)::lon(nx)
real*4, intent(out)::lat(ny)
real*4, intent(out)::type_map(nx, ny)

integer, intent(in)::nx, ny
real*4::dlon, dlat
integer::funit, i, j
real*4 ::map(360, 180)


funit=36
dlon=360.0/nx
lon(1)=-180+0.5*dlon
do i=2, nx
   lon(i)=lon(i-1)+dlon
end do

dlat=180.0/ny
lat(1)=-90+0.5*dlat
do j=2, ny
   lat(j)=lat(j-1)+dlat
end do


print *, nx, ny

open(funit, file='smoothmap.fix.2.bin', form='unformatted')
! print *, 'try to read'
read(funit) ((map(i, j), i=1, 360), j=1,180)
do i=1, nx
   do j=1, ny
      type_map(i,j)=map(i,j)
   end do
end do

end subroutine read_std_map




subroutine do_flux_by_coef(nx, ny, nz, flux, nm, coef, flux_out)

! calculate matrix multiplication: flux*coef
! inputs
!----------------------------------------------
! integer::nx, ny, nz # model grid size 
! real*8 ::flux(nx, ny, nz) # fluxes to be multiplied
! real*8, intent(in)::coef(nz, nm) # coefficent matrix  

! outputs
! -----------------------------------------------
! real*8, intent(out)::flux_out(nx, ny, nm) # 


implicit none 

integer, intent(in)::nx, ny, nz, nm
real*8, intent(in)::flux(nx, ny, nz)
real*8, intent(in)::coef(nz, nm)
real*8, intent(out)::flux_out(nx, ny, nm)

! local 
integer::i, j, k, m
do i=1, nx
   do j=1, ny
      do m=1, nm
         flux_out(i,j,m)=0.0
         do k=1, nz
            flux_out(i,j,m)=flux_out(i,j,m)+ &
                 flux(i, j, k)*coef(k, m)
         end do ! k
         
      end do ! m
   end do  ! j
end do ! i
         


end subroutine do_flux_by_coef


subroutine do_flux_by_coef_wgt(nx, ny, nz, flux, wgt, nm, coef, flux_out)

! calculate matrix multiplication: flux*wgt*coef
! inputs
!----------------------------------------------
! integer::nx, ny, nz # model grid size, nz is the layer  
! real*8 ::flux(nx, ny, nz) # fluxes to be multiplied
! real*8 ::coef(nz, nm) # coefficent matrix
! real*8 ::wgt(nx, ny) # scaling weight     

! outputs
! -----------------------------------------------
! real*8, intent(out)::flux_out(nx, ny, nm) # 


implicit none 

integer, intent(in)::nx, ny, nz, nm
real*8, intent(in)::flux(nx, ny, nz)
real*8, intent(in)::wgt(nx, ny)
real*8, intent(in)::coef(nz, nm)
real*8, intent(out)::flux_out(nx, ny, nm)

! local 
integer::i, j, k, m

flux_out=0.0

do i=1, nx
   do j=1, ny
      do m=1, nm
         flux_out(i,j,m)=0.0
         do k=1, nz
            flux_out(i,j,m)=flux_out(i,j,m)+ &
                 wgt(i,j)*flux(i, j, k)*coef(k, m)
         end do ! k
         
      end do ! m
   end do  ! j
end do ! i
         


end subroutine do_flux_by_coef_wgt




subroutine add_val_to_map(nx, ny, inmap, npt, val, idx, idy, outmap)

! add values for selected points to map location  
! inputs
!----------------------------------------------
! integer::nx, ny, nz # model grid size 
! real*8 ::flux(nx, ny, nz) # fluxes to be multiplied
! real*8, intent(in)::coef(nz, nm) # coefficent matrix  

! outputs
! -----------------------------------------------
! real*8, intent(out)::outmap(nx, ny, nm) # outvalues

implicit None

integer, intent(in)::nx, ny, npt
real*8, intent(in)::inmap(nx, ny)
real*8, intent(in):: val(npt)
integer*8, intent(in):: idx(npt), idy(npt)
! output
real*8, intent(out)::outmap(nx, ny)


! local


integer::ix, jy, ip

outmap(:,:)=inmap(:,:)


do ip=1, npt
   ix=idx(ip)+1
   jy=idy(ip)+1
   outmap(ix, jy)=outmap(ix, jy)+val(ip)
end do

end subroutine add_val_to_map






subroutine get_obs_wgt_map(nlon, nlat,&
     grd_lon, grd_lat, &
     nobs,  &
     olon, olat, &
     dist_table, wlen, lon_period, wgt)
  
  !  Weigting observations according to distance at map grid boxes
  !  
  !  Inputs:
  ! ----------------------------------------------

  ! integer ::nlon # size of x grid
  ! integer ::nlat # size of y grid 
  ! real*8 ::grd_lon(nlon) # longitude
  ! real*8 ::grd_lat(nlat) # latitude
  ! real*8 ::dist_table(nlon, nlat, nlat)  # table for distances for points at (0, grd_lats) to points at map
  
  ! integer ::nobs # observation number 
  ! real*8  ::olon(nobs), olat(nobs) # observation locations
  ! 
  ! real*8  :: wlen  # weighting length 
  ! real*8  ::lon_period  # cycle of longitude (such as 360 etc)
  ! 
  
  !  Outputs: 
  !  ----------------------------------------------
  !  real*8:: wgt(nlon, nlat) # weigted observation impacts
  

  implicit none
  ! inputs
  ! map and distance between each point 
  
  integer, intent(in)::nlon
  integer, intent(in)::nlat
  real*8, intent(in)::grd_lon(nlon)
  real*8, intent(in)::grd_lat(nlat)
  real*8, intent(in)::dist_table(nlon, nlat, nlat) 
  
  ! observations 
  
  integer, intent(in)::nobs
  real*8, intent(in)::olon(nobs)
  real*8, intent(in)::olat(nobs)
  
  ! length

  real*8, intent(in)::wlen 
  
  ! period 

  real*8, intent(in)::lon_period 

  ! output 
  real*8, intent(out)::wgt(nlon, nlat)

  ! local 
  
  integer, allocatable:: plon(:), plat(:)
  real*8, allocatable::dlon(:)
  integer::ilon, ilat, iobs
  integer::ix 
  real*8::exp_wgt
  
  
  allocate(plon(nlon))
  allocate(dlon(nlon))
  ! plat is for olat 
  
  allocate(plat(nobs))
  ! search for closest lat grid for olat 
  
  call get_closest_point(nlat, grd_lat, nobs, olat, plat)

  wgt=0.0
  
  do iobs=1, nobs  ! loop over nobs
     
     dlon=grd_lon-olon(iobs)
     ! shift points
     
     do ix=1, nlon
        if (dlon(ix)<grd_lon(1)) then 
           ! shift by one cycle right 
           dlon(ix)=dlon(ix)+lon_period
        else if (dlon(ix)>grd_lon(nlon)) then
           ! shift by one cycle right 
           dlon(ix)=dlon(ix)-lon_period
        end if
     
     end do ! ix 
     
     call get_closest_point(nlon, grd_lon, nlon, dlon, plon)
     ! fill the map point by point 
     
     do ilon=1, nlon
        do ilat=1, nlat
           ! table index 
           ! plat(iobs)=postion in grd_lat of olat  
           ! ilat--> target postion in grd_lat
           ! plon(ilon)--> target position in grd_lon (after being shifted by olon)
           
           exp_wgt=dist_table(plon(ilon), ilat, plat(iobs))
           exp_wgt=exp_wgt/wlen
           
           if (exp_wgt<12.0) then
              exp_wgt=exp(-exp_wgt)
              wgt(ilon, ilat)=wgt(ilon, ilat)+exp_wgt
           end if
           
           
        end do ! ilat
     end do ! ilon
     
  end do ! iobs 
  
  deallocate(plon)
  deallocate(plat)
  deallocate(dlon)
  
  
end subroutine get_obs_wgt_map


