subroutine get_circle_distance(lons, lats, gc_dist, n)
! this is a routine to calculate the distance of a set of points 
! 
! Inputs:
! ------------------------------------------------------
! integer::n # Set size 
! 
! real*8::lons(n), lats(n) # longitude and latitudes of  the set  
! 
! Outputs
!---------------------------------------------------------
! real*8 ::gc_dist(n,n)  # distance bewteen points in km


implicit none
real*8, intent(in)::lons(n), lats(n)
real*8,  intent(out)::gc_dist(n,n)
integer, intent(in)::n

real*8, parameter::   g0 = 9.80665,&
     pi                       =  3.141592653589793115997963,&
     earth_major_axis         = 6378137.0, & 
     earth_minor_axis         = 6356752.3141, &
     earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
     (earth_minor_axis*earth_minor_axis)), &
     crad=pi/180.0
! local variables
real*8, allocatable::rlon(:), rlat(:)
real*8::ran, rearth, dlon
integer:: ix, iy 

allocate(rlon(n), rlat(n))
! convert to radius

rlon(:)=crad*lons(:)
rlat(:)=crad*lats(:)
ran=0.0
rearth=0.5*(earth_major_axis+earth_minor_axis)
gc_dist(:,:)=0.0
do ix=1, n
   do iy=ix+1, n
      dlon=rlon(ix)-rlon(iy)
      
      ran=sin(rlat(ix))*sin(rlat(iy))+cos(rlat(ix))*cos(rlat(iy))*cos(dlon)
      ran=acos(ran)
      ! print *, ran/crad
      ! to km 
      gc_dist(ix, iy)=rearth*ran/1000.0
      gc_dist(iy, ix)=gc_dist(ix, iy)
   end do
end do
! remove arrays
deallocate(rlon, rlat)
end subroutine get_circle_distance


subroutine set_circle_distance_table(lons, lats, gc_dist, nx, ny)
! this is a routine to calculate the distance of a set of points 
! 
! Inputs:
! ------------------------------------------------------
! integer::n # Set size 
! 
! real*8::lons(n), lats(n) # longitude and latitudes of  the set  
! 
! Outputs
!---------------------------------------------------------
! real*8 ::gc_dist(n,n)  # distance bewteen points

implicit none
integer, intent(in)::nx, ny

real*8, intent(in)::lons(nx), lats(ny)
real*8,  intent(out)::gc_dist(nx, ny, ny)


real*8, parameter::   g0 = 9.80665,&
     pi                       =  3.141592653589793115997963,&
     earth_major_axis         = 6378137.0, & 
     earth_minor_axis         = 6356752.3141, &
     earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
     (earth_minor_axis*earth_minor_axis)), &
     crad=pi/180.0
! local variables

real*8, allocatable::rlon(:), rlat(:)
real*8::ran, rearth, lat_py
integer:: ix, iy, ipy 

allocate(rlon(nx), rlat(ny))
! convert to radius
rlon(:)=crad*lons(:)
rlat(:)=crad*lats(:)
ran=0.0
rearth=0.5*(earth_major_axis+earth_minor_axis)
gc_dist(:,:,:)=0.0

do ipy=1, ny ! 
   ! start point 
   lat_py=rlat(ipy) 
   do ix=1, nx
      do iy=1, ny

         ran=sin(rlat(iy))*sin(lat_py)+cos(rlat(iy))*cos(lat_py)*cos(rlon(ix))
         ran=acos(ran)
         ! print *, ran/crad
         gc_dist(ix, iy, ipy)=rearth*ran/1000.0
      end do ! ix
   end do ! iy
end do ! ipy

! remove arrays
deallocate(rlon, rlat)
end subroutine set_circle_distance_table


subroutine get_circle_distance_xy(lonx, latx, lony, laty, gc_dist, nx, ny)
  implicit none
  ! get circle distance for two set of points
  !
  ! Inputs:
  ! ----------------------------
  ! integer::nx, ny # Set X size and Set Y size
  ! 
  ! real*8::lonx(nx), latx(nx) # longitude and latitudes of  Set X 
  ! real*8, intent(in)::lony(ny), laty(ny) #  longitude and latitudes of  Set X
  ! 
  ! Outputs
  !-----------------------------------
  ! real*8 ::gc_dist(nx,ny)  # distance bewteen points in Set X and Set Y  in km



  integer, intent(in)::nx, ny
  real*8, intent(in)::lonx(nx), latx(nx)
  real*8, intent(in)::lony(ny), laty(ny)
  
  real*8,  intent(out)::gc_dist(nx,ny)
  
  real*8, parameter::   g0 = 9.80665,&
       pi                       =  3.141592653589793115997963,&
       earth_major_axis         = 6378137.0, & 
       earth_minor_axis         = 6356752.3141, &
       earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
       (earth_minor_axis*earth_minor_axis)), &
       crad=pi/180.0
  ! local variables
  real*8, allocatable::rlonx(:), rlatx(:),rlony(:), rlaty(:)
  real*8::ran, rearth
  integer:: ix, iy 
  allocate(rlonx(nx), rlatx(nx), rlony(ny), rlaty(ny))
  ! convert to radius
  rlonx(:)=crad*lonx(:)
  rlatx(:)=crad*latx(:)
  
  rlony(:)=crad*lony(:)
  rlaty(:)=crad*laty(:)
  
  ran=0.0
  rearth=0.5*(earth_major_axis+earth_minor_axis)
  gc_dist(:,:)=0.0
  do ix=1, nx
     do iy=1, ny
        ran=sin(rlatx(ix))*sin(rlaty(iy))+cos(rlatx(ix))*cos(rlaty(iy))*cos(rlonx(ix)-rlony(iy))
        ran=acos(ran)
        ! print *, ran/crad
        gc_dist(ix, iy)=rearth*ran/1000.0
     end do
  end do
  ! remove arrays
  deallocate(rlonx, rlatx, rlony, rlaty)
end subroutine get_circle_distance_xy



subroutine gen_conv(cor_len, lwi, sf_err, dist, err_conv,n)
  ! calculate error correlation for one set of  points 
  ! the correlation between point A and B are calculated Err_A*Err_B*exp(-dist_AB/cor_len)
  !
  ! Inputs: 
  !------------------------------------------------------
  !  integer::n # size of set X and Y 
  ! 
  !  real*8::lwi  # correlation length 
  !  real*8::sf_err(n) # error   
  !  real*8::dist(n,n) # distance between Set X and Y 
  !  
  !  Outs
  !-----------------------------------------------------
  ! real*8::err_conv(n,n)  # error covariance 
  
  implicit none
  real*8, intent(in)::cor_len(3)
  integer*4, intent(in)::lwi(n)
  real*8,    intent(in)::sf_err(n)
  real*8,    intent(in)::dist(n,n)
  real*8,    intent(out)::err_conv(n,n)
  integer*4, intent(in)::n
  ! local variables
  integer::ix, iy, sfx,sfy
  real*8::errx, erry, corl 
  err_conv=0.0
  do ix=1, n
     sfx=lwi(ix)
     errx=sf_err(ix)
     do iy=ix,n 
        sfy=lwi(iy)
        erry=sf_err(iy)
        if (sfx==sfy) then
           corl=cor_len(sfx+1)
           err_conv(ix, iy)=sqrt(errx*erry)*exp(-dist(ix,iy)/corl)
           err_conv(iy, ix)=err_conv(ix,iy)
        end if
     end do
  end do
end subroutine gen_conv



subroutine gen_conv_xy(errx, nx, erry, ny, dist, cor_len, err_conv)
  ! calculate error correlation for two sets of points X and Y
  !  the correlation between point A and B are calculated Err_A*Err_B*exp(-dist_AB/cor_len)
  !  Inputs: 
  !
  !  integer::nx, ny # size of set X and Y 
  ! 
  !  real*8::cor_len  # correlation length 
  !  real*8::errx(nx) # error for Set X  
  !  real*8::erry(ny)  # error for Set Y 
  
  !  real*8::dist(nx,ny) # distance between Set X and Y 
  !  real*8::err_conv(nx,ny)  # error covariance matrix 
  !  
  !  Outs
  !-----------------------------------------------------
  ! real*8::err_conv(nx,ny)  # error covariance 
  

  implicit none
  
  integer*4, intent(in)::nx, ny
  real*8, intent(in)::cor_len
  real*8,    intent(in)::errx(nx)
  real*8,    intent(in)::erry(ny)
  
  real*8,    intent(in)::dist(nx,ny)
  real*8,    intent(out)::err_conv(nx,ny)
  
  ! local variables
  
  integer::ix,iy
  real*8::xerr, yerr, sf 
  ! initialization
  
  err_conv=0.0
  
  do ix=1, nx
     xerr=errx(ix)
     do iy=ix,ny 
        yerr=erry(iy)
        sf=dist(ix,iy)/cor_len
        
        if (sf>10.0) then
           sf=0.0
        else
           sf=exp(-sf)
        end if
        
        err_conv(ix, iy)=sqrt(xerr*yerr)*sf
        err_conv(iy, ix)=err_conv(ix,iy)
        
        
     end do ! iy
  end do  ! ix
  
end subroutine gen_conv_xy



 
