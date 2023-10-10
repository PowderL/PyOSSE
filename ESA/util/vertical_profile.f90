! functions for interpolation 
! 1. getpos
! 2. getwgt #  get weighting factors  
! 3. getwgt_mask  # version of getwgt for masked array 
! 4. get_vertical_wgt_0d  # wrapper for getwgt_mask for single pressure profile. 
! 5. get_vertical_wgt_1d  # As get_vertical_wgt_0d but for pressure profiles along a track
! 6. get_vertical_wgt_2d  #  As get_vertical_wgt_0d but for profiles over (x,y) region
! 7. prof_vertical_intpl_0d # interpolate single tracer profile using position and weight arrays from get_vertical_wgt.
! 8. prof_vertical_intpl_1d # interpolate tracer profiles along a track using position and weight arrays from get_vertical_wgt.

! 9. prof_vertical_intpl_2d # 2D version of prof_vertical_intpl_1d
! 10.prof_vertical_intpl_0d_em: # ensemble version of prof_vertical_intpl_0d

! 11 prof_vertical_intpl_1d_em :  #  ensemble version of prof_vertical_intpl_1d
! grid_ob

subroutine getpos(x0, x, nx0, nx, pl, pr)
! find the lower (pl) and upper (pr) boundary of value x in axis x0
! x0 ------in--------- array of size nx0 in ascending order. 
! x  ------in--------- array of size nx
! pl ------out-------- array of size nx. Location of lower boundary
! pr ------out-------- array of size nx. Location of upper boundary


implicit none 
integer, intent(in)::nx, nx0  
real*8, intent(in)::x0(nx0), x(nx)
integer*4, intent(out)::pl(nx), pr(nx)   
integer::i, j, l

j=1
do i=1, nx
   do l=j, nx0
      if (x0(l)>=x(i)) exit
   end do
   
   if (l==1) then 
      pl(i)=0
      pr(i)=0
   else if (l>nx0) then  
      pl(i)=nx0-1
      pr(i)=nx0-1
   else
      pl(i)=pr(i)-1
      pr(i)=l-1
   end if
   j=l
end do

end subroutine getpos

subroutine getwgt(nx0,x0, nx, x,  pl, pr,wgt, xmask)

! Find the lower (pl) and upper (pr) boundary of value x at coordinate x0, 
! and the weighting factor for interpolate. 
! Inputs:
! 1. Intetger::nx0 # size of x 
! 2.     real*8:: x0(nx0)   # original coordiantes
! 3.     Integer::nx   # size of the new  coordinate.
! 4.     real*8 ::x(nx) # new coordinate
! 5.     real*8 :: mask_val # missing and bad data

! Outputs:
! 1. integer*4::pl(nx)  # lower boundary
! 2. integer*4::pr(nx)  # upper boundary
! 3. real*8   ::wgt(nx) # weight for lower boundary

implicit none 
! in
integer, intent(in)::nx0
real*8, intent(in)::x0(nx0)
integer, intent(in)::nx
real*8, intent(in)::x(nx)
real*8, intent(in), optional::xmask

! out
integer*4, intent(out)::pl(nx), pr(nx)
real*8, intent(out)::wgt(nx)

! loca; variables

integer::i, l
real*8::dx
real*8::fmask
fmask=-999.0
if (present(xmask)) fmask=xmask

pl(:)=fmask
pr(:)=fmask

wgt(:)=0.0


! loop over the nx 


do i=1, nx
   if (x(i)==fmask) then
      ! bad value
      wgt(i)=0.0
      pl(i)=fmask
      pr(i)=fmask
      
   else
      ! seek boundary 
      do l=1, nx0
         if (x0(l)>=x(i)) exit
      end do
   
      ! if l is the  head or  bottom point 
      
      if (l==1) then ! x(i)<x0(1) 
         pl(i)=0
         pr(i)=0
         wgt(i)=1.0
      
      else if (l>nx0) then  ! x(i)>x0(nx0)  
         pl(i)=nx0-1
         pr(i)=nx0-1
         wgt(i)=1.0
      
      else   
         !  if  x0(1) < x(i) < x0(nx)
         
         pr(i)=l-1  ! right boundary;  -1 to python index 
         
         pl(i)=pr(i)-1 ! left 
         dx=x0(pr(i)+1)-x0(pl(i)+1)
         wgt(i)=0.0
         
         if (abs(dx)>0) wgt(i)=(x0(pr(i)+1)-x(i))/dx
      
      end if
   end if
end do

end subroutine getwgt

subroutine getwgt_mask(nx0, x0, nx, x,pl, pr,wgt, mask_val)

  ! Find the lower (pl) and upper (pr) boundary for value x in the original coordinate x0, 
  ! Inputs:
  !-------------------------------------------------
  ! 1.     Intetger::nx0  # size of x 
  ! 2.     real*8:: x0(nx0)   # original coordiantes
  ! 3.     Integer::nx   # size of the new  coordinate.
  ! 4.     real*8 ::x(nx) # new coordinate
  ! 5.     real*8 ::mask_val  #  mask for bad or missing value 
  !
  ! Outputs:
  !----------------------------------------------------
  ! 1. integer*4::pl(nx)  # lower boundary
  ! 2. integer*4::pr(nx)  # upper boundary
  ! 3. real*8   ::wgt(nx) # weight for lower boundary

  implicit none 
  integer, intent(in)::nx, nx0
  real*8, intent(in)::x0(nx0), x(nx)
  real*8, intent(in), optional::mask_val
  ! out
  integer*4, intent(out)::pl(nx), pr(nx)
  real*8, intent(out)::wgt(nx)
  real*8::fmask
  
  ! locat
  integer::i, j, l

  integer::ist, iend

  
  fmask=-999.0
  if (present(mask_val)) fmask=mask_val
  

  pl(:)=fmask
  pr(:)=fmask
  wgt(:)=0.0

  
  ! decide lower and higher boundary in use
  
  
  ist=nx0
  iend=0
  
  do i=1, nx0
     if (x0(i)/=fmask) then
        if (i<ist) ist=i
        if (i>iend) iend=i
     end if
  end do
  
  ! calculate left-right  boundary and interpolation weight for each point
  ! the weighting factor is for left-boundary
  
  
  do i=1, nx ! points
     
     if (x(i)==fmask) then
        pl(i)=fmask
        pr(i)=fmask
        wgt(i)=0.0
        
     else
        ! get boundary 
        
        j=0
        do l=ist, iend ! 
           
           if (x0(l)/=fmask) then 
              ! reach right  boundary
              
              if (x0(l)>=x(i)) exit
              
              ! otherwise move left boundary
              
              j=l
           end if
        end do
        ! if l is head or bottom
        
        if (l==ist) then

           pl(i)=l-1
           pr(i)=l-1
           wgt(i)=1.0
           
        else if (l>iend) then  
            
           pl(i)=iend-1
           pr(i)=iend-1
           wgt(i)=1.0
           
           
        else
           ! l is upper boundary, j is the lower boundary 
           ! -1 to python index
           
           pr(i)=l-1 
           pl(i)=j-1

           wgt(i)=(x0(l)-x(i))/(x0(l)-x0(j))
        end if
     end if
  end do
  
end subroutine getwgt_mask


subroutine get_vertical_wgt_0d(nz0, z0, nz, z, pl, pr, wgt, mask_val)
  
  ! Position and weight for vertical interpolation of profile at  nx points
  ! Inputs: 
  ! 1.   integer*4 ::nx  #  number of points 
  ! 2    integer*4 ::nz0 #  number of the original vertical levels      
  ! 3.   real*8    ::z0(nx, nz0) # original altitudes (or log(pressure)) in ascending order
  ! 4    integer*4 ::nz #  number of the new vertical levels
  ! 5.    real*8    ::z(nx, nz) # altitudes (or log(pressure))   
  ! Outputs: 
  ! 1. integer4:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 2. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old griwd
  ! 3. real*8  :: wgt(nx, nz)  #  weighting factors for the lower boundaries.
  
  implicit none
  integer, intent(in)::nz0, nz
  real*8, intent(in)::z0(nz0)
  real*8, intent(in)::z(nz)
  integer*4, intent(out)::pl(nz)
  integer*4, intent(out)::pr(nz)
  real*8, intent(out)::wgt(nz)
  real*8, optional, intent(in)::mask_val
  
  ! local variables
  integer::i
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  real*8::fmask
  
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  if (present(mask_val)) then
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  call getwgt_mask(nz0, z0(:), nz, z(:), & 
       pl1d, pr1d,wgt1d, fmask)
  pl(:)=pl1d(:)
  pr(:)=pr1d(:)
  wgt(:)=wgt1d(:)
  
  deallocate(pl1d, pr1d, wgt1d)
  
  
end subroutine get_vertical_wgt_0d




subroutine get_vertical_wgt_1d(nx, nz0, z0, nz, z, pl, pr, wgt, mask_val)
  
  ! Position and weight for vertical interpolation of profile at  nx points
  ! Inputs: 
  ! 1.   integer*4 ::nx  #  number of points 
  ! 2    integer*4 ::nz0 #  number of the original vertical levels      
  ! 3.   real*8    ::z0(nx, nz0) # original altitudes (or log(pressure)) in ascending order
  ! 4    integer*4 ::nz #  number of the new vertical levels
  ! 5.    real*8    ::z(nx, nz) # altitudes (or log(pressure))   
  ! Outputs: 
  ! 1. integer4:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 2. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old griwd
  ! 3. real*8  :: wgt(nx, nz)  #  weighting factors for the lower boundaries.
  
  implicit none
  integer, intent(in)::nx, nz0, nz
  real*8, intent(in)::z0(nx, nz0)
  real*8, intent(in)::z(nx, nz)
  integer*4, intent(out)::pl(nx, nz)
  integer*4, intent(out)::pr(nx, nz)
  real*8, intent(out)::wgt(nx,nz)
  real*8, optional, intent(in)::mask_val
  
  ! local variables
  integer::i
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  real*8::fmask
  
  
  
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  if (present(mask_val)) then
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  ! print *,'get vertical wgt add' , nz0, nz, fmask
  
  do i=1, nx
     call getwgt_mask(nz0, z0(i,:), & 
          nz, z(i, :), pl1d, pr1d,wgt1d, fmask)
     pl(i,:)=pl1d(:)
     pr(i,:)=pr1d(:)
     wgt(i,:)=wgt1d(:)
  end do
  
  
  if (allocated(pl1d)) deallocate(pl1d)
  if (allocated(pr1d)) deallocate(pr1d)
  if (allocated(wgt1d)) deallocate(wgt1d)
  
  ! deallocate(pl1d, pr1d, wgt1d)
  
  ! print *,'get vertical wgt after remove' 
  
end subroutine get_vertical_wgt_1d


subroutine get_vertical_wgt_2d(nx, ny, nz0, z0, nz, z, pl, pr, wgt, mask_val)
  ! Position and weight for vertical interpolation of profile at  nx points
  ! Inputs: 
  ! 1,2.   integer*4 ::nx, ny  #  number of points 
  ! 3    integer*4 ::nz0 #  number of the original vertical levels      
  ! 4.   real*8    ::z0(nx, ny, nz0) # original altitudes (or log(pressure)) in ascending order
  ! 5    integer*4 ::nz #  number of the new vertical levels
  ! 6.    real*8    ::z(nx, ny, nz) # altitudes (or log(pressure))   
  ! Outputs: 
  ! 1. integer4:: pl(nx, ny, nz)  #  lower boundaries of each new points at old grid  
  ! 2. integer4:: pr(nx, ny, nz)  #  upper boundaries of each new points at old griwd
  ! 3. real*8  :: wgt(nx, ny, nz)  #  weighting factors for the lower boundaries.
  
  implicit none
  integer, intent(in)::nx, ny, nz0, nz
  real*8, intent(in)::z0(nx, ny, nz0)
  real*8, intent(in)::z(nx, ny, nz)
  real*8, intent(in), optional::mask_val

  integer*4, intent(out)::pl(nx, ny, nz)
  integer*4, intent(out)::pr(nx, ny, nz)
  real*8, intent(out)::wgt(nx,ny, nz)
  ! local variables
  integer::i, j
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  real*8::fmask
  
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  if (present(mask_val)) then
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  do i=1, nx
     do j=1, ny
        call getwgt_mask(nz0, z0(i,j, :), & 
             nz, z(i,j, :), pl1d, pr1d,wgt1d, fmask)
        pl(i,j,:)=pl1d(:)
        pr(i,j,:)=pr1d(:)
        wgt(i,j,:)=wgt1d(:)
     end do
  end do
  
  deallocate(pl1d, pr1d, wgt1d)
  
  
end subroutine get_vertical_wgt_2d

subroutine prof_vertical_intpl_0d(nz, pl, pr, wgt,& 
     nz0, prof, out_prof, mask_val)

  ! Interpolate profiles at nx points to new vertical levels. 
  
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points in  x and y directions 
  ! 2. integer*4 ::nz #  number of new levels 
  ! 3. real*8:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 4. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old grid   ! 5. real*8  :: wgt(nx,nz)  #  weighting factors for the lower boundaries.
  ! 6. integer*4 ::nz0 #  number of new levels
  ! 7. real*8:: prof(nx, nz0)  #  original profiles
  ! 8. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, nz)
  

  implicit none
  ! in
  integer, intent(in)::nz
  real*8, intent(in)::wgt(nz)
  integer*4, intent(in)::pl(nz), pr(nz)
  integer, intent(in)::nz0
  real*8, intent(in)::prof(nz0)
  real*8, intent(in), optional::mask_val
  ! out
  
  real*8, intent(out)::out_prof(nz)
  
  ! local
  integer::i, k, kl,kr 
  integer::imask
  real*8:: fmask

  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if
  
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  
  do k=1, nz
     kl=pl(k)
     kr=pr(k)
     ! print *, 'k', kl, kr
     if ((kl/=imask).and.(kr/=imask)) then 
        kl=kl+1
        kr=kr+1
        
        out_prof(k)=wgt(k)*prof(kl)+(1.0-wgt(k))*prof(kr)
        ! put the filled value if original data are invalid m
        if ((prof(kl)==fmask).or.(prof(kr)==fmask)) out_prof(k)=fmask  
     else
        out_prof(k)=fmask
     end if
     
  end do
  
end subroutine prof_vertical_intpl_0d


subroutine prof_vertical_intpl_1d(nx, nz, pl, pr, wgt,& 
     nz0, prof, out_prof, mask_val)

  ! Interpolate profiles at nx points to new vertical levels. 
  
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points in  x and y directions 
  ! 2. integer*4 ::nz #  number of new levels 
  ! 3. real*8:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 4. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old grid   ! 5. real*8  :: wgt(nx,nz)  #  weighting factors for the lower boundaries.
  ! 6. integer*4 ::nz0 #  number of new levels
  ! 7. real*8:: prof(nx, nz0)  #  original profiles
  ! 8. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, nz)
  

  implicit none
  ! in
  integer, intent(in)::nx, nz
  real*8, intent(in)::wgt(nx, nz)
  integer*4, intent(in)::pl(nx, nz), pr(nx, nz)
  integer, intent(in)::nz0
  real*8, intent(in)::prof(nx, nz0)
  real*8, intent(in), optional::mask_val
  ! out
  
  real*8, intent(out)::out_prof(nx, nz)
  
  ! local
  integer::i, k, kl,kr 
  integer::imask
  real*8:: fmask

  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if
  
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  
  do i=1, nx
     do k=1, nz
        kl=pl(i,k)
        kr=pr(i,k)
        ! print *, 'k', kl, kr
        if ((kl/=imask).and.(kr/=imask)) then 
           kl=kl+1
           kr=kr+1
           
           out_prof(i, k)=wgt(i,k)*prof(i, kl)+(1.0-wgt(i,k))*prof(i,kr)
           ! put the filled value if original data are invalid m
           if ((prof(i, kl)==fmask).or.(prof(i, kr)==fmask)) out_prof(i, k)=fmask  
        else
           out_prof(i, k)=fmask
        end if
        
     end do
  end do
  
end subroutine prof_vertical_intpl_1d

subroutine prof_vertical_intpl_2d(nx, ny, nz, pl, pr, wgt,& 
     nz0, prof, out_prof, mask_val)
  
  ! Interpolate profiles at nx X ny points to new vertical levels. 
  
  ! Inputs: 
  ! 1,2.  integer*4 ::nx, ny #  number of points in  x and y directions 
  ! 3.  integer*4 ::nz #  number of new levels 
  ! 4. real*8:: pl(nx, ny, nz)  #  lower boundaries of each new points at old grid  
  ! 5. integer4:: pr(nx, ny, nz)  #  upper boundaries of each new points at old grid   
  ! 6. real*8  :: wgt(nx,ny, nz)  #  weighting factors for the lower boundaries.
  ! 7. integer*4 ::nz0 #  number of new levels
  ! 8. real*8:: prof(nx, ny, nz0)  #  original profiles
  ! 9. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, ny, nz) # new profiles
  
     
     implicit none
     integer, intent(in)::nx, ny, nz
     real*8, intent(in)::wgt(nx, ny, nz)
     integer*4, intent(in)::pl(nx, ny, nz), pr(nx, ny, nz)
     integer, intent(in)::nz0
     real*8, intent(in)::prof(nx, ny, nz0)
     real*8, intent(in), optional::mask_val
     
     ! out 
     real*8, intent(out)::out_prof(nx, ny, nz)
     
     ! local
     
     integer::i, j, k, kl,kr 
     
     integer::imask  ! mask values 
     real*8:: fmask  ! mask values
     
     if (present(mask_val)) then
        fmask=mask_val
        imask=mask_val
     else
        fmask=-999.0
        imask=-999
     end if
     
     do i=1, nx
        do j=1, ny
           do k=1, nz
              kl=pl(i,j, k)
              kr=pr(i,j, k)
              ! print *, 'k', kl, kr
              if ((kl/=imask).and.(kr/=imask)) then
                 kl=kl+1
                 kr=kr+1
                 
                 out_prof(i, j,k)=wgt(i,j,k)*prof(i, j, kl)+(1.0-wgt(i,j,k))*prof(i,j, kr)
                 ! put the filled value if original data are invalid m
                 if ((prof(i, j, kl)==fmask).or.(prof(i, j, kr)==fmask)) out_prof(i, j,k)=fmask  
              else
                 out_prof(i, j, k)=fmask
                 
              end if
           end do
        end do
     end do
     
end subroutine prof_vertical_intpl_2d


subroutine prof_vertical_intpl_0d_em(nz, pl, pr, wgt, nz0,ne, & 
     prof, out_prof, mask_val)
  
  ! Interpolate (ne) profiles at (nx) points to new vertical levels. 
  
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points in  x and y directions 
  ! 2. integer*4 ::nz #  number of new levels 
  ! 3. real*8:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 4. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old grid   ! 5. real*8  :: wgt(nx,nz)  #  weighting factors for the lower boundaries.
  ! 6,7. integer*4 ::ne, nz0 #  number of new levels
  ! 8. real*8:: prof(nx, nz0,ne)  #  original profiles
  ! 9. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, nz,ne)
  

  
  implicit none
  
  integer*4, intent(in)::nz 
  integer*4, intent(in)::pl(nz), pr(nz)
  real*8, intent(in)::wgt(nz)
  integer*4, intent(in)::nz0
  integer*4, intent(in)::ne
  real*8, intent(in)::prof(nz0,ne)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::out_prof(nz, ne)
  ! local
  
  integer::i, k, kl,kr, ie
  integer::mask_val_int
  
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  integer::imask  ! mask values 
  real*8:: fmask  ! mask values
     
  
  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if

  
  do k=1, nz
     kl=pl(k)
     kr=pr(k)
     ! print *, 'k', kl, kr
     if ((kl/=imask).and.(kr/=imask)) then
        kl=kl+1
        kr=kr+1
        
        do ie=1, ne
           if ((prof(kl,ie)==fmask).or.(prof(kr,ie)==fmask)) then
              out_prof(k,ie)=fmask
           else
              out_prof(k,ie)=wgt(k)*prof(kl,ie)+(1.0-wgt(k))*prof(kr,ie)
           end if
        end do
     else
        out_prof(k, :)=fmask
     end if
  end do
  
end subroutine prof_vertical_intpl_0d_em


subroutine prof_vertical_intpl_1d_em(nx, nz, pl, pr, wgt, & 
     nz0, ne,  prof, out_prof, mask_val)

  ! Interpolate (ne) profiles at (nx) points to new vertical levels. 
  
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points in  x and y directions 
  ! 2. integer*4 ::nz #  number of new levels 
  ! 3. real*8:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 4. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old grid   ! 5. real*8  :: wgt(nx,nz)  #  weighting factors for the lower boundaries.
  ! 6,7. integer*4 ::ne, nz0 #  number of new levels
  ! 8. real*8:: prof(nx, nz0,ne)  #  original profiles
  ! 9. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, nz,ne)
  

  
  implicit none
  integer*4, intent(in)::nx
  integer*4, intent(in)::nz 
  integer*4, intent(in)::pl(nx, nz), pr(nx, nz)
  real*8, intent(in)::wgt(nx, nz)
  integer*4, intent(in)::nz0
  integer*4, intent(in)::ne
  real*8, intent(in)::prof(nx, nz0,ne)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::out_prof(nx, nz, ne)
  ! local
  
  integer::i, k, kl,kr, ie
  integer::mask_val_int
  
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  integer::imask  ! mask values 
  real*8:: fmask  ! mask values
     
  
  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if

  
  do i=1, nx
     do k=1, nz
        kl=pl(i, k)
        kr=pr(i,k)
        ! print *, 'k', kl, kr
        if ((kl/=imask).and.(kr/=imask)) then
           kl=kl+1
           kr=kr+1
        
           do ie=1, ne
              if ((prof(i, kl,ie)==fmask).or.(prof(i, kr,ie)==fmask)) then
                 out_prof(i, k,ie)=fmask
              else
                 out_prof(i,k,ie)=wgt(i,k)*prof(i, kl,ie)+ & 
                      (1.0-wgt(i,k))*prof(i, kr,ie)
              end if
           end do
        else
           out_prof(i,k, :)=fmask
        end if
     end do
  end do
  
end subroutine prof_vertical_intpl_1d_em



subroutine prof_vertical_intpl_2d_em(nx, ny, nz, pl, pr, wgt, nz0, ne,  & 
     prof, out_prof, mask_val)

  ! Interpolate (ne) profiles at (nx) points to new vertical levels. 
  
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points in  x and y directions 
  ! 2. integer*4 ::nz #  number of new levels 
  ! 3. real*8:: pl(nx, nz)  #  lower boundaries of each new points at old grid  
  ! 4. integer4:: pr(nx, nz)  #  upper boundaries of each new points at old grid   ! 5. real*8  :: wgt(nx,nz)  #  weighting factors for the lower boundaries.
  ! 6,7. integer*4 ::ne, nz0 #  number of new levels
  ! 8. real*8:: prof(nx, nz0,ne)  #  original profiles
  ! 9. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_prof(nx, nz,ne)
  

  
  implicit none
  integer*4, intent(in)::nx
  integer*4, intent(in)::ny
  integer*4, intent(in)::nz 
  integer*4, intent(in)::pl(nx, ny, nz), pr(nx, ny, nz)
  real*8, intent(in)::wgt(nx, ny, nz)

  integer*4, intent(in)::nz0
  integer*4, intent(in)::ne

  real*8, intent(in)::prof(nx, ny, nz0,ne)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::out_prof(nx, ny, nz, ne)
  ! local
  
  integer::i, j, k, kl,kr, ie
  integer::mask_val_int
  
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  integer::imask  ! mask values 
  real*8:: fmask  ! mask values
     
  
  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if

  
  do i=1, nx
     do j=1, ny
        do k=1, nz
           kl=pl(i,j, k)
           kr=pr(i,j,k)
           ! print *, 'k', kl, kr
           if ((kl/=imask).and.(kr/=imask)) then
               kl=kl+1
               kr=kr+1
        
               do ie=1, ne
                  if ((prof(i, j,kl,ie)==fmask).or.(prof(i, j,kr,ie)==fmask)) then
                     out_prof(i, j,k,ie)=fmask
                  else
                     out_prof(i,j,k,ie)=wgt(i,j,k)*prof(i, j,kl,ie)+(1.0-wgt(i,j,k))*prof(i, j,kr,ie)
                  end if
              end do
           else
              out_prof(i,j,k, :)=fmask
           end if
        end do
     end do
  end do
  
end subroutine prof_vertical_intpl_2d_em


subroutine grid_profile_1d(nx, x, nlvl, prof, nx0, x0, use_intpl, mask_val, & 
     counts, sum_val,  sum_square)

! Allocate profiles at poistions x to grid x0 
!
! Inputs:
!-----------------------------------------------
! 1. integer, intent(in)::nx #  number of locations
! 2. real*8,  intent(in)::x(nx) # locations
! 3. integer, intent(in)::nlvl  # number of profile levels
! 5. real*8,  intent(in):: prof(nob, nlvl) # profiles 
! 6. integer, intent(in)::nx0  # size of the grid 
! 7. real*8,  intent(in)::x0(nx0)       # grid 
! 8. integer, intent(in)::use_intpl  # if use_intpl==0, profiles will be allocated to the closest grid point 
! 9. real*8, intent(in)::mask_val  # fillings for bad or missing values

! Outputs:
!--------------------------------------------------------
! 1. real*8, intent(out)::sum_val(nx0, nlvl)   # summary of the profile at each grid boxes 
! 2. real*8, intent(out)::sum_square(nx0, nlvl) # summary of the profile^2 at each grid boxes
! 3. real*8, intent(out)::counts(nx0, nlvl) # number of profiles at each grid boxes. 


Implicit none
integer, intent(in)::nx ! number of observations
real*8, intent(in)::x(nx) ! locations 

integer, intent(in)::nlvl ! number of observations
real*8,  intent(in)::prof(nx, nlvl) ! obs profiles 

integer, intent(in)::nx0 
real*8,  intent(in)::x0(nx0)

integer, intent(in)::use_intpl
real*8, intent(in)::mask_val

! Out

real*8, intent(out)::sum_val(nx0, nlvl)
real*8, intent(out)::sum_square(nx0, nlvl)
real*8, intent(out)::counts(nx0, nlvl)

! local

integer, allocatable::pl(:), pr(:)
real*8, allocatable::wgt(:)
integer::i, l


! allocate memory
allocate(pl(nx), pr(nx), wgt(nx))

! initialize counts, sum_val and sum_square to zero 

counts=0.0
sum_val=0.0
sum_square=0.0

! get weighting of x in grid x0 

call getwgt_mask(nx0, x0, nx, x,pl, pr,wgt, mask_val)

if (use_intpl==0) then

   ! no interpolation
   
   do i=1, nx ! location 
      if ((pl(i)/=mask_val).and.(pr(i)/=mask_val).and.(wgt(i)/=mask_val)) then
         do l=1, nlvl  ! level 
            if (prof(i, l)/=mask_val) then
               if (wgt(i)>=0.5) then 
                  ! allocated to left point
                  counts(pl(i), l)=counts(pl(i), l)+1
                  sum_val(pl(i), l)=sum_val(pl(i), l)+prof(i,l)
                  sum_square(pl(i), l)=sum_square(pl(i), l)+prof(i,l)*prof(i,l)
               else
                 ! allocated to right point
                  
                  counts(pr(i), l)=counts(pr(i), l)+1
                  sum_val(pr(i), l)=sum_val(pr(i), l)+prof(i,l)
                  sum_square(pr(i), l)=sum_square(pr(i), l)+prof(i,l)*prof(i,l)
               end if
            end if
         end do ! l
      end if
   end do ! i
else
   ! intepolation (i.e, sharing between two different locations
   
   do i=1, nx ! location 
      if ((pl(i)/=mask_val).and.(pr(i)/=mask_val).and.(wgt(i)/=mask_val)) then
         do l=1, nlvl ! level 
            if (prof(i, l)/=mask_val) then
               ! contribution to left point 
               counts(pl(i), l)=counts(pl(i), l)+wgt(i)
               sum_val(pl(i), l)=sum_val(pl(i), l)+prof(i,l)*wgt(i)
               sum_square(pl(i), l)=sum_square(pl(i), l)+prof(i,l)*prof(i,l)*wgt(i)
               
               ! contribution to right point 
               
               counts(pr(i), l)=counts(pr(i), l)+(1.0-wgt(i))
               sum_val(pr(i), l)=sum_val(pr(i), l)+prof(i,l)*(1.0-wgt(i))
               sum_square(pr(i), l)=sum_square(pr(i), l)+prof(i,l)*prof(i,l)*(1.0-wgt(i))
            end if
         end do ! level 
      end if
   end do  ! location

end if

! deallocate memory 
deallocate(pl)
deallocate(pr)
deallocate(wgt)


end subroutine grid_profile_1d


subroutine grid_profile_2d(nob, x, y, nlvl, prof, nx0, x0, &
     ny0, y0, use_intpl, mask_val, & 
     counts, sum_val, sum_square)


! Allocate profiles at poistions x to grid x0 
!
! Inputs:
!-----------------------------------------------
! 1. integer, intent(in)::nob #  number of locations
! 2. real*8,  intent(in)::x(nob) # locations
! 3. integer, intent(in)::y(nob) # locations
! 4. integer, intent(in)::nlvl  # number of profile levels
! 5. real*8,  intent(in):: prof(nob, nlvl) # profiles 
! 6. integer, intent(in)::nx0  # size of x grid 
! 7. real*8,  intent(in)::x0(nx0)       # x grid 
! 8. integer, intent(in)::ny0  # size of y grid 
! 9. real*8,  intent(in)::y0(ny0)       # y grid
! 10. integer, intent(in)::use_intpl  # if use_intpl==0, profiles will be allocated to the closest grid point 
! 11. real*8, intent(in)::mask_val  # fillings for bad or missing values

! Outputs:
!--------------------------------------------------------
! 1. real*8, intent(out)::sum_val(nx0, ny0, nlvl)   # summary of the profile at each grid boxes 
! 2. real*8, intent(out)::sum_square(nx0, ny0, nlvl) # summary of the profile^2 at each grid boxes
! 3. real*8, intent(out)::counts(nx0, ny0, nlvl) # number of profiles at each grid boxes. 


Implicit none
integer, intent(in)::nob ! number of locations 
real*8, intent(in)::x(nob) ! x locations
real*8, intent(in)::y(nob) ! y locations


integer, intent(in)::nlvl ! number of observations
real*8,  intent(in)::prof(nob, nlvl) ! obs profiles 

integer, intent(in)::nx0 
real*8,  intent(in)::x0(nx0)
integer, intent(in)::ny0 
real*8,  intent(in)::y0(ny0)

integer, intent(in)::use_intpl
real*8, intent(in)::mask_val

! Out

real*8, intent(out)::sum_val(nx0, ny0, nlvl)
real*8, intent(out)::sum_square(nx0, ny0, nlvl)
real*8, intent(out)::counts(nx0, ny0, nlvl)

! local

integer, allocatable::xpl(:), xpr(:)
real*8, allocatable::xwgt(:)
integer, allocatable::ypl(:), ypr(:)
real*8, allocatable::ywgt(:)

integer::i, j, l, px, py


! allocate memory
allocate(xpl(nob), xpr(nob), xwgt(nob))
allocate(ypl(nob), ypr(nob), ywgt(nob))



! initialize counts, sum_val and sum_square to zero 

counts=0.0
sum_val=0.0
sum_square=0.0

! get weighting of x in grid x0 

call getwgt_mask(nx0, x0, nob, x,xpl, xpr,xwgt, mask_val)

! get weighting of y in grid y0 

call getwgt_mask(ny0, y0, nob, y,ypl, ypr,ywgt, mask_val)

if (use_intpl==0) then
   
   ! no interpolation
   
   do i=1, nob ! location 
      if ((xpl(i)/=mask_val).and.(xpr(i)/=mask_val).and.(xwgt(i)/=mask_val)) then
         if (xwgt(i)>=0.5) then 
            ! allocated to left

            px=xpl(i)  
         else
            ! allocated to right

            px=xpr(i)
         end if
         
         if ((ypl(i)/=mask_val).and.(ypr(i)/=mask_val).and.(ywgt(i)/=mask_val)) then
            
            if (ywgt(i)>=0.5) then
               ! allocated to left
               
               py=ypl(i)
            else
               ! allocated to right
               py=ypr(i)
            end if
            
            
            do l=1, nlvl  ! level 
               if (prof(i, l)/=mask_val) then
                  ! allocated to left point
                  counts(px,py,l)=counts(px, py, l)+1
                  sum_val(px,py, l)=sum_val(px, py, l)+prof(i,l)
                  sum_square(px, py, l)=sum_square(px, py, l)+prof(i,l)*prof(i,l)
               end if
            end do  ! l 
            
         end if
      end if
      
   end do ! i

else
   ! intepolation (i.e, sharing between two different locations
 
   do i=1, nob ! location 
      if ((xpl(i)/=mask_val).and.(xpr(i)/=mask_val).and.(xwgt(i)/=mask_val)) then
         if ((ypl(i)/=mask_val).and.(ypr(i)/=mask_val).and.(ywgt(i)/=mask_val)) then
            
            do l=1, nlvl  ! level 
               if (prof(i, l)/=mask_val) then

                  ! 0,0

                  px=xpl(i)  
                  py=ypl(i)
                  counts(px,py,l)=counts(px, py, l)+xwgt(i)*ywgt(i)
                  sum_val(px,py, l)=sum_val(px, py, l)+prof(i,l)*xwgt(i)*ywgt(i)
                  sum_square(px, py, l)=sum_square(px, py, l)+ & 
                       prof(i,l)*prof(i,l)*xwgt(i)*ywgt(i)        

                  ! 0,1
                  
                  px=xpl(i)  
                  py=ypr(i)
                  counts(px,py,l)=counts(px, py, l)+xwgt(i)*(1.0-ywgt(i))
                  sum_val(px,py, l)=sum_val(px, py, l)+prof(i,l)*xwgt(i)*(1.0-ywgt(i))
                  sum_square(px, py, l)=sum_square(px, py, l)+ & 
                       prof(i,l)*prof(i,l)*xwgt(i)*(1.0-ywgt(i))        
                  

                  ! 1,0
                  
                  px=xpr(i)  
                  py=ypl(i)
                  counts(px,py,l)=counts(px, py, l)+(1.0-xwgt(i))*ywgt(i)
                  sum_val(px,py, l)=sum_val(px, py, l)+prof(i,l)*(1.0-xwgt(i))*ywgt(i)
                  sum_square(px, py, l)=sum_square(px, py, l)+ & 
                       prof(i,l)*prof(i,l)*(1.0-xwgt(i))*(1.0-ywgt(i))

                  
                  

                  ! 1, 1
                  
                  px=xpr(i)  
                  py=ypr(i)
                  counts(px,py,l)=counts(px, py, l)+(1.0-xwgt(i))*(1.0-ywgt(i))
                  sum_val(px,py, l)=sum_val(px, py, l)+prof(i,l)*(1.0-xwgt(i))*(1.0-ywgt(i))
                  sum_square(px, py, l)=sum_square(px, py, l)+ & 
                       prof(i,l)*prof(i,l)*(1.0-xwgt(i))*(1.0-ywgt(i))
                  
                  
               end if
            end do ! l
         end if
      end if
   end do ! i
end if


! deallocate memory 
deallocate(xpl)
deallocate(xpr)
deallocate(xwgt)

deallocate(ypl)
deallocate(ypr)
deallocate(ywgt)


end subroutine grid_profile_2d


subroutine slice_profile(idx, idy, nob, nx, ny, nlvl, field, prof)
! get profiles along along a track
 
!
! Inputs:
!-----------------------------------------------
! 1. integer, intent(in)::nob #  number of locations
! 2. integer, intent(in)::idx(nob) # x locations
! 3. integer, intent(in)::idy(nob) # y locations
! 4. integer, intent(in)::nlvl  # number of profile levels
! 5. real*8,  intent(in):: field(nx, ny, nlvl) # 3D field 

! Outputs:
!--------------------------------------------------------
! 1. real*8, intent(out)::profile(nobs,nlvl)   # profile at each location


! inputs

integer, intent(in)::nob !  number of locations
integer, intent(in)::idx(nob) ! x locations
integer, intent(in)::idy(nob) ! y locations
integer, intent(in)::nlvl  ! number of profile levels
integer, intent(in)::nx, ny ! x and y size of field grid 
real*8,  intent(in):: field(nx, ny, nlvl) ! 3D field 

! Outputs:
!--------------------------------------------------------
real*8, intent(out)::prof(nob,nlvl)   ! profile at each location


! local variables

integer::ii
integer::ll

do ii=1, nob
   do ll=1, nlvl
      prof(ii, ll)=field(idx(ii), idy(ii), ll)
   end do
end do

end subroutine slice_profile

  




subroutine prof_replace_1d(pres, prof, sp, nx, nz, & 
     out_pres,out_prof,  mask_val)
  
  ! replace the value below surface pressure with mask vale
  ! 
  ! Inputs: 
  ! 1. real*8:: press(nx, nz)  #  pressure in ascending order
  ! 2. real*8:: prof(nx, nz)  #  profile   
  ! 3. real*8:: sp(nx)  # surface pressure 
  ! 4. integer ::nx #  number of points in  x and y directions 
  ! 5. integer ::nz #  number of new levels 
  ! 6. real*8:: mask_val #  mask for bad or missing values
  
  ! outputs:
  ! 1. real*8::out_pres(nx, nz) # pressure 
  ! 2. real*8::out_prof(nx, nz) # profile 
  
  

  implicit none
  ! in
  real*8, intent(in)::pres(nx, nz)
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::sp(nx)
  integer, intent(in)::nx, nz
  real*8, intent(in), optional::mask_val
  
  ! out
  real*8, intent(out)::out_prof(nx, nz)
  real*8, intent(out)::out_pres(nx, nz)
  
  ! out
  
  
  ! local
  integer::i, k, kl,kr 
  integer::imask
  real*8:: fmask
  real*8::dp, dval

  ! set filling value 
  
  if (present(mask_val)) then
     fmask=mask_val
     imask=mask_val
  else
     fmask=-999.0
     imask=-999
  end if
  
  out_pres=fmask
  out_prof=fmask
  
  ! check each pressure profile 
  
  do i=1, nx
     ! find the locations
     do k=1, nz
        kl=k
        kr=k+1
        if (pres(i, k)>fmask.and.(prof(i, k)>fmask)) then
           out_pres(i,k)=pres(i, k)
           out_prof(i,k)=prof(i,k)
           
           
           if (pres(i, k)>=sp(i)) then
              kl=k
              kr=k+1
              out_pres(i, k)=sp(i)
              if (kr<=nz) then
                 if ((pres(i,kr)>fmask).and.(prof(i, kr)>fmask)) then
                    ! a simple interpolation 
                    dp=pres(i,kr)-pres(i, k)
                    dval=prof(i,kr)-prof(i, k)
                    if (dp>0) out_prof(i, k)=out_prof(i, k)+(sp(i)-pres(i,k))*dval/dp
                 end if
              end if
              exit
              
           else
              if (k==nz) out_pres(i, k)=sp(i)
              
           end if ! pres(i,k)>sp
           
        end if ! pressure and profile data are valid 
        
     end do ! nz 
     
  end do  ! nx

  
  
end subroutine prof_replace_1d

