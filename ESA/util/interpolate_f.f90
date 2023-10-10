! functions for allocating masked array  in one masked axis 
! 1. search_upper_boundary:  right edge of x in a axis
! 2. getpos: get left-right boundaries of array x in a axis
! 3. getwgt:get left-right boundaries, and weighting factor
! 4. getwgt_nask: wrapper for getwgt

subroutine search_upper_boundary(nx0, x0, x, pr)
! find the right (upper) boundary for x within an array x0
!   
! Inputs:
!--------------------------------------------------------
! 1. integer::nx0: size of axis
! 2. real*8:: x0(nx0): axis array in ascending order. 
! 3. real::  x: value to be allocated 
! Outputs:
! 1. integer*4:: pr : left boundary -1 
! 
! Notes
! ------------------------------
! if x > x0(nx0),  pr=nx0
! if x<x0(1),   pr=-1
! otherwise  x0(pr)<x<x0(pr+1), pr is in python index

implicit none

integer, intent(in)::nx0 !  size of axis
real*8, intent(in):: x0(nx0) !  axis array in ascending order. 
real*8, intent(in) ::  x ! value to be allocated 
! Outputs:

integer*4, intent(out):: pr ! left boundary -1 

! 
! locals 
integer:: ll, lmid, lr  ! left, right ,mid 

if (x<x0(1)) then
   pr=0
else if (x>x0(nx0)) then
   pr=nx0
   
else
   ll=1
   lr=nx0
   pr=lr
   
   do 
!      print *, 'll,lr', ll, lr
      
      if (lr<=ll+1) then
         pr=lr-1
         exit
      end if
      
      lmid=0.5*(ll+lr)
      
      if (x0(lmid)>x) then 
         lr=lmid  ! right-boundary moved to middle
      else if (x0(lmid)==x) then
         pr=lmid-1
         exit
      else
         ll=lmid  ! left-boundary moved to middle 
         
      end if
   end do
   

end if

end subroutine search_upper_boundary

subroutine getpos(x0, x, nx0, nx, pl, pr, mask_val)
! find the lower (pl) and upper (pr) boundary of value x in axis x0
! Inputs:
! ------------------------------------------------------------
! 
! x0 ------in--------- array of size nx0 in ascending order. 
! x  ------in--------- array of size nx
! pl ------out-------- array of size nx. Location of lower boundary
! pr ------out-------- array of size nx. Location of upper boundary


implicit none 
integer, intent(in)::nx, nx0  
real*8, intent(in)::x0(nx0), x(nx)
real*8, intent(in), optional::mask_val
integer*4, intent(out)::pl(nx), pr(nx)   
integer::i, j
! local variables
real*8::fmask
real*8, allocatable::x0use(:)
integer, allocatable::idx(:)
integer::ist, iend, nuse
integer::cpr

fmask=-999.0

IF (present(mask_val)) then
   fmask=mask_val
!   PRINT *, 'MASK_VAL', MASK_VAL
END IF

! allocate memory 
allocate(x0use(nx0), idx(nx0))

! remove filling or bad values


ist=nx0
iend=1
nuse=0

do i=1, nx0
   if (x0(i)/=fmask) then
      if (i<ist) ist=i
      if (i>iend) iend=i
      nuse=nuse+1
      x0use(nuse)=x0(i)
      idx(nuse)=i
!      print *, i, ist, iend, fmask

   end if

end do

! print *, 'ist', 'iend', 'nuse', ist, iend, nuse

j=1
do i=1, nx
   if (x(i)==fmask) then
      pl(i)=fmask
      pr(i)=fmask
   else
      call search_upper_boundary(nuse, x0use(1:nuse), x(i), cpr)
      if (cpr==0) then  ! head
         pl(i)=ist    
         pr(i)=ist
      else if (cpr==nuse) then ! end 
         pl(i)=iend
         pr(i)=iend
      else
         pr(i)=idx(cpr+1)   ! changed to FORTRAN index
         pl(i)=idx(cpr)
      end if
      
      pl(i)=pl(i)-1      ! move to python idex
      pr(i)=pr(i)-1
   end if
end do
! release memory 

if (allocated(x0use)) deallocate(x0use)
if (allocated(idx)) deallocate(idx)

end subroutine getpos

subroutine getpos_mask_outsider(x0, x, nx0, nx, pl, pr, mask_val)
! find the lower (pl) and upper (pr) boundary of value x in axis x0
! Inputs:
! ------------------------------------------------------------
! 
! x0 ------in--------- array of size nx0 in ascending order. 
! x  ------in--------- array of size nx
! pl ------out-------- array of size nx. Location of lower boundary
! pr ------out-------- array of size nx. Location of upper boundary


implicit none 
integer, intent(in)::nx, nx0  
real*8, intent(in)::x0(nx0), x(nx)
real*8, intent(in), optional::mask_val
integer*4, intent(out)::pl(nx), pr(nx)   
integer::i, j
! local variables
real*8::fmask
real*8, allocatable::x0use(:)
integer, allocatable::idx(:)
integer::ist, iend, nuse
integer::cpr

fmask=-999.0

IF (present(mask_val)) then
   fmask=mask_val
!   PRINT *, 'MASK_VAL', MASK_VAL
END IF

! allocate memory 
allocate(x0use(nx0), idx(nx0))

! remove filling or bad values


ist=nx0
iend=1
nuse=0

do i=1, nx0
   if (x0(i)/=fmask) then
      if (i<ist) ist=i
      if (i>iend) iend=i
      nuse=nuse+1
      x0use(nuse)=x0(i)
      idx(nuse)=i
!      print *, i, ist, iend, fmask

   end if

end do

! print *, 'ist', 'iend', 'nuse', ist, iend, nuse

j=1
do i=1, nx
   if (x(i)==fmask) then
      pl(i)=fmask
      pr(i)=fmask
   else
      call search_upper_boundary(nuse, x0use(1:nuse), x(i), cpr)
      if (cpr==0) then  ! head
         pl(i)=ist    
         pr(i)=ist
         if (x(i)<x0use(0)) then
            pl(i)=fmask    
            pr(i)=fmask
         else
             pl(i)=pl(i)-1      ! move to python idex
             pr(i)=pr(i)-1
         end if
      else if (cpr==nuse) then ! end 
         pl(i)=iend
         pr(i)=iend
         if (x(i)>x0use(nuse)) then
            pl(i)=fmask    
            pr(i)=fmask
         else
             pl(i)=pl(i)-1      ! move to python idex
             pr(i)=pr(i)-1
         end if
      else
         pr(i)=idx(cpr+1)   ! changed to FORTRAN index
         pl(i)=idx(cpr)
         pl(i)=pl(i)-1      ! move to python idex
         pr(i)=pr(i)-1
      end if
      
     
   end if
end do
! release memory 

if (allocated(x0use)) deallocate(x0use)
if (allocated(idx)) deallocate(idx)

end subroutine getpos_mask_outsider

subroutine getwgt(nx0,x0, nx, x,  pl, pr,wgt, mask_val)

! Find the lower (pl) and upper (pr) boundary of value x at coordinate x0, 
! and the weighting factor for interpolate. 
! Inputs:
!------------------------------------------------
! 1. Intetger::nx0 # size of x0 
! 2.     real*8:: x0(nx0)   # original coordiantes
! 3.     Integer::nx   # size of the new  coordinate.
! 4.     real*8 ::x(nx) # new coordinate
! 5.     real*8 :: mask_val # missing and bad data

! Outputs:
!------------------------------------------------------
! 1. integer*4::pl(nx)  # lower boundary
! 2. integer*4::pr(nx)  # upper boundary
! 3. real*8   ::wgt(nx) # weight for lower boundary

implicit none 
! in
integer, intent(in)::nx0
real*8, intent(in)::x0(nx0)
integer, intent(in)::nx
real*8, intent(in)::x(nx)
real*8, intent(in), optional::mask_val

! out
integer*4, intent(out)::pl(nx), pr(nx)
real*8, intent(out)::wgt(nx)

! loca; variables

! local variables
real*8::fmask
real*8, allocatable::x0use(:)
integer, allocatable::idx(:)
integer::i,j,ist, iend, nuse
integer::cpr
real*8::dx, dx0


fmask=-999.0

if (present(mask_val)) fmask=mask_val
! allocate memory 
allocate(x0use(nx0), idx(nx0))

! remove filling or bad values


ist=nx0
iend=1
nuse=0

do i=1, nx0
   if (x0(i)/=fmask) then
      if (i<ist) ist=i
      if (i>iend) iend=i
      nuse=nuse+1
      x0use(nuse)=x0(i)
      idx(nuse)=i


   end if

end do


j=1
do i=1, nx
   if (x(i)==fmask) then
      pl(i)=fmask
      pr(i)=fmask
      wgt(i)=0.0
      
   else
      
      
      call search_upper_boundary(nuse, x0use(1:nuse), x(i), cpr)
      if (cpr==0) then  ! head
         pl(i)=ist    
         pr(i)=ist
         wgt(i)=1.0
      else if (cpr==nuse) then ! end 
         pl(i)=iend
         pr(i)=iend
         wgt(i)=0.0
      else
         pr(i)=idx(cpr+1)   ! changed to FORTRAN index
         pl(i)=idx(cpr)
         dx0=x0use(cpr+1)-x0use(cpr)
         dx=x0use(cpr+1)-x(i)
         wgt(i)=0.0
         if (dx0>0) wgt(i)=dx/dx0
         
         
      end if
      
      pl(i)=pl(i)-1      ! move to python idex
      pr(i)=pr(i)-1
   end if
      
end do
! release memory 

if (allocated(x0use)) deallocate(x0use)
if (allocated(idx)) deallocate(idx)

end subroutine getwgt


subroutine getwgt_mask_outsider(nx0,x0, nx, x,  pl, pr,wgt, mask_val)

! Find the lower (pl) and upper (pr) boundary of value x at coordinate x0, 
! and the weighting factor for interpolate. The outsiders will be masked
! Inputs:
!------------------------------------------------
! 1. Intetger::nx0 # size of x0 
! 2.     real*8:: x0(nx0)   # original coordiantes
! 3.     Integer::nx   # size of the new  coordinate.
! 4.     real*8 ::x(nx) # new coordinate
! 5.     real*8 :: mask_val # missing and bad data

! Outputs:
!------------------------------------------------------
! 1. integer*4::pl(nx)  # lower boundary
! 2. integer*4::pr(nx)  # upper boundary
! 3. real*8   ::wgt(nx) # weight for lower boundary

implicit none 
! in
integer, intent(in)::nx0
real*8, intent(in)::x0(nx0)
integer, intent(in)::nx
real*8, intent(in)::x(nx)
real*8, intent(in), optional::mask_val

! out
integer*4, intent(out)::pl(nx), pr(nx)
real*8, intent(out)::wgt(nx)

! loca; variables

! local variables
real*8::fmask
real*8, allocatable::x0use(:)
integer, allocatable::idx(:)
integer::i,j,ist, iend, nuse
integer::cpr
real*8::dx, dx0


fmask=-999.0

if (present(mask_val)) fmask=mask_val
! allocate memory 
allocate(x0use(nx0), idx(nx0))

! remove filling or bad values


ist=nx0
iend=1
nuse=0

do i=1, nx0
   if (x0(i)/=fmask) then
      if (i<ist) ist=i
      if (i>iend) iend=i
      nuse=nuse+1
      x0use(nuse)=x0(i)
      idx(nuse)=i


   end if

end do


j=1
do i=1, nx
   if (x(i)==fmask) then
      pl(i)=fmask
      pr(i)=fmask
      wgt(i)=0.0
      
   else
      call search_upper_boundary(nuse, x0use(1:nuse), x(i), cpr)
      if (cpr==0) then  ! head
         pl(i)=ist    
         pr(i)=ist
         wgt(i)=1.0
         ! outsider
         if (x(i)<x0use(0)) then
            pl(i)=fmask    
            pr(i)=fmask
            wgt(i)=0.0
         else
             pl(i)=pl(i)-1      ! move to python idex
             pr(i)=pr(i)-1
         end if
         
      else if (cpr==nuse) then ! end 
         pl(i)=iend
         pr(i)=iend
         wgt(i)=0.0
         if (x(i)>x0use(nuse)) then ! outsider
            pl(i)=fmask    
            pr(i)=fmask
         else
             pl(i)=pl(i)-1      ! move to python idex
             pr(i)=pr(i)-1
         end if

      else
         pr(i)=idx(cpr+1)   ! changed to FORTRAN index
         pl(i)=idx(cpr)
         dx0=x0use(cpr+1)-x0use(cpr)
         dx=x0use(cpr+1)-x(i)
         wgt(i)=0.0
         if (dx0>0) wgt(i)=dx/dx0
         pl(i)=pl(i)-1      ! move to python idex
         pr(i)=pr(i)-1
         
      end if
      
  
   end if
      
end do
! release memory 

if (allocated(x0use)) deallocate(x0use)
if (allocated(idx)) deallocate(idx)

end subroutine getwgt_mask_outsider




subroutine getwgt_mask(nx0, x0, nx, x,pl, pr,wgt, mask_val)

  ! wrapper of getwet
  ! Inputs:
  !-------------------------------------------------
  ! 1.     Intetger::nx0  # size of x0 
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
  call getwgt(nx0,x0, nx, x,  pl, pr,wgt, fmask)
  
end subroutine getwgt_mask


subroutine get_interpl(nx0,f0, nx, pl, pr,wgt, f, mask_val)
  
! interpolate f0 to locations defined by pl, pr and wgt 
! Inputs:
!------------------------------------------------
! 1. Intetger::nx0 # size of f0 
! 2.     real*8:: f0(nx0)   # original values
! 3.     Integer::nx   # size of the new  coordinate.
! 4. integer*4::pl(nx)  # lower boundary
! 5. integer*4::pr(nx)  # upper boundary
! 6. real*8   ::wgt(nx) # weight for lower boundary
! 7.     real*8 :: mask_val # missing and bad data

! Outputs:
! 1.     real*8 ::f(nx) # f0 at new coordinate
  implicit none

  integer, intent(in)::nx, nx0
  real*8, intent(in)::f0(nx0)
  integer*4, intent(in)::pl(nx), pr(nx)
  real*8, intent(out)::f(nx)
  real*8, intent(in), optional::mask_val
  
  ! out

  real*8, intent(in)::wgt(nx)
  real*8::fmask
  integer::i, ll, pp
  fmask=-999.0
  if (present(mask_val)) fmask=mask_val
  do i=1, nx
     ll=pl(i)
     pp=pr(i)
     if ((ll==fmask).or.(pp==fmask)) then
        f(i)=fmask
     else
        ll=ll+1
        pp=pp+1
        
        if (f0(ll)==fmask) then
           f(i)=fmask
        else if (f0(pp)==fmask) then
           f(i)=fmask
        else if (wgt(i)==fmask) then
           f(i)=fmask
        else
!           print *, ll, pp, wgt(i)
!           print *, f0(ll), f0(pp)
           
           f(i)=wgt(i)*f0(ll)+(1.0-wgt(i))*f0(pp)
        end if
       
     end if ! ll, pp 
  end do ! nx 
end subroutine get_interpl



subroutine get_interpl_1d(nx0,nz,f0, nx, pl, pr,wgt, f, mask_val)

  
! interpolate f0(nx0, nz) to locations defined by pl, pr and wgt 
! Inputs:
!------------------------------------------------
! 1. Intetger::nx0, nz # size of f0 
! 2.     real*8:: f0(nx0, nz)   # original values
! 3.     Integer::nx   # size of the new  coordinate.
! 4. integer*4::pl(nx)  # lower boundary
! 5. integer*4::pr(nx)  # upper boundary
! 6. real*8   ::wgt(nx) # weight for lower boundary
! 7.     real*8 :: mask_val # missing and bad data

! Outputs:
! 1.     real*8 ::f(nx) # f0 at new coordinate
  implicit none

  integer, intent(in)::nx, nz, nx0
  real*8, intent(in)::f0(nx0,nz)
  integer*4, intent(in)::pl(nx), pr(nx)
  real*8, intent(out)::f(nx, nz)
  real*8, intent(in), optional::mask_val
  
  ! out

  real*8, intent(in)::wgt(nx)
  real*8::fmask
  integer::i, j, ll, pp
  fmask=-999.0
  if (present(mask_val)) fmask=mask_val
  do i=1, nx
     ll=pl(i)
     pp=pr(i)
     if ((ll==fmask).or.(pp==fmask)) then
        f(i,:)=fmask
     else
        ll=ll+1
        pp=pp+1
        do j=1, nz
           
           if (f0(ll, j)==fmask) then
              f(i,j)=fmask
           else if (f0(pp,j)==fmask) then
              f(i,j)=fmask
           else if (wgt(i)==fmask) then
              f(i,j)=fmask
           else
              f(i,j)=wgt(i)*f0(ll,j)+(1.0-wgt(i))*f0(pp,j)
           end if
        end do ! j
        
     end if ! ll, pp 
  end do ! nx 
end subroutine get_interpl_1d

