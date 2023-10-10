! function for sample model grid 
! 
! 
subroutine sample_field(nlat, nlon, nl, gp, & 
     nxp, lonp1, lonp2, latp1, latp2, & 
     w1, w2, w3, w4, & 
     gp_prof, mask_val)

  ! Get model profiles for multiple tracers  along the track 
  !
  ! Inputs
  !-------------------------------------------
  ! 1. integer*4::  nlat: size of latitude 
  ! 2. integer*4::  nlon: size of longitude 
  ! 3. integer*4::  nl:  vertical levels
  ! 4. real*8   :: gp(nlon, nlat, nl): gridded model data
  ! 5 .integer*4::  nxp: Point number along the track
  ! 6,7. integer*4::  lonp1(nxp), lonp2(nxp):  left-right longitude indexes
  ! 8, 9, integer*4:: latp1(nxp), lonp2(nxp):  left-right latitude indexes
  ! 9, 10, 11, 12.  real*8 :: w1(nxp), w2(nxp), w3(nxp), w4(nxp): 
  ! ---weighting factors for four corners around each points. 
  ! 13. real*8:: mask_val: filling value for missing or bad data 

  ! Outputs
  !--------------------------------------
  ! real*8, intent(out)    :: gp_prof(nxp, nl) : profiles along the track  !
  !
  ! Notes:
  !------------------------------------------------------
  ! 1. if the impact of bad or missing model values larger than 50% to a track point, 
  ! ---model value at this point will be set to mask value 
  
  
  
  
  implicit none
  
  integer*4, intent(in)  :: nlat, nlon, nl
  real*8, intent(in)     :: gp(nlon, nlat, nl)
  integer*4, intent(in)  :: nxp
  integer*4, intent(in)  :: lonp1(nxp), latp1(nxp)
  integer*4, intent(in)  :: lonp2(nxp), latp2(nxp)
  real*8, intent(in)     :: w1(nxp), w2(nxp), w3(nxp), w4(nxp)
  real*8, intent(out)    :: gp_prof(nxp, nl)
  real*8, intent(in), optional::mask_val
  
  ! local variables
  
  integer:: ilon1, ilon2, ilat1, ilat2

  real*8::fmask
  integer::i, k
  real*8::total_w
  
  fmask=-999.0
  
  if (present(mask_val)) fmask=mask_val
  
  ! print *, 'fmask', fmask
  ! print *, 'maxval', maxval(gp)
  
  
  do i=1, nxp  ! points
     
     ilon1=lonp1(i)
     ilat1=latp1(i)
     
     ilon2=lonp2(i)
     ilat2=latp2(i)
     if ((ilon1==fmask).or.(ilon2==fmask).or.(ilat1==fmask).or.(ilat2==fmask)) then
        
        gp_prof(i,:)=fmask
     
     else if ((w1(i)==fmask).or. (w2(i)==fmask)) then
        gp_prof(i, :) =fmask
     else if ((w3(i)==fmask).or.(w4(i)==fmask)) then 
        gp_prof(i, :) =fmask
     else
        
        ilon1=lonp1(i)+1
        ilat1=latp1(i)+1
        
        ilon2=lonp2(i)+1
        ilat2=latp2(i)+1
        do k=1, nl ! vertical level 
           
           ! filter out the missing values
           gp_prof(i, k)=0.0
           total_w=0.0
           ! p1 
           if (gp(ilon1, ilat1,k)/=fmask) then
              gp_prof(i, k)=gp_prof(i, k)+ w1(i)*gp(ilon1, ilat1,k)
              total_w= total_w+ w1(i)
           end if
           
           ! p2
           if (gp(ilon1, ilat2,k)/=fmask) then
              gp_prof(i, k)=gp_prof(i, k)+ & 
                   w2(i)*gp(ilon1, ilat2,k)
              total_w= total_w+ w2(i)
           end if
           
           
           ! p3
           if (gp(ilon2, ilat1,k)/=fmask) then
              gp_prof(i, k)=gp_prof(i, k)+ & 
                   w3(i)*gp(ilon2, ilat1,k)
              total_w= total_w+w3(i)
           end if
           
           
           ! p4
           if (gp(ilon2, ilat2,k)/=fmask) then
              gp_prof(i, k)=gp_prof(i, k)+ & 
                   w4(i)*gp(ilon2, ilat2,k)
              total_w= total_w+ w4(i)
           end if
     !       print *, ilon1, ilon2, ilat1, ilat2, gp(ilon2, ilat2,k)
           
     !      print *, total_w, gp_prof(i,k)
           
           
           if (total_w>0.5) then
              gp_prof(i, k)=gp_prof(i, k)/total_w
           else 
              ! under heavy influences of bad points 
              gp_prof(i, k)=fmask
           end if
        end do
     end if
  end do
  
end subroutine sample_field


subroutine sample_field_em(nlat, nlon, nl, ne, gp, nxp, lonp1, &
     lonp2, latp1, latp2, w1, w2, w3, w4, gp_prof, mask_val)
  
  ! Get model profiles for multiple tracers 
  !
  ! Inputs
  !-------------------------------------------
  ! 1. integer*4::  nlat: size of latitude 
  ! 2. integer*4::  nlon: size of longitude 
  ! 3. integer*4::  nl:  vertical levels
  ! 4. integer*4::  ne :  number of  tracers (time)
  ! 5. real*8   :: gp(nlon, nlat, nl,ne): gridded model data
  ! 6 .integer*4::  nxp: Point number along the track
  ! 7,8. integer*4::  lonp1(nxp), lonp2(nxp):  left-right longitude indexes
  ! 8, 9, integer*4:: latp1(nxp), lonp2(nxp):  left-right latitude indexes
  ! 9, 10, 11, 12.  real*8 :: w1(nxp), w2(nxp), w3(nxp), w4(nxp): 
  ! ---weighting factors for four corners around each points. 
  ! 13. real*8:: mask_val: filling value for missing or bad data 
  ! 
  ! Outputs
  !--------------------------------------
  ! real*8, intent(out)    :: gp_prof(nxp, nl, ne) : profiles along slice
  !
  

implicit none

integer*4, intent(in)  :: nlat, nlon, nl, ne
real*8, intent(in)     :: gp(nlon, nlat, nl,ne)
integer, intent(in)    :: nxp
integer*4, intent(in)  :: lonp1(nxp), latp1(nxp), lonp2(nxp), latp2(nxp)
real*8, intent(in)     :: w1(nxp), w2(nxp), w3(nxp), w4(nxp)
real*8, intent(out)    :: gp_prof(nxp, nl,ne)
real*8, intent(in), optional::mask_val
! local variables

integer:: ilon1, ilon2, ilat1, ilat2

real*8::fmask
integer::i, j, k
real*8::total_w

fmask=-999.0

if (present(mask_val)) fmask=mask_val



do i=1, nxp  ! points
   ilon1=lonp1(i)
   ilat1=latp1(i)
   
   ilon2=lonp2(i)
   ilat2=latp2(i)
   

 
   if ((ilon1==fmask).or.(ilon2==fmask).or.(ilat1==fmask).or.(ilat2==fmask)) then
      gp_prof(i, :,:) =fmask
   else if ((w1(i)==fmask).or. (w2(i)==fmask)) then
      gp_prof(i, :,:) =fmask
   else if ((w3(i)==fmask).or.(w4(i)==fmask)) then 
      gp_prof(i, :,:) =fmask
   else
      ilon1=lonp1(i)+1
      ilat1=latp1(i)+1
      
      ilon2=lonp2(i)+1
      ilat2=latp2(i)+1
      
      do k=1, nl  ! level 
         
         do j=1,ne  ! tracer
            
            gp_prof(i, k,j)=0.0
            total_w=0.0
            ! p1 
            if (gp(ilon1, ilat1,k, j)/=fmask) then
               gp_prof(i, k, j)=gp_prof(i, k,j)+ & 
                    w1(i)*gp(ilon1, ilat1,k,j)
               total_w= total_w+ w1(i)
            end if
         
            ! p2
            
            if (gp(ilon1, ilat2,k,j)/=fmask) then
               gp_prof(i, k,j)=gp_prof(i, k,j)+ & 
                    w2(i)*gp(ilon1, ilat2,k,j)
               total_w= total_w+ w2(i)
            end if
            
         
            ! p3
         
            if (gp(ilon2, ilat1,k,j)/=fmask) then
               gp_prof(i, k,j)=gp_prof(i, k,j)+ & 
                    w3(i)*gp(ilon2, ilat1,k,j)
               total_w= total_w+w3(i)
            end if
            
         
            ! p4
            if (gp(ilon2, ilat2,k,j)/=fmask) then
               gp_prof(i, k,j)=gp_prof(i, k,j)+& 
                    w4(i)*gp(ilon2, ilat2,k,j)
               total_w= total_w+ w4(i)
            end if
            
            if (total_w>0.5) then
               gp_prof(i, k,j)=gp_prof(i, k,j)/total_w
            else
               ! under heavy influences of bad points
               gp_prof(i, k,j)=fmask
            end if
         
         end do ! j
      end do  ! k
   end if
end do   ! i


end subroutine sample_field_em




subroutine regrid_field(nlat, nlon, nl, gp, & 
     ob_nlon,  lonp1, lonp2, w1, & 
     ob_nlat, latp1, latp2, w3, & 
     ob_gp, mask_val)
  

  ! regrid model field 
  !
  ! Inputs
  !-------------------------------------------
  ! 1. integer*4::  nlat: size of latitude 
  ! 2. integer*4::  nlon: size of longitude 
  ! 3. integer*4::  nl:  vertical levels
  ! 4. real*8   :: gp(nlon, nlat, nl): gridded model data
  ! 5. integer*4:: ob_nlon: size of new longitude (x)
  ! 6,7. integer*4::  lonp1(ob_nlon), lonp2(ob_nlon):  left-right longitude indexes
  ! 8.  real*8::     w1(ob_nlon): weighting for left longitude boundary
  ! 9. integer*4::  ob_nlat: size of new latitude (y)
  ! 10, 11, integer*4:: latp1(ob_nlat), lat2(ob_nlat):  left-right latitude indexes
  ! 12.  real*8 :: w3(ob_nlat) :  weighting for left latitude boundary
  ! 13. real*8:: mask_val: filling value for missing or bad data 
  ! 
  ! Outputs
  !--------------------------------------
  ! real*8, intent(out)    :: gp_prof(ob_nlon, ob_nlat, nl) : regridded fields 
  ! 
  ! Notes:
  !------------------------------------------------------
  ! 1. if the impact of bad or missing model values larger than 50% to a track point, 
  ! ---model value at this point will be set to mask value 


implicit none

integer*4, intent(in)  :: nlat, nlon, nl 
real*8, intent(in)     :: gp(nlon, nlat, nl)

!  olon 
integer*4, intent(in)  :: ob_nlon
integer*4, intent(in)  :: lonp1(ob_nlon), lonp2(ob_nlon)
real*8, intent(in)     :: w1(ob_nlon)


! olat 
integer*4, intent(in)  :: ob_nlat
integer*4, intent(in)  :: latp1(ob_nlat), latp2(ob_nlat)
real*8, intent(in)     :: w3(ob_nlat)
real*8, intent(in), optional :: mask_val 
real*8, intent(out)    :: ob_gp(ob_nlon, ob_nlat, nl)

! local variables

integer:: ilon1, ilon2, ilat1, ilat2

real*8::fmask
integer::i, j, k
real*8::total_w

fmask=-999.0
if (present(mask_val)) fmask=mask_val


do i=1, ob_nlon  ! lon 
   
   ilon1=lonp1(i)
   ilon2=lonp2(i)
   do j=1, ob_nlat  ! lat 
      ilat1=latp1(j)
      ilat2=latp2(j)
      if ((ilon1==fmask).or.(ilon2==fmask).or.(ilat1==fmask).or.(ilat2==fmask)) then
           ob_gp(i, j,:) =fmask
      
      else if ((w1(i)==fmask).or.(w3(j)==fmask)) then ! check filling value 
         ob_gp(i, j, :) =fmask
      else
         ! shift to FORTRAN index
         
         ilon1=ilon1+1
         ilon2=ilon2+1
         ilat1=ilat1+1
         ilat2=ilat2+1
         
         do k=1, nl ! lz 
            ! filter out the missing values 
            ob_gp(i, j,k)=0.0
            total_w=0.0
            ! p1 
            if (gp(ilon1, ilat1,k)/=fmask) then
               ob_gp(i, j,k)=ob_gp(i, j,k)+ & 
                    w1(i)*w3(j)*gp(ilon1, ilat1,k)
               total_w= total_w+ w1(i)*w3(j)
            end if
            
            ! p2
            if (gp(ilon1, ilat2,k)/=fmask) then
               ob_gp(i, j,k)=ob_gp(i, j,k)+ & 
                    w1(i)*(1.0-w3(j))*gp(ilon1, ilat2,k)
               total_w= total_w+ w1(i)*(1.0-w3(j))
            end if
            
            
            ! p3
            
            if (gp(ilon2, ilat1,k)/=fmask) then
               ob_gp(i, j,k)=ob_gp(i, j,k)+ & 
                    (1.0-w1(i))*w3(j)*gp(ilon2, ilat1,k)
               total_w= total_w+ (1.0-w1(i))*w3(j)
            end if
            
            
            ! p4
            
            if (gp(ilon2, ilat2,k)/=fmask) then
               ob_gp(i, j,k)=ob_gp(i, j,k)+ & 
                    (1.0-w1(i))*(1.0-w3(j))*gp(ilon2, ilat2,k)
               total_w= total_w+ (1.0-w1(i))*(1.0-w3(j))
            end if
!            print *, i, j, k, total_w, ob_gp(i,j,k)
            
            if (total_w>0.5) then
               ob_gp(i, j,k)=ob_gp(i, j,k)/total_w
            else
               ! under heavy influences of bad points
               ob_gp(i, j,k)=fmask
            end if
         end do ! k
      end if
   end do  ! j 
end do  ! i


end subroutine regrid_field
