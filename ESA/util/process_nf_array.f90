subroutine process_nf_array(nlon, nlat, nlvl, ntime, data, &
     scaling, offset, filling_value, new_filling, outdata)
  

! post-processing cloud data read from ECMWF netcdf file 
! Real cloud coverage in percentage (0-1) is given by 
! data=scaling*data+offset 
! 
  
implicit none
! inputs
! -------------------
!! dimension 

integer, intent(in)::nlon, nlat, nlvl, ntime   
real*8, intent(in)::scaling ! scaling factor 
real*8, intent(in)::offset ! offset 
real*8, intent(in)::filling_value ! bad or missing value 
real*8, intent(in)::new_filling ! new bad or missing value 

! data
! --------------------------
real*8, intent(in)::data(nlon, nlat, nlvl, ntime)


! outputs 
! --------------------------
real*8, intent(out)::outdata(nlon, nlat, nlvl, ntime)
 
! local variable 
integer::i, j, l, k 

!$OMP DO PRIVATE(i,j,l, k)

do i=1, nlon ! lon 
   do j=1, nlat ! lat
      do l=1, nlvl
         do k=1, ntime ! time
            if (data(i, j, l, k)==filling_value) then
               outdata(i, j, l, k)=new_filling
            else
               ! only process valid data
               outdata(i,j,l, k)=scaling*data(i,j,l, k)+offset
            end if
         end do ! k
      end do ! l
   end do ! j
end do !i



!$OMP END DO

end subroutine process_nf_array



subroutine fill_nf_3d_int(nlon, nlat, ntime, data, &
     filling_value, new_filling)

! subroutine process_nf_3d(nlon, nlat, ntime, data, &
!     scaling, offset, filling_value, new_filling)
  

! post-processing cloud data read from ECMWF netcdf file 
! Real cloud coverage in percentage (0-1) is given by 
! data=scaling*data+offset 
! 
  
implicit none
! inputs
! -------------------
!! dimension 

integer, intent(in)::nlon, nlat, ntime   
integer, intent(in)::filling_value ! bad or missing value 
integer, intent(in)::new_filling ! new bad or missing value 

! data
! --------------------------
integer, intent(inout)::data(nlon, nlat, ntime)


!out
! real*8, intent(out)::outdata(nlon, nlat, ntime)

 
! local variable 
integer::i, j, k 

!$OMP DO PRIVATE(i,j,k)

do i=1, nlon ! lon 
   do j=1, nlat ! lat
      do k=1, ntime ! time
         if (data(i, j, k)==filling_value) then
            ! outdata(i, j, k)=new_filling
            data(i, j, k)=new_filling
         end if
      end do ! k
   end do ! j
end do !i



!$OMP END DO

end subroutine fill_nf_3d_int



subroutine process_nf_3d(nlon, nlat, ntime, data, &
     scaling, offset, filling_value, new_filling, outdata)

! subroutine process_nf_3d(nlon, nlat, ntime, data, &
!     scaling, offset, filling_value, new_filling)
  

! post-processing cloud data read from ECMWF netcdf file 
! Real cloud coverage in percentage (0-1) is given by 
! data=scaling*data+offset 
! 
  
implicit none
! inputs
! -------------------
!! dimension 

integer, intent(in)::nlon, nlat, ntime   
real*8, intent(in)::scaling ! scaling factor 
real*8, intent(in)::offset ! offset 
real*8, intent(in)::filling_value ! bad or missing value 
real*8, intent(in)::new_filling ! new bad or missing value 

! data
! --------------------------
real*8, intent(in)::data(nlon, nlat, ntime)


!out
real*8, intent(out)::outdata(nlon, nlat, ntime)

 
! local variable 
integer::i, j, k 

!$OMP DO PRIVATE(i,j,k)

do i=1, nlon ! lon 
   do j=1, nlat ! lat
         do k=1, ntime ! time
            if (data(i, j, k)==filling_value) then
               outdata(i, j, k)=new_filling
               ! data(i, j, k)=new_filling
            else
               ! only process valid data
               outdata(i,j,k)=scaling*data(i,j,k)+offset
               ! data(i,j,k)=scaling*data(i,j,k)+offset
               
            end if
         end do ! k
      end do ! j
end do !i



!$OMP END DO

end subroutine process_nf_3d
