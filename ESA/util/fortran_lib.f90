! functions called by python
! Interpolation 
! --------------------------------------
! 1.  getpos :  find the lower (pl) and upper (pr) boundary of value x in a array x0
! 2.  getwgt :  find the lower (pl) and upper (pr) boundary of value x in a array x0, 
!     and the weigting factor for interpolate. 
! 3.  get_vertical_wgt_1d: position and wgt for vertical interpolation over nx horizontal points
! 4.  get_vertical_wgt_1d_sel: 
!


subroutine getpos(x0, x, nx0, nx, pl, pr)
! find the lower (pl) and upper (pr) boundary of value x in a array x0
! Inputs
! ----------------------------------------------------
!  integer::nx0, nx # size of old and new grid 
!
!  real*8::x0(nx0)  # grid 
!  real*8::x(nx)   # new grid
!
!  Outputs
! -----------------------------------------
!  integer::pl(nx) # Location of lower boundary for each x 
!  integer::pl(nx) # Location of upper boundary for each x


implicit none 
! inputs
integer, intent(in)::nx, nx0  
real*8, intent(in)::x0(nx0), x(nx)

! outputs
integer*4, intent(out)::pl(nx), pr(nx)   

! local
integer::i, j, l

j=1
do i=1, nx
   ! looking for upper boundary
   
   do l=j, nx0
      if (x0(l)>=x(i)) exit
   end do ! l
   
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
end do  ! i 

end subroutine getpos

subroutine getwgt(x0, x, nx0_in, nx, pl, pr,wgt)

! find the lower (pl) and upper (pr) boundary of value x in a array x0, 
!  and the weigting factor for interpolate. 

! Inputs
! ----------------------------------------------------
!  integer::nx0, nx # size of old and new grid 
!
!  real*8::x0(nx0)  # grid 
!  real*8::x(nx)   # new grid
!
!  Outputs
! -----------------------------------------
!  integer::pl(nx) # Location of lower boundary for each x 
!  integer::pl(nx) # Location of upper boundary for each x
!  real*8::wgt(nx) # Weighting factor for lower boundary 
! 
   

implicit none 
integer, intent(in)::nx, nx0_in
real*8, intent(in)::x0(nx0_in), x(nx)
integer*4, intent(out)::pl(nx), pr(nx)
real*8, intent(out)::wgt(nx)
integer::i, j, l
integer::nx0, nx_st


pl(:)=-999
pr(:)=-999
wgt(:)=0.0

nx0=0

do i=1, nx0_in
   if (x0(i)<-990.0) exit
   nx0=nx0+1
end do

j=1



do i=1, nx
   if (x(i)<-990.0) exit
   
   do l=j, nx0
      if (x0(l)>=x(i)) exit
   end do
   
   if (l==1) then 
      pl(i)=0
      pr(i)=0
      wgt(i)=1.0
      
   else if (l>nx0) then  
      pl(i)=nx0-1
      pr(i)=nx0-1
      wgt(i)=1.0
      
   else
      pr(i)=l-1
      pl(i)=pr(i)-1
      wgt(i)=(x0(pr(i)+1)-x(i))/(x0(pr(i)+1)-x0(pl(i)+1))
   end if
   j=l
end do

end subroutine getwgt


subroutine get_vertical_wgt_1d(z0, z, nx, nz0, nz, pl, pr, wgt)
  ! position and wgt for vertical interpolation over nx horizontal points
  ! Inputs
  !------------------------------------------------------------
  ! integer::nx # horizontal point number
  ! integer::nz0 # size of the old vertical coordinate
  ! integer::nz # size of the new vertical coordinate
  
  ! real*8:: z0(nx,nz0) #  old vertical coordinate 
  ! real*8:: z(nx, nz) # new vertical coordiante
  !
  ! Outputs:
  !-------------------------------------------- 
  ! integer:: pl(nx, nz)  # lower boundary 
  ! integer::pr(nx, nz) # upper boundary
  ! real*8::wgt(nx, nz) # weighting factor at lower boundary
  ! 

  implicit none
  ! Inputs
  integer, intent(in)::nx, nz0, nz
  real*8, intent(in)::z0(nx, nz0)
  real*8, intent(in)::z(nx, nz)
  ! outputs
  integer*4, intent(out)::pl(nx, nz)
  integer*4, intent(out)::pr(nx, nz)
  real*8, intent(out)::wgt(nx,nz)
  ! locals
  integer::i
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  
  do i=1, nx
     call getwgt(z0(i, :), z(i, :), nz0, nz, pl1d, pr1d,wgt1d)
     pl(i,:)=pl1d(:)
     pr(i,:)=pr1d(:)
     wgt(i,:)=wgt1d(:)
  end do
  
  deallocate(pl1d, pr1d, wgt1d)
  
  
end subroutine get_vertical_wgt_1d

subroutine get_vertical_wgt_1d_sel(x, nx, nz, sel_prof, nsel, xout, nz_out, pl, pr, wgt)
  
  ! position and wgt for vertical interpolation over nx horizontal points
  ! 
  ! Inputs
  !------------------------------------------------------------
  ! integer::nprof # profile number
  ! integer::nz # size of the old vertical coordinate
  ! integer::nz_out # size of the new vertical coordinate
  ! integer::nsel # selected points
  ! real*8:: x(nx, nz) # old vertical coordinate
  !
  ! Outputs:
  !-------------------------------------------- 
  ! real*8:: xout(nsel, nz_out) # x values
  ! integer:: pl(nx, nz)  # lower boundary 
  ! integer::pr(nx, nz) # upper boundary
  ! real*8::wgt(nx, nz) # weighting factor at lower boundary
  ! 

! find the lower (pl) and upper (pr) boundary of value x in array x0, 
! and the weigting factor for interpolate. 
! x ------in--------- array of size (nx, nz) order.
! sel_prof ------in--------- selected position of the profiles. 
! nsel  --------in --------- the number of selected profiles
 
! xout--------out------------  array of (nsel, nz_out) 
! pr ------out-------- array of size (nsel, nz_out)  vectical location of the upper boundary
! pl ------out-------- array of size (nsel, nz_out)  vectical location of the lower boundary
! wgt ------out ------- array of size nx. The weighting factor 
! get the lower (pl) and upper (pr) boundary of value x0 in a array x

implicit none
real*8, intent(in)::x(nx, nz)
integer*4, intent(in)::sel_prof(nsel)
integer, intent(in)::nx, nz, nsel, nz_out
real*8, intent(in)::xout(nsel, nz_out)
integer*4, intent(out)::pl(nsel, nz_out)
integer*4, intent(out)::pr(nsel, nz_out)
real*8, intent(out)::wgt(nsel, nz_out)

real*8, allocatable::x0(:), x1(:), wgt1d(:)
integer, allocatable::pl1d(:), pr1d(:)

integer::i, iobs, j, k


allocate(x0(nz))
allocate(x1(nz_out),pl1d(nz_out),pr1d(nz_out),wgt1d(nz_out))

! print*, size(x(1,:)), size(x(:,1))
! print*, size(sel_prof(:))
! print*, x(1,1), x(1,2), x(1,3), xout(1,1), xout(1,2), xout(1,3)





do i=1, nsel
   iobs=sel_prof(i)+1
   x0=x(iobs, :)
   x1=xout(iobs,:)
   call getwgt(x0, x1, nz, nz_out, pl1d, pr1d,wgt1d)
   pl(i, :)=pl1d(:)
   pr(i, :)=pr1d(:)
   wgt(i, :)=wgt1d(:)
end do
! print *, 'pl', pl(1, 1:5)
! print *, 'pr', pr(1, 1:5)

deallocate(x0)
deallocate(x1, pl1d, pr1d, wgt1d)

end subroutine get_vertical_wgt_1d_sel


subroutine get_vertical_wgt_2d(z0, z, nx, ny, nz0, nz, pl, pr, wgt)
  ! position and wgt for vertical interpolation over a mesh  of (nx) X (jy)
  ! z0 -----in-----  old vertical coordinate 
  ! z ------in------ new vertical coordiante
  ! nx ------in------- size in x direction 
  ! ny ------in-------- size in y direction 
  ! nz0 ----in--------- size in z direction 
  ! nz -----in --------- size in z direction for out-put 
  ! pl ----out---------  lower boundary 
  ! pr ----out --------- higher boundary
  ! wgt ---out --------  weeighting at lower boundary
  
  implicit none
  integer, intent(in)::nx, ny, nz0, nz
  real*8, intent(in)::z0(nx, ny, nz0)
  real*8, intent(in)::z(nx, ny, nz)
  integer*4, intent(out)::pl(nx, ny, nz)
  integer*4, intent(out)::pr(nx, ny, nz)
  real*8, intent(out)::wgt(nx, ny, nz)
  ! local variables
  integer::i, j,k
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  do i=1, nx
     do j=1, ny
        call getwgt(z0(i, j, :), z(i, j,:), nz0, nz, pl1d, pr1d,wgt1d)
        pl(i,j, :)=pl1d(:)
        pr(i,j, :)=pr1d(:)
        wgt(i,j, :)=wgt1d(:)
     end do
  end do
  
  deallocate(pl1d, pr1d, wgt1d)
  
  
end subroutine get_vertical_wgt_2d


subroutine prof_vertical_intpl_1d(pl, pr, wgt, prof, out_prof, nx, nz,nz0)
  ! vertical interpolation of profiles over x direction
  ! pl ----in---------  lower boundary 
  ! pr ----in --------- higher boundary
  ! wgt ---in --------  weighting at lower boundary
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x direction
  
 
  implicit none
  integer, intent(in)::nx, nz
  integer, intent(in)::nz0
  
  real*8, intent(in)::prof(nx, nz0)
  real*8, intent(in)::wgt(nx, nz)
  integer*4, intent(in)::pl(nx, nz), pr(nx, nz)
  real*8, intent(out)::out_prof(nx, nz)
  
  integer::i, k, kl,kr 
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  
  do i=1, nx
     do k=1, nz
        kl=pl(i,k)+1
        kr=pr(i,k)+1
        ! print *, 'k', kl, kr
        if ((kl>-990).and.(kr>-990)) then 
           out_prof(i, k)=wgt(i, k)*prof(i, kl)+(1.0-wgt(i,k))*prof(i, kr)
           if ((prof(i, kl)<-990.0).or.(prof(i, kr)<-990.0)) out_prof(i, k)=-999.0  
        else
           out_prof(i, k)=-999.0
        end if
     end do
  end do
end subroutine prof_vertical_intpl_1d

subroutine prof_vertical_intpl_em(pl, pr, wgt, wv,prof, out_prof, nx, ny, nz,nz0)
  
  ! vertical interpolation of profiles over mesh of (nx) X (ny)
  ! pl ----in---------  lower boundary 
  ! pr ----in --------- higher boundary
  ! wgt ---in --------  weighting at lower boundary
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x and y direction
  
  
  implicit none
  integer, intent(in)::nx, ny, nz
  integer, intent(in)::nz0
  
  real*8, intent(in)::prof(nx, ny, nz0)
  real*8, intent(in)::wgt(ny, nz), wv(ny, nz)
  integer*4, intent(in)::pl(ny, nz), pr(ny, nz)
  real*8, intent(out)::out_prof(nx, ny, nz)
  
  integer::i, j, k, kl,kr 
  real*8::scaled_water_r
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  do j=1, ny
     do k=1, nz
        kl=pl(j, k)+1
        kr=pr(j, k)+1
        scaled_water_r=wv(j,k)
        ! print *, 'k', kl, kr
        if ((kl>-990).and.(kr>-990).and.(scaled_water_r>=0.0)) then 
           do i=1, nx
              out_prof(i, j,k)=wgt(j,k)*prof(i, j, kl)+(1.0-wgt(j,k))*prof(i,j, kr)
              out_prof(i, j,k)=out_prof(i, j,k)/(1.0+scaled_water_r)
              ! put the filled value if original data are invalid 
              if ((prof(i, j, kl)<-990.0).or.(prof(i, j, kr)<-990.0)) out_prof(i, j,k)=-999.0  
           end do
        else
           out_prof(:, j, k)=-999.0
        end if
     end do
  end do
  
  
end subroutine prof_vertical_intpl_em


subroutine prof_vertical_intpl_wv(pl, pr, wgt, wv,prof, out_prof, ny, nz,nz0)
  
  ! vertical interpolation of profiles over mesh of (nx) X (ny)
  ! pl ----in---------  lower boundary 
  ! pr ----in --------- higher boundary
  ! wgt ---in --------  weighting at lower boundary
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x and y direction
  
  
  implicit none
  integer, intent(in)::ny, nz
  integer, intent(in)::nz0
  
  real*8, intent(in)::prof(ny, nz0)
  real*8, intent(in)::wgt(ny, nz), wv(ny, nz)
  integer*4, intent(in)::pl(ny, nz), pr(ny, nz)
  real*8, intent(out)::out_prof(ny, nz)
  
  integer::j, k, kl,kr 
  real*8::scaled_water_r
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  do j=1, ny
     do k=1, nz
        kl=pl(j, k)+1
        kr=pr(j, k)+1
        scaled_water_r=wv(j,k)
        ! print *, 'k', kl, kr
        if ((kl>-990).and.(kr>-990).and.(scaled_water_r>=0.0)) then 
           
           out_prof(j,k)=wgt(j,k)*prof(j, kl)+(1.0-wgt(j,k))*prof(j, kr)
           out_prof(j,k)=out_prof(j,k)/(1.0+scaled_water_r)
           ! put the filled value if original data are invalid 
           if ((prof(j, kl)<-990.0).or.(prof(j, kr)<-990.0)) out_prof(j,k)=-999.0  
           
        else
           out_prof(j, k)=-999.0
        end if
     end do
  end do
  
  
end subroutine prof_vertical_intpl_wv

subroutine ak_col_int_em(colwgt, ak, prof, nx, ny, nz, colval)
  ! calculate the column values for a list of profiles using the weighting factors and average kernel.
  ! prof -------in-------- array of (nx, nz) volume mixing ratios
  ! colwagt ------in ------- array of (nx, nz) weighting factor for column integration 
  !  ak --------in ------------ array of (nx, nz) the averaging kernel . 
  !
  ! colval -------out ------- the column values at nx locations. 
  

  implicit none
  integer, intent(in)::nx, ny,nz
  real*8, intent(in)::prof(nx, ny, nz)
  real*8, intent(in)::ak(ny, nz)
  real*8, intent(in)::colwgt(ny, nz)
  real*8, intent(out)::colval(nx, ny)
  integer::i, j, k
  
  do i=1, nx
     do j=1, ny
        colval(i,j)=0.0
        do k=1, nz
           if ((prof(i, j, k)>-990.0).and.(ak(j,k)>-990.0)) then 
              colval(i,j)=colval(i,j)+colwgt(j,k)*ak(j,k)*prof(i, j, k)
           end if
        end do
     end do
  end do
  
  
end subroutine ak_col_int_em



subroutine prof_vertical_intpl_2d(pl, pr, wgt, prof, out_prof, nx, ny, nz,nz0)
  
  ! vertical interpolation of profiles over mesh of (nx) X (ny)
  ! pl ----in---------  lower boundary 
  ! pr ----in --------- higher boundary
  ! wgt ---in --------  weighting at lower boundary
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x and y direction
  
  
  implicit none
  integer, intent(in)::nx, ny, nz
  integer, intent(in)::nz0
  
  real*8, intent(in)::prof(nx, ny, nz0)
  real*8, intent(in)::wgt(nx, ny, nz)
  integer*4, intent(in)::pl(nx, ny, nz), pr(nx, ny, nz)
  real*8, intent(out)::out_prof(nx, ny, nz)
  
  integer::i, j, k, kl,kr 
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
  do i=1, nx
     do j=1, ny
        do k=1, nz
           kl=pl(i,j, k)+1
           kr=pr(i,j, k)+1
           ! print *, 'k', kl, kr
           if ((kl>-990).and.(kr>-990)) then 
              out_prof(i, j,k)=wgt(i,j,k)*prof(i, j, kl)+(1.0-wgt(i,j,k))*prof(i,j, kr)
              ! put the filled value if original data are invalid 
              if ((prof(i, j, kl)<-990.0).or.(prof(i, j, kr)<-990.0)) out_prof(i, j,k)=-999.0  
           else
              out_prof(i, j, k)=-999.0
           end if
        end do
     end do
  end do

end subroutine prof_vertical_intpl_2d


subroutine prof_horizontal_intpl_2d(xpl, xpr, xwgt, ypl, ypr, ywgt, prof, out_prof, nx, ny, nx0, ny0, nz)
  
  ! interpolating profiles along x and y direction
  ! xpl ----in---------  lower boundary along x  
  ! xpr ----in --------- higher boundary
  ! xwgt ---in --------  weighting at x lower boundary
  
  ! ypl ----in---------  lower boundary along y  
  ! ypr ----in --------- higher boundary along y
  ! ywgt ---in --------  weighting at y lower boundary
  
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x and y direction
  
  
  implicit none
  integer, intent(in)::nx, ny, nx0, ny0, nz
  
  
  real*8, intent(in)::prof(nx0, ny0, nz)
  real*8, intent(in)::xwgt(nx)
  integer*4, intent(in)::xpl(nx), xpr(nx)
  
  real*8, intent(in)::ywgt(ny)
  integer*4, intent(in)::ypl(ny), ypr(ny)
  
  real*8, intent(out)::out_prof(nx, ny, nz)
  
  integer::i, j, k, nv
  real*8::xw, yw
  integer::xl, xr, yl, yr
  real*8::wgt(4)
  real*8::v(4)
  real*8::avg_v, avg_w
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
             
  do i=1, nx
     
     xl=xpl(i)+1
     xr=xpr(i)+1
     xw=xwgt(i)
     
     do j=1, ny
        yl=ypl(j)+1
        yr=ypr(j)+1
        yw=ywgt(j)
        
        do k=1, nz
           ! point 1
           v(1)=prof(xl,yl, k)
           wgt(1)=xw*yw
           ! point 2
           v(2)=prof(xr,yl, k)
           wgt(2)=(1.0-xw)*yw
           ! point 3
           v(3)=prof(xl,yr, k)
           wgt(3)=xw*(1.0-yw)
           ! point 4
           v(4)=prof(xr,yr, k)
           wgt(4)=(1.0-xw)*(1.0-yw)
           
           avg_v=0.0
           avg_w=0.0
           
           do nv=1, 4
              if (v(nv)>-990.0) then 
                 avg_v=avg_v+wgt(nv)*v(nv)
                 avg_w=avg_w+wgt(nv)
              end if
           end do
           if (avg_w>0.0) then 
              out_prof(i, j, k)=avg_w/avg_v
           else
              out_prof(i, j, k)=-999.0
           end if
        
        end do
     end do
  end do
  
end subroutine prof_horizontal_intpl_2d


subroutine prof_horizontal_intpl_1d(xpl, xpr, xwgt, ypl, ypr, ywgt, prof, & 
     out_prof, nx, nx0, ny0, nz)
  
  ! interpolating profiles along x and y direction
  ! xpl ----in---------  lower boundary along x  
  ! xpr ----in --------- higher boundary
  ! xwgt ---in --------  weighting at x lower boundary
  
  ! ypl ----in---------  lower boundary along y  
  ! ypr ----in --------- higher boundary along y
  ! ywgt ---in --------  weighting at y lower boundary
  ! nx ---- in -------   the number of profiles 
  ! prof -------in ----- the profile  along x direction 
  ! out_prof -------out ----- interpolated profiles along x and y direction
  
  
  implicit none
  integer, intent(in)::nx, nx0,ny0, nz
  
  
  real*8, intent(in)::prof(nx0, ny0, nz)
  real*8, intent(in)::xwgt(nx)
  integer*4, intent(in)::xpl(nx), xpr(nx)
  
  real*8, intent(in)::ywgt(nx)
  integer*4, intent(in)::ypl(nx), ypr(nx)
  
  real*8, intent(out)::out_prof(nx, nz)
  
  integer::i, j, k, nv
  real*8::xw, yw
  integer::xl, xr, yl, yr
  real*8::wgt(4)
  real*8::v(4)
  real*8::avg_v, avg_w
  
  ! print *, pl(1,:)
  ! print *, pr(1,:)
  
             
  do i=1, nx
     
     xl=xpl(i)+1
     xr=xpr(i)+1
     xw=xwgt(i)
     
     yl=ypl(j)+1
     yr=ypr(j)+1
     yw=ywgt(j)
     
     do k=1, nz
        ! point 1
        v(1)=prof(xl,yl, k)
        wgt(1)=xw*yw
        ! point 2
        v(2)=prof(xr,yl, k)
        wgt(2)=(1.0-xw)*yw
        ! point 3
        v(3)=prof(xl,yr, k)
        wgt(3)=xw*(1.0-yw)
        ! point 4
        v(4)=prof(xr,yr, k)
        wgt(4)=(1.0-xw)*(1.0-yw)
        
        avg_v=0.0
        avg_w=0.0
        
        do nv=1, 4
           if (v(nv)>-990.0) then 
              avg_v=avg_v+wgt(nv)*v(nv)
              avg_w=avg_w+wgt(nv)
           end if
        end do
        if (avg_w>0.0) then 
           out_prof(i, k)=avg_w/avg_v
        else
           out_prof(i, k)=-999.0
        end if
        
        
     end do
  end do
  
end subroutine prof_horizontal_intpl_1d


subroutine get_col_wgt(pres_in, nz, wgt)
! the weighting factor for calculation of the total column from vertical profiles, 
! pres_in  ------in--------- array of nz, pressure in hPa or log10(pressure) in descending order
! wgt ----------out --------- array of nz, the weighting factor for column value calculation.

implicit none

integer, intent(in)::nz
real*8, intent(in) ::pres_in(nz)
real*8, intent(out)::wgt(nz)
real*8, allocatable::pres_a(:), pres(:), fsum(:), logp(:)
real*8  :: a=log(10.0)
real*8::w1, w2
integer::i, nz_real


wgt(:)=0.0
nz_real=0
! 
do i=1, nz
   if (pres_in(i)<-990.0) exit
   nz_real=nz_real+1
end do

allocate(pres_a(nz_real), logp(nz_real), pres(nz_real))
allocate(fsum(nz_real))

pres(1:nz_real)=pres_in(1:nz_real)
! guess it is log10(pressure) or pressure in hPa

if (maxval(pres)>10) then 
   pres_a=pres
   logp=log10(pres)
else
   pres_a=10.0**(pres)
   logp=pres
end if

fsum=0.0
! log to log 10  
a=log(10.0) 
! weighting factor for single layer 

do i=2, nz_real
   w1=pres_a(i)-pres_a(i-1)
   w2=a*(logp(i)-logp(i-1))
   
   if (w2==0) then 
      ! remove this level
      fsum(i-1)=0.0
   else
      fsum(i-1)=w1/w2
   end if
end do

do i=2, nz_real-1
   wgt(i)=fsum(i)-fsum(i-1)
end do

if (logp(2)==logp(1)) then 
   ! the starting point moving to second 
   wgt(1)=0.0
   wgt(2)=-pres_a(2)+(pres_a(3)-pres_a(2))/(a*(logp(3)-logp(2)))
else
   wgt(1)=-pres_a(1)+(pres_a(2)-pres_a(1))/(a*(logp(2)-logp(1)))
end if
   
i=nz_real
if (logp(i-1)==logp(i)) then
   ! the ending point move to next 
   wgt(i)=0.0
   i=nz_real-1
   wgt(i)=pres_a(i)-(pres_a(i-1)-pres_a(i))/(a*(logp(i-1)-logp(i)))
   
else 
   wgt(i)=pres_a(i)-(pres_a(i-1)-pres_a(i))/(a*(logp(i-1)-logp(i)))
end if

wgt(1:nz_real)=wgt(1:nz_real)/(pres_a(nz_real)-pres_a(1))

deallocate(fsum)
deallocate(pres_a, logp, pres)

end subroutine get_col_wgt


subroutine get_col_wgt_1d(pres, nx, nz, wgt)
! the weighting factor for calculation of the total column from vertical profil, 
! pres   ------in--------- array of size (nx, nz), pressure in hPa or log10(pressure) in descending order
! wgt ----------out --------- array of (nselm nz), the weighting factor for column value calculation.


implicit none

integer, intent(in)::nz, nx
real*8, intent(in) ::pres(nx, nz)

real*8, intent(out)::wgt(nx, nz)
! local variables
real*8, allocatable::cur_pres(:), cur_wgt(:)
integer::i
allocate(cur_pres(nz), cur_wgt(nz))

do i=1, nx
   cur_pres(:)=pres(i,:)
   call get_col_wgt(cur_pres, nz, cur_wgt)
   wgt(i,:)=cur_wgt(:)
end do

deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_1d


subroutine get_col_wgt_1d_bt(pres, nx, nz, wgt)
! the weighting factor for calculation of the total column from vertical profil, 
! pres   ------in--------- array of size (nx, nz), pressure in hPa or log10(pressure) in descending order
! wgt ----------out --------- array of (nselm nz), the weighting factor for column value calculation.


implicit none

integer, intent(in)::nz, nx
real*8, intent(in) ::pres(nx, nz)

real*8, intent(out)::wgt(nx, nz)
! local variables
real*8, allocatable::cur_pres(:), cur_wgt(:)
integer::i
allocate(cur_pres(nz), cur_wgt(nz))

do i=1, nx
   cur_pres(:)=pres(i,:)
   call get_col_wgt(cur_pres, nz, cur_wgt)
   wgt(i,:)=cur_wgt(:)
end do

deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_1d_bt





subroutine get_col_wgt_1d_sel(pres, nx, nz, sel_prof, nsel, wgt)
! the weighting factor for calculation of the total column from vertical profil, 
! pres   ------in--------- array of size (nx, nz), pressure in descending order
! sel_prof -----in-------- selection of profiles
! wgt ----------out --------- array of (nselm nz), the weighting factor for column value calculation.


implicit none

integer, intent(in)::nz, nx, nsel
real*8, intent(in) ::pres(nx, nz)
integer*4, intent(in)::sel_prof(nsel)

real*8, intent(out)::wgt(nsel, nz)
! local variables
real*8, allocatable::cur_pres(:), cur_wgt(:)
integer::iobs, i
allocate(cur_pres(nz), cur_wgt(nz))

do i=1, nsel
   ! python to fortran +1 
   iobs=sel_prof(i)+1 
   cur_pres(:)=pres(iobs,:)
   call get_col_wgt(cur_pres, nz, cur_wgt)
   wgt(i,:)=cur_wgt(:)
end do

deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_1d_sel

subroutine get_col_wgt_2d(pres, nx, ny,nz, wgt)
! the weighting factor for calculation of the total column from vertical profil, 
! pres   ------in--------- array of size (nx, nz), pressure in hPa or log10(pressure) in descending order
! wgt ----------out --------- array of (nselm nz), the weighting factor for column value calculation.


implicit none

integer, intent(in)::nz, nx, ny
real*8, intent(in) ::pres(nx, ny, nz)

real*8, intent(out)::wgt(nx, ny, nz)
! local variables
real*8, allocatable::cur_pres(:), cur_wgt(:)
integer::i, j
allocate(cur_pres(nz), cur_wgt(nz))

do i=1, nx
   do j=1, ny
      cur_pres(:)=pres(i,j, :)
      call get_col_wgt(cur_pres, nz, cur_wgt)
      wgt(i,j, :)=cur_wgt(:)
   end do
end do

deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_2d



subroutine get_col_wgt_sel(pres, nx, nz, sel_prof, nsel, wgt)

 ! the weighting factor for calculation of the total column from vertical profil, 
! pres-------in--------- array of (nx, nz), pressure in descending order
! sel_prof-----in-------- selection of profiles
! wgt ----------out --------- array of (nselm nz), the weighting factor for column value calculation.


implicit none

integer, intent(in)::nz, nx, nsel
real*8, intent(in) ::pres(nx, nz)
integer, intent(in)::sel_prof(nsel)

real*8, intent(out)::wgt(nsel, nz)
! local variables
real*8, allocatable::cur_pres(:), cur_wgt(:)
integer::iobs, i
allocate(cur_pres(nz), cur_wgt(nz))

do i=1, nsel
   ! python to fortran +1 
   iobs=sel_prof(i)+1 
   cur_pres(:)=pres(iobs,:)
   call get_col_wgt(cur_pres, nz, cur_wgt)
   wgt(i,:)=cur_wgt(:)
end do

deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_sel


subroutine col_int_1d(prof, colwgt, nx, nz, colval)
  
  ! calculate the column values for list of profiles using the weighting factor
  ! prof   -------in-------- array of (nx, nz)==volume mixing ratios
  ! colwagt ------in ------- array of (nx, nz)== weighting factor for column integration 
  ! colval -------out ------- array of (nx) == column mixing values at nx locations. 
    
  implicit none
  integer, intent(in)::nx, nz
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::colwgt(nx, nz)
  real*8, intent(out)::colval(nx)
  ! local variable
  
  integer::i, k
  
  do i=1, nx
     colval(i)=0
     do k=1, nz
        colval(i)=colval(i)+colwgt(i,k)*prof(i, k)
     end do
  end do
end subroutine col_int_1d


subroutine col_int_2d(prof, colwgt, nx, ny, nz, colval)
  
  ! calculate the column values for list of profiles using the weighting factor
  ! prof    -------in--------array of (nx, ny, nz)== the volume mixing ratios
  ! colwagt -------in --------array of (nx, ny, nz) ==weighting factor for column integration 
  ! colval  -------out -------array of (nx, ny)== the column mixing ratio values at mesh of (nx) X (ny). 
  
  
  implicit none
  integer, intent(in)::nx, ny, nz
  real*8, intent(in)::prof(nx, ny, nz)
  real*8, intent(in)::colwgt(nx, ny, nz)
  real*8, intent(out)::colval(nx, ny)
  integer::i, j, k
  
  do i=1, nx
     do j=1, ny
        colval(i,j)=0
        do k=1, nz
           colval(i,j)=colval(i,j)+colwgt(i,j,k)*prof(i,j, k)
        end do
     end do
  end do

end subroutine col_int_2d


subroutine ak_col_int_1d(colwgt, ak, prof, nx, nz, colval)
  ! calculate the column values for a list of profiles using the weighting factors and average kernel.
  ! prof -------in-------- array of (nx, nz) volume mixing ratios
  ! colwagt ------in ------- array of (nx, nz) weighting factor for column integration 
  !  ak --------in ------------ array of (nx, nz) the averaging kernel . 
  !
  ! colval -------out ------- the column values at nx locations. 
  

  implicit none
  integer, intent(in)::nx, nz
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::ak(nx, nz)
  real*8, intent(in)::colwgt(nx, nz)
  real*8, intent(out)::colval(nx)
  integer::i, j, k
  
  do i=1, nx
     colval(i)=0
     do k=1, nz
        colval(i)=colval(i)+colwgt(i,k)*ak(i,k)*prof(i, k)
     end do
  end do
end subroutine ak_col_int_1d


subroutine linear_itpl(x0, y0, x, nx0_in, nx, pl, pr,wgt,y)
! linear interpolate the values  
! x0 -------in------ array of (nx0_in)
! y0 ------in -------- array of nx0_in. Y values  
! x  ---- in --------- array of nx. The x values for interpolation
! y -------out ------- array of nx. The y values a fter interpolation 
! pl ------out ----------- array of nx. The lower boundary of the x in x0
! pr ------out ----------- array of nx. The higher boundary of the x in x0
! x0 & x in asending oder
  

 
implicit none 
integer, intent(in)::nx, nx0_in
real*8, intent(in)::x0(nx0_in), y0(nx0_in), x(nx)
integer*4, intent(out)::pl(nx), pr(nx)
real*8, intent(out)::wgt(nx), y(nx)
integer::i, j, l,nx0


j=1

y=0.0
wgt=0.0
pl(:)=-999
pr(:)=-999

nx0=0

do i=1, nx0_in
   if (x0(i)<-990.0) exit
   nx0=nx0+1
end do




do i=1, nx
   if (x(i)<-990.0) exit
   
   do l=j, nx0
      if (x0(l)>=x(i)) exit
      
   end do

   
   if (l==1) then 
      pl(i)=0
      pr(i)=0
      wgt(i)=1.0
      y(i)=y0(pl(i)+1)
      
      
   else if (l>nx0) then  
      pl(i)=nx0-1
      pr(i)=nx0-1
      wgt(i)=1.0
      y(i)=y0(pl(i)+1)
      
   else
      pr(i)=l-1
      pl(i)=pr(i)-1
      wgt(i)=(x0(pr(i)+1)-x(i))/(x0(pr(i)+1)-x0(pl(i)+1))
      y(i)=wgt(i)*y0(pl(i)+1)+(1.0-wgt(i))*y0(pr(i)+1)
      
   end if
   j=l
end do
end subroutine linear_itpl


subroutine reverse_matrix(m, nx, ny, mout, idx)
  implicit none

  integer, intent(in)::nx, ny
  integer, intent(in)::idx
  real*8, intent(in)::m(nx, ny)
  real*8, intent(out)::mout(nx, ny)

  integer::i
  

  if (idx==0) then
     do i=1, nx
        mout(i,:)=m(nx-i+1,:)
     end do
  else
     do i=1, ny
        mout(:,i)=m(:,ny-i+1)
     end do
  end if
  

end subroutine reverse_matrix


subroutine get_zonal_mean(y, pres, out_pres, nx, ny, nz0, nz, zm)
  ! calculate the zonal mean for the 3D field given at (lon, lat, pressure)
  ! pres -------in --------- pressure field in hPa and descending order
  ! nx ------in ---------  x size
  ! ny ------in ---------  y size
  ! nz0 ---- in ---------- z size
  ! nz -------in --------  z size for output 
  ! out_pres -------in ---------- pressure grid for output 
  ! 
  implicit none
  
  integer, intent(in)::nx, ny, nz0, nz
  real*8, intent(in)::y(nx, ny, nz0)
  real*8, intent(in)::pres(nx, ny, nz0)
  
  real*8, intent(in)::out_pres(nz)
  
  real*8, intent(out)::zm(ny, nz)
  
  ! local variables
  real*8, allocatable::prof(:),  logp0(:), logp(:), wgt(:)
  integer*4, allocatable::count_lvl(:, :), pl(:), pr(:)
  integer::i, j, k
  zm=0.0
  
  ! we need first interpolate and then do average if necessary 
  allocate(prof(nz), logp(nz), logp0(nz0), pl(nz), pr(nz), wgt(nz), count_lvl(ny, nz))
  count_lvl=0
  logp=-log10(out_pres)
  
  do i=1, nx
     do j=1, ny
        logp0=-log10(pres(i,j,:))
        call linear_itpl(logp0, y(i,j,:), logp, nz0, nz, pl, pr,wgt,prof)
        do k=1, nz
           if ((pres(i, j, 1)>=0.99*out_pres(k)).and.(prof(k)>-990.0)) then
              count_lvl(j, k)=count_lvl(j, k)+1
              zm(j, k)=zm(j, k)+prof(k)
           end if
        end do
     end do
  end do
  
  ! do average 
  
  do j=1, ny
     do k=1, nz
        if (count_lvl(j, k)>0) then
           zm(j, k)=zm(j, k)/count_lvl(j,k)
        else
           zm(j,k)=-999.0
        end if
     end do
  end do

end subroutine get_zonal_mean

subroutine do_svd_2(a, m,n, u, w, v, real_nx_ny, info)

  implicit none
  integer, intent(in)::m,n 
  real*8, intent(in)::a(m,n)
  real*8, intent(out)::u(m, m)
  integer, intent(in)::real_nx_ny
  
  real*8, intent(out)::v(real_nx_ny, n)

  real*8, intent(out)::w(m+n)
  integer, intent(out)::info
  ! local variable 
  CHARACTER(len=1)::          JOBU, JOBVT
  INTEGER::            LDA, LDU, LDVT, LWORK
  real*8, allocatable::a_tmp(:,:), WORK(:)
  info=-1
  ! allocate(a_tmp(m,n))
  ! a_tmp=a
  w=0.0
  JOBU='A'
  JOBVT='S'
  LDA=m
  LDU=m
  LDVT=real_nx_ny
  
  LWORK = MAX(1,3*MIN(m,n)+MAX(m,n),5*MIN(m,n))
  allocate(WORK(LWORK))
  

  call  DGESVD( JOBU, JOBVT, m, n, a, LDA, w, u, LDU, v, LDVT,&
       WORK, LWORK, info)
  
  ! deallocate(a_tmp)
  deallocate(WORK)
end subroutine do_svd_2

subroutine do_svd(a, m,n, u, w, v, info)

  implicit none
  integer, intent(in)::m,n 
  real*8, intent(in)::a(m,n)
  real*8, intent(out)::u(m, m)
  real*8, intent(out)::v(n, n)
  real*8, intent(out)::w(m+n)
  integer, intent(out)::info
  ! local variable 
  CHARACTER(len=1)::          JOBU, JOBVT
  INTEGER::            LDA, LDU, LDVT, LWORK
  real*8, allocatable::a_tmp(:,:), WORK(:)
  info=-1
  allocate(a_tmp(m,n))
  a_tmp=a
  w=0.0
  JOBU='A'
  JOBVT='A'
  LDA=m
  LDU=m
  LDVT=n
  LWORK = MAX(1,3*MIN(m,n)+MAX(m,n),5*MIN(m,n))
  allocate(WORK(LWORK))
  

  call  DGESVD( JOBU, JOBVT, m, n, a_tmp, LDA, w, u, LDU, v, LDVT,&
       WORK, LWORK, info)
  
  deallocate(a_tmp)
  deallocate(WORK)
end subroutine do_svd

subroutine do_ltrm_inv(a, n, b)
! inverse lower-triangle matrix 

  implicit none
  integer, intent(in)::n
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::b(n, n)
  real*8::sum_prev
  integer::ii, jj, kk
  ! last row 
  b(n,1:n)=0.0
  b(n, n)=1.0/a(n,n)
  do ii=n-1, 1, -1
     do jj=1, n
        sum_prev=sum(-a(ii, ii+1:n)*b(ii+1:n, jj))
        b(ii,jj)=sum_prev/a(ii,ii)
        if (abs(b(ii,jj))<1.0e-10) b(ii,jj)=0.0
        
     end do
     b(ii,ii)=b(ii,ii)+1.0/a(ii,ii)
  end do

end subroutine do_ltrm_inv


subroutine do_sym_eig(a, n,w, u, info)

  implicit none
  integer, intent(in)::n 
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::w(n)
  real*8, intent(out)::u(n, n)
  integer, intent(out)::info
  ! local variable 
  CHARACTER(len=1)::   JOBZ, UPLO
  INTEGER::            LDA, LDU, LDVT, LWORK
  real*8, allocatable::a_tmp(:,:), WORK(:)
  info=-1
  allocate(a_tmp(n,n))
  a_tmp=a
  w=0.0
  JOBZ='V'
  UPLO='U'
  
  LDA=n
  LWORK = max(1,3*n-1)
  allocate(WORK(LWORK))
  
  
  call  DSYEV( JOBZ, UPLO, n, a, LDA, w, WORK, LWORK, INFO )
  u=a
  deallocate(a_tmp)
  deallocate(WORK)
  
end subroutine do_sym_eig

subroutine do_sys_sqrt(a,n, sqrt_a)
  integer, intent(in)::n
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::sqrt_a(n,n)
  ! local variable
  real*8, allocatable::w(:), u(:,:), ah(:,:)
  integer:: info, ii, jj, kk
  allocate(w(n), u(n,n), ah(n,n))
  
  ah=a
  call do_sym_eig(ah, n,w, u, info)
  if (info.ne.0) then 
     print *,'error in do_sym_eig', info
     stop
  end if
  w=sqrt(w)
  ! print *, w
  ah=0.0
  
  do ii=1, n
     do jj=1,n
        do kk=1, n
           ah(ii, jj)=ah(ii, jj)+u(ii,kk)*w(kk)*u(jj, kk)
        end do
     end do
  end do
  sqrt_a=ah
  deallocate(w)
  deallocate(u)
  deallocate(ah)
end subroutine do_sys_sqrt

subroutine do_sys_sqrt_inv(a,n, sqrt_a)
  integer, intent(in)::n
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::sqrt_a(n,n)
  ! local variable
  real*8, allocatable::w(:), u(:,:), ah(:,:)
  integer:: info, ii, jj, kk
  allocate(w(n), u(n,n), ah(n,n))
  
  ah=a
  call do_sym_eig(ah, n,w, u, info)
  if (info.ne.0) then 
     print *,'error in do_sym_eig', info
     stop
  end if
  w=sqrt(1.0/w)
  ! print *, w
  ah=0.0
  
  do ii=1, n
     do jj=1,n
        do kk=1, n
           ah(ii, jj)=ah(ii, jj)+u(ii,kk)*w(kk)*u(jj, kk)
        end do
     end do
  end do
  sqrt_a=ah
  deallocate(w)
  deallocate(u)
  deallocate(ah)
end subroutine do_sys_sqrt_inv

subroutine do_svd_sqrt(a,n, sqrt_a)
  integer, intent(in)::n
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::sqrt_a(n,n)
  ! local variable
  real*8, allocatable::w(:), u(:,:), v(:,:), ah(:,:)
  integer:: info, ii, jj, kk
  allocate(w(n+n), u(n,n), v(n,n),ah(n,n))
  ah=0.0
  w=0.0
  ah=a
  
!  call do_sym_eig(a, n,w, u, info)
  call do_svd(ah, n, n, u, w, v, info)
  if (info.ne.0) then 
     print *,'error in do_sym_eig', info
     stop
  end if
  ! w=1.0/sqrt(w)
  ! print *, w
  ah=0.0
  
  do ii=1, n
     do jj=1,n
        do kk=1, n
           ah(ii, jj)=ah(ii, jj)+u(ii,kk)*v(jj, kk)/sqrt(w(kk))
        end do
     end do
  end do
  sqrt_a=ah
  deallocate(w)
  deallocate(u)
  deallocate(v)
  deallocate(ah)
end subroutine do_svd_sqrt

subroutine do_sym_inv_sqrt(a,n, sqrt_inv_a, info)
  integer, intent(in)::n
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::sqrt_inv_a(n,n)
  integer, intent(out)::info
  
  ! local variable
  real*8, allocatable::w(:), u(:,:), tmp1(:), tmp2(:)
  
  integer:: ii, jj, kk

  
  allocate(w(2*n), u(n,n), tmp1(n), tmp2(n))
  
  ! call do_svd(a, m,n, u, w, v, info)

  call do_sym_eig(a, n,w, u, info)
  if (info.ne.0) then 
     print *, 'error in do_sym_eig', info
     stop
  end if
  
  
  do kk=1, n
     if (w(kk)<=0) then 
        print *, 'error with eigenvalue', kk, w(kk)
        stop
     else
        w(kk)=1.0/sqrt(w(kk))
     end if
  end do
  
  print *, 'eig done' 

  do ii=1, n
     tmp1(1:n)=w(1:n)*u(jj,1:n)
     do jj=ii,n
        tmp2(1:n)=tmp1(1:n)*u(ii,1:n)
        sqrt_inv_a(ii, jj)=sum(tmp2)
        sqrt_inv_a(jj, ii)=sqrt_inv_a(ii, jj)
     end do
  end do
  
  print *, 'mutiply_done_inv' 
  
  deallocate(w)
  deallocate(u)
  deallocate(tmp1)
  deallocate(tmp2)
  
  
end subroutine do_sym_inv_sqrt


subroutine do_sym_sqrt_inv(a,n, sqrt_inv_a, nout, info)
  integer, intent(in)::n, nout
  real*8, intent(in)::a(n,n)
  real*8, intent(out)::sqrt_inv_a(nout,nout)
  integer, intent(out)::info
  
  ! local variable
  real*8, allocatable::w(:), u(:,:), v(:,:), ah(:,:)
  integer:: ii, jj, kk
  
  allocate(w(2*n), u(n,n), ah(n,n))
  ah=0.0
  ! call do_svd(a, m,n, u, w, v, info)

  call do_sym_eig(a, n,w, u, info)
  if (info.ne.0) then 
     print *, 'error in do_sym_eig', info
     stop
  end if
  
  do kk=1, n
     if (w(kk)<=0) then 
        print *, 'error with eigenvalue', kk, w(kk)
        stop
     else
        w(kk)=1.0/sqrt(w(kk))
     end if
  end do
  do ii=1, n
     do jj=1,n
        do kk=1, n
           ah(ii, jj)=ah(ii, jj)+u(ii,kk)*w(kk)*u(jj, kk)
        end do
     end do
  end do
  sqrt_inv_a(1:n, 1:n)=ah(1:n,1:n)

  deallocate(w)
  deallocate(u)
  deallocate(ah)
end subroutine do_sym_sqrt_inv

subroutine refill_bad_points_1d(prof, out_prof, nx, nz)
  integer, intent(in)::nx, nz
  
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(out)::out_prof(nx, nz)
  integer::i, k
  do i=1, nx
     do k=1, nz
        out_prof(i,k)=prof(i,k)
        if (prof(i,k)<-990.) out_prof(i,k)=out_prof(i,k-1)
        
     end do
  end do

end subroutine refill_bad_points_1d

subroutine get_spatial_err_cov(cor_len, error, lons, & 
     lats, err_cov, err_cov_inv, n)
! calculate spatial correlation between points defined by lons and lats. 

!f2py threadsafe
!f2py intent(in):: cor_len
!f2py intent(in):: n
!f2py intent(in):: error
!f2py intent(in):: lons
!f2py intent(in):: lats
!f2py intent(out):: err_cov
!f2py intent(out):: err_cov_inv

!f2py depend(error) n


implicit none
real*8::cor_len ! correlation length in km 
real*8::error(n) ! error of individule points 

real*8::lons(n), lats(n) ! location of the points 
real*8::err_cov(n,n), err_cov_inv(n,n)

integer::n

real*8, parameter::   g0 = 9.80665,&
     pi                       =  3.141592653589793115997963,&
     earth_major_axis         = 6378137.0, & 
     earth_minor_axis         = 6356752.3141, &
     earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
     (earth_minor_axis*earth_minor_axis)), &
     crad=pi/180.0
! local variables
real*4::rlon(8000), rlat(8000)

real*8::ran, rearth, sum_x, cri, cor_factor

integer:: ix, iy, kk


!  call do_sym_eig(a, n,w, u, info)


print *, n
! used to make a cut-off for point-point distance larger than 0.5*rearth
cri=cos(0.5)


! convert to radius

rlon(1:n)=crad*lons(1:n)
rlat(1:n)=crad*lats(1:n)

ran=0.0

rearth=0.5*(earth_major_axis+earth_minor_axis)/1000.0
rearth=rearth/cor_len

do ix=1, n
   err_cov(ix, ix)=error(ix)
   err_cov_inv(ix,ix)=1.0/sqrt(error(ix))
   
   do iy=ix+1, n
      err_cov(ix, iy)=0.0
      err_cov(iy, ix)=0.0
      ran=sin(rlat(ix))*sin(rlat(iy))+cos(rlat(ix))*cos(rlat(iy))*cos(rlon(ix)-rlon(iy))
      if (abs(ran)>1) ran=ran/ran
      

      if (abs(ran)>cri) then
         ran=acos(ran)
         cor_factor=exp(-rearth*ran)
         if (cor_factor>0.8) cor_factor=0.8
         if (cor_factor<5.0e-3) cor_factor=0.0
         
         err_cov(ix, iy)=sqrt(error(ix)*error(iy))*cor_factor
         err_cov(iy, ix)= err_cov(ix, iy)
         sum_x=sum_x+err_cov(ix,iy)*err_cov(ix,iy)
      end if
      
   end do
   sum_x=sqrt(error(ix)/sum_x)
   
   
  
   
end do
print *, 'err_cov', err_cov(1,1), err_cov(2,2)

print *, 'Done construnt'

! call do_sys_sqrt_inv(err_cov,n, err_cov_inv)
! call do_svd_sqrt(err_cov,n, err_cov_inv)
! call  do_ltrm_inv(err_cov, n, err_cov_inv)

print *, 'Done inv'
print *, 'err_cov', err_cov(1,1), err_cov(2,2)


end subroutine get_spatial_err_cov


subroutine get_spatial_err_cov_xy(cor_len, err_x, lon_x, lat_x, & 
     err_y, lon_y, lat_y, err_cor, nx, ny)
! calculate spatial correlation between two set of points defined by lon_x and lat_x, lon_y, and lat_y
 

real*8,  intent(in):: cor_len
integer, intent(in):: nx, ny
real*8,  intent(in):: err_x(nx), lon_x(nx), lat_x(nx)
real*8,  intent(in):: err_y(ny), lon_y(ny), lat_y(ny)
real*8,  intent(out):: err_cor




real*8, parameter::   g0 = 9.80665,&
     pi                       =  3.141592653589793115997963,&
     earth_major_axis         = 6378137.0, & 
     earth_minor_axis         = 6356752.3141, &
     earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
     (earth_minor_axis*earth_minor_axis)), &
     crad=pi/180.0
! local variables
real*4::rlon1(4000), rlat1(4000)
real*4::rlon2(4000), rlat2(4000)

real*8::ran, rearth, cri, cor_factor,err1, err2

integer:: ix, iy


!  call do_sym_eig(a, n,w, u, info)



cri=cos(0.5)


! convert to radius

rlon1(1:nx)=crad*lon_x(1:nx)
rlat1(1:nx)=crad*lat_x(1:nx)

rlon2(1:ny)=crad*lon_y(1:ny)
rlat2(1:ny)=crad*lat_y(1:ny)


ran=0.0

rearth=0.5*(earth_major_axis+earth_minor_axis)/1000.0
rearth=rearth/cor_len
err1=sum(err_x)
err2=sum(err_y)

do ix=1, nx
   do iy=1, ny
      ran=sin(rlat1(ix))*sin(rlat2(iy))+cos(rlat1(ix))*cos(rlat2(iy))*cos(rlon1(ix)-rlon2(iy))
      
      if (abs(ran)>1) ran=ran/ran
      if (abs(ran)>cri) then
         ran=acos(ran)
         cor_factor=exp(-rearth*ran)
         
         if (cor_factor>0.99) cor_factor=0.99
         if (cor_factor<5.0e-2) cor_factor=0.0
         
         err_cor=err_cor+err_x(ix)*err_y(iy)*cor_factor
         
      end if
   end do
end do
err_cor=err_cor/(err1*err2)


end subroutine get_spatial_err_cov_xy



subroutine get_spatial_err_cov_new(cor_len, total_error, error, lons, & 
     lats, err_cov, err_cov_inv, n)

! calculate spatial correlation between points defined by lons and lats. 
! this version separates the correlated and uncorrelated parts 
  

!f2py intent(in):: cor_len
!f2py intent(in):: n
!f2py intent(in):: error
!f2py intent(in):: total_error
!f2py intent(in):: lons
!f2py intent(in):: lats
!f2py intent(out):: err_cov
!f2py intent(out):: err_cov_inv
  


implicit none
real*8::cor_len ! correlation length in km 
real*8::error(n), total_error(n) ! error of individule points 

real*8::lons(n), lats(n) ! location of the points 
real*8::err_cov(n,n), err_cov_inv(n,n)

integer::n

real*8, parameter::   g0 = 9.80665,&
     pi                       =  3.141592653589793115997963,&
     earth_major_axis         = 6378137.0, & 
     earth_minor_axis         = 6356752.3141, &
     earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
     (earth_minor_axis*earth_minor_axis)), &
     crad=pi/180.0
! local variables
real*4::rlon(8000), rlat(8000)

real*8::ran, rearth, cri, cor_factor

integer:: ix, iy, kk


!  call do_sym_eig(a, n,w, u, info)


print *, n

cri=cos(0.5)


! convert to radius

rlon(1:n)=crad*lons(1:n)
rlat(1:n)=crad*lats(1:n)

ran=0.0

rearth=0.5*(earth_major_axis+earth_minor_axis)/1000.0
rearth=rearth/cor_len

do ix=1, n
   err_cov(ix, ix)=total_error(ix)
   err_cov_inv(ix,ix)=1.0/sqrt(total_error(ix))
   
   do iy=ix+1, n
      err_cov(ix, iy)=0.0
      err_cov(iy, ix)=0.0
      ran=sin(rlat(ix))*sin(rlat(iy))+cos(rlat(ix))*cos(rlat(iy))*cos(rlon(ix)-rlon(iy))
      
      if (abs(ran)>1) ran=ran/ran
      
      
      if (abs(ran)>cri) then
         ran=acos(ran)
         cor_factor=exp(-rearth*ran)
         if (cor_factor>0.8) cor_factor=0.8
         if (cor_factor<5.0e-3) cor_factor=0.0
         
         err_cov(ix, iy)=sqrt(error(ix)*error(iy))*cor_factor
         err_cov(iy, ix)= err_cov(ix, iy)
         ! sum_x=sum_x+err_cov(ix,iy)*err_cov(ix,iy)
      end if
      
   end do
   ! sum_x=sqrt(error(ix)/sum_x)
   
end do


print *, 'err_cov', err_cov(1,1), err_cov(2,2)

print *, 'Done construnt'

! call do_sys_sqrt_inv(err_cov,n, err_cov_inv)
! call do_svd_sqrt(err_cov,n, err_cov_inv)
! call  do_ltrm_inv(err_cov, n, err_cov_inv)

print *, 'Done inv'
print *, 'err_cov', err_cov(1,1), err_cov(2,2)


end subroutine get_spatial_err_cov_new


subroutine check_time_cov(time_wd, otime, error, error_md, nx) 
!remove correlation outside the time window implicit none integer,
implicit none
integer, intent(in)::nx 
real*8, intent(in)::time_wd
real*8, intent(in)::otime(nx), error(nx, nx)
real*8, intent(out)::error_md(nx,nx)
integer::i, j
real*8::dt 
error_md=error 
do i=1, nx
   do j=1,nx
      dt=abs(otime(j)-otime(i))
      if (dt>time_wd) error_md(i,j)=0.0
   end do
end do

end subroutine check_time_cov



subroutine get_idx(x0, dx, x, idx)
implicit none
real*8, intent(in)::x0
real*8, intent(in)::dx
real*8, intent(in)::x
integer, intent(out)::idx
idx=int((x-x0)/dx)+1
end subroutine get_idx

subroutine form_rh_rdy(hm, bm, xm, lvls, nobs, nlvl, nx, obs, & 
     oerr, oap_r, opres, pres_st, pres_end, rh, rdy, np)
! 1) scaling h to hm_scaled=R^(-1/2)*h*b^(-1/2)
! 2) reform the matrix  from hm_scaled(nobs, nlvl, nx) to  rh(nobs*nlvl, nx)
! 3) scaling dy to dy_scaled=R^(-1/2)*dy
! 4) reform the matrix  from dy_scaled(nobs, nlvl) to  rdy(nobs*nlvl)
! hm ------in------array of (nobs, nlvl, nx) Jacobian matrix 
! bm ------in----- array of (nx, nx)   square root of B 
! xm  ---- in ---- array of (nx) prior x  values
! obs ---- in ---- array of (nobs, nlvl) observation values
! oerr ----in----- array of (nobs, nlvl, nlvl) observation error covaraince
! oap_r  ---- in ---- array of (nobs, nlvl) 'transformed' apriori vector oap_r=(1-A)*oap
! opres_r -----in ------ array of (nobs, nlvl)  observation pressure 
! pres_st, pres_end -----in --------- selected pressure levels  
! rh        -------out -------- array of (np, nx) reformed jacobian 
! rdy    -----out ----------- array of (np)   reformed dy 
! 

  
implicit none
integer, parameter::max_lvl=7
integer, parameter::max_obs=5000
integer, intent(in)::nobs, nlvl, nx
real*8, intent(in)::hm(nobs, nlvl, nx), bm(nx,nx), xm(nx), opres(nobs, nlvl)
real*8, intent(in)::obs(nobs, nlvl), oerr(nobs, nlvl, nlvl), oap_r(nobs, nlvl)
integer, intent(in)::lvls(nobs)
real*8, intent(out)::rh(max_obs*max_lvl, nx)
real*8, intent(out)::rdy(max_obs*max_lvl)
real*8, intent(in)::pres_st, pres_end
integer, intent(out)::np

! local variables
integer::iobs, ii, jj, kk, ilvl
integer::usd_lvl, ip
integer::lvl_st, lvl_end
logical::use_it
integer::info

real*8 ::rr(max_lvl, max_lvl), arr(max_lvl, max_lvl)
real*8, allocatable::new_hm(:,:), new_y(:)

! usd_lvl=lvl_end-lvl_st+1


allocate(new_hm(max_lvl, nx), new_y(max_lvl))


ip=0
do iobs=1, nobs
   ! adjustment according to new_values
   lvl_st=nlvl
   lvl_end=1
   usd_lvl=0
   do ilvl=1, nlvl
      if (opres(iobs,ilvl)>0) then
         if (abs(opres(iobs,ilvl)-pres_st)<0.1*opres(iobs, ilvl)) lvl_st=ilvl
         usd_lvl=usd_lvl+1
      end if
      
      if (opres(iobs,nlvl-ilvl+1)>0) then
         ! print *, iobs, abs(opres(iobs, nlvl-ilvl+1)-pres_end), pres_end
         if (abs(opres(iobs, nlvl-ilvl+1)-pres_end)<0.1*opres(iobs, nlvl-ilvl+1)) then 
            ! print *, iobs, nlvl-ilvl+1, abs(opres(iobs, nlvl-ilvl+1)-pres_end), opres(iobs, nlvl-ilvl+1)
            lvl_end=nlvl-ilvl+1
         end if
      end if
      
   end do
   
!    print *, 'iobs, lvl_st, lvl_end, usd_lvl', iobs,  lvl_st, lvl_end, usd_lvl
   use_it=.FALSE.
   if ((usd_lvl>=lvl_end-lvl_st+1).and.(lvl_end>=lvl_st)) then 
      
      rr(1:usd_lvl, 1:usd_lvl)=oerr(iobs, 1:usd_lvl, 1:usd_lvl)
      ! 
      call do_sym_sqrt_inv(rr(1:usd_lvl, 1:usd_lvl),usd_lvl,& 
           arr, max_lvl, info)
      if (info==0) then 
         use_it=.TRUE.
      else
         stop
      end if
      
   end if
   
   if (use_it) then 
      ! new_hm=hm*bm
      
      do kk=1, usd_lvl
         do jj=1, nx
            new_hm(kk,jj)=0.0         
            do ii=1, nx
               new_hm(kk,jj)=new_hm(kk,jj)+hm(iobs, kk, ii)*bm(ii, jj)
            end do
         end do
         new_y(kk)=0.0
      !  new_y   =  h* x
         
         do jj=1, nx
            new_y(kk)=new_y(kk)+hm(iobs, kk, jj)*xm(jj)
         end do

      end do
   


      ! usd_lvl=lvl_end-lvl_st+1
      
      do ilvl=lvl_st, lvl_end
         
         ip=ip+1
         rdy(ip)=0.0
         ! dy=a*(obs-oap_r-new_y)
         ! rdy=r^(-1/2)*dy
         
         do kk=1, usd_lvl
            rdy(ip)=rdy(ip)+arr(ilvl, kk)*(obs(iobs, kk)-oap_r(iobs, kk)-new_y(kk))
         end do
         
         ! dy=a*(obs-oap_r-new_y)
         ! rdy=r^(-1/2)*dy
      
         do jj=1, nx
            rh(ip, jj)=0.0
            do kk=1, usd_lvl
               rh(ip, jj)=rh(ip, jj)+arr(ilvl, kk)*new_hm(kk, jj)
            end do
         end do
      end do
   end if
end do
np=ip
deallocate(new_hm)
deallocate(new_y)

end subroutine form_rh_rdy

 
subroutine get_dx_tm(bm, rh, rdy, nx, nobs, dx, tm)
! calculate the analysis increment dx and transform matrix tm
! bm --------in------- array  of (nx, nx)  bm=sqrt(B) the square root of priori error covariance
! rh -------in ------  array of rh(nobs, nx) pre-weighted jacobian rh=R^(-1/2)*H*B^(1/2)
! rdy ------ in ------ array of (nobs)  scaled dy  dy=R^(-1/2)*dy
! dx  ------ out ------- array of (nx) analysis increment 
!                dx=bm*transpose(rh)[rh*transpose(rh)+1]^(-1)*rdy
! tm ------- out ------- array of (nx, nx) transform matrix to calculate 
! the sqare root of posteriori error covariance via bm(a)=bm(f)*tm
!                tm=sqrt(1-transpose(rh)[rh*transpose(rh)+1]^(-1)*rh)
   
implicit none
integer, intent(in)::nx, nobs
real*8, intent(in)::bm(nx,nx), rh(nobs, nx), rdy(nobs)
real*8, intent(out)::dx(nx), tm(nx,nx)
! local variables
real*8, allocatable::uy(:,:), vy(:,:), wy(:), wd(:), rht(:,:),  vrdy(:),uwy(:)
integer::ny
integer::ii, jj, kk,info, nusd

ny=nobs
allocate(rht(nx, ny))

rht=transpose(rh(1:ny, 1:nx))
allocate(uy(nx,nx), wy(nx+ny), vy(ny, ny))

call do_svd(rht, nx, ny, uy, wy, vy, info)
if (info.ne.0) then
   print *, 'error in svd', info
   stop
end if

! choose the usd number of singular vaules (the columns in U, and rows in V)
 
nusd=0

do ii=1, min(nx,ny)
   if (wy(ii)>1.0e-30) nusd=ii
end do

allocate(wd(ny))
wd=1.0

! wd=1.0/(1+transpose(wy)*wy) wd is in  size ny

do ii=1, nusd
   wd(ii)=wd(ii)+wy(ii)*wy(ii)
   wd(ii)=1.0/wd(ii)
end do


allocate(vrdy(nusd)) ! VR*DY 


vrdy=0.0
! only 1:nusd is useful when  mutippled by wd & then wy

!$OMP PARALLEL DO 

do ii=1, nusd
   do jj=1, ny
      vrdy(ii)=vrdy(ii)+vy(ii, jj)*rdy(jj)
   end do
end do
!$OMP END PARALLEL DO

! effectively reduced vrdy from ny to nx with only 1:nusd not equal to zeros
! by multiplied with wy 
vrdy(1:nusd)=wy(1:nusd)*wd(1:nusd)*vrdy(1:nusd) 

allocate(uwy(nx))

uwy=0.0

! u*wy*wd*v*rdy 
!$OMP PARALLEL DO 

do ii=1, nx
   do jj=1, nusd
      uwy(ii)=uwy(ii)+uy(ii, jj)*vrdy(jj)
   end do
end do
!$OMP END PARALLEL DO

dx=0.0

do ii=1, nx
   do jj=1, nx
      dx(ii)=dx(ii)+bm(ii, jj)*uwy(jj)
   end do
end do

! deallocate(uyw)

! calculate transform matrix tm 
tm=0.0

! allocate(uyw(nx))

! for the terms of nusd: 1.0-wy*[1+transpose(wy)*wy]^(-1)*transpose(wy)
! effectively only 1:nusd terms are different from 1

uwy=1.0 

do jj=1, nusd
   uwy(jj)=sqrt(1.0-wy(jj)*wd(jj)*wy(jj))
end do

! tm=uy*uwy*transpose(uy)
!$OMP PARALLEL DO 

do ii=1, nx
   do jj=1, nx
      tm(ii, jj)=0.0
      do kk=1, nx
         tm(ii,jj)=tm(ii,jj)+uy(ii, kk)*uwy(kk)*uy(jj, kk)
      end do
      
   end do
end do
!$OMP END PARALLEL DO

deallocate(uy)
deallocate(wy)
deallocate(vy)

deallocate(wd)
deallocate(vrdy)
deallocate(uwy)
deallocate(rht)


end subroutine get_dx_tm
  

subroutine reform_h(h, lvls, idx, nobs, nlvl, nx0, hm, nx)
! reduce jacobian h with respect to nx0 varaibles to  
! hm with respect to nx (nx<nx0) variables
! h ------in------- array of (nobs, nx0, nlvl) for jacobians 
! idx -----in ------ array of nx0: the position of the variables in the new order
! lvls -----in------ array of (nobs)  number of used levels
! hm -----out ------array of (nobs, nlvl, nx0)  the reformed jacobian matrix 

implicit none
integer, intent(in)::nobs, nlvl, nx0, nx
real*8, intent(in)::h(nobs, nx0, nlvl)
integer, intent(in)::idx(nx0)
integer, intent(in)::lvls(nobs)
real*8, intent(out)::hm(nobs, nlvl,nx)
integer::iobs, ix, jj, kk, ml
hm=0.0
!$OMP PARALLEL DO 

do kk=1, nx0
   ix=idx(kk)  
   if (ix>0) then 
      do iobs=1, nobs
         ml=lvls(nobs)
         do jj=1, ml
            hm(iobs, jj, ix)=hm(iobs, jj, ix)+h(iobs, kk, jj)
         end do
      end do
   end if
end do
!$OMP END PARALLEL DO
end subroutine reform_h

subroutine get_vertical_wgt_half(z0, z, nx, nz0, nz, pl, pr, wgt)
  ! position and wgt for vertical interpolation over nx points
    
  ! z0 -----in-----  old vertical coordinate 
  ! z ------in------ new vertical coordiante
  ! nx ------in------- size in x direction 
  ! nz0 ----in--------- size in z direction 
  ! nz -----in --------- size in the out direction
  implicit none
  integer, intent(in)::nx, nz0, nz
  real*8, intent(in)::z0(nx, nz0)
  real*8, intent(in)::z(nz)
  integer*4, intent(out)::pl(nx, nz)
  integer*4, intent(out)::pr(nx, nz)
  real*8, intent(out)::wgt(nx,nz)
  ! local variables
  integer::i
  integer*4, allocatable::pl1d(:), pr1d(:)
  real*8, allocatable::wgt1d(:)
  
  
  allocate(pl1d(nz), pr1d(nz), wgt1d(nz))
  
  
  do i=1, nx
     call getwgt(z0(i, :), z(:), nz0, nz, pl1d, pr1d,wgt1d)
     pl(i,:)=pl1d(:)
     pr(i,:)=pr1d(:)
     wgt(i,:)=wgt1d(:)
  end do
  
  deallocate(pl1d, pr1d, wgt1d)
  
  
end subroutine get_vertical_wgt_half

subroutine col_int_half(prof, colwgt, nx, nz, colval)
  
  ! calculate the column values for list of profiles using the weighting factor
  ! prof   -------in-------- array of (nx, nz)==volume mixing ratios
  ! colwagt ------in ------- array of (nx, nz)== weighting factor for column integration 
  ! colval -------out ------- array of (nx) == column mixing values at nx locations. 
    
  implicit none
  integer, intent(in)::nx, nz
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::colwgt(nz)
  real*8, intent(out)::colval(nx)
  ! local variable
  real*8 :: sum_wgt
  integer::i, k
  
  do i=1, nx
     colval(i)=0
     sum_wgt=0.0
     
     do k=1, nz
        if (prof(i, k)>-990.0) then 
           colval(i)=colval(i)+colwgt(k)*prof(i, k)
           sum_wgt=sum_wgt+colwgt(k)
        end if
     end do

     if (sum_wgt>0) then
        colval(i)= colval(i)/sum_wgt
     else
        colval(i)=-999.0
     end if
  end do
end subroutine col_int_half

subroutine get_co_locate_points(lon, lat, n, pcor,puse)
! find points with same co-ordinates 
!
! lat lon ------in------array of size n 
!   ------in--------- array of size nx
! pcor ------out-------- array of size n.  which point is co-located 
! puse ------out-------- array of size n.  -1 equal to with co-located
 


implicit none 
integer, intent(in)::n
real*8, intent(in)::lon(n), lat(n)
integer*4, intent(out)::pcor(n), puse(n)   
integer::i, j

j=1
do i=1, n
   pcor(i)=0
   puse(i)=1
   do j=i+1, n
      if ((lon(i)==lon(j)).and.(lat(i)==lat(j))) then 
         pcor(i)=j
         puse(i)=-1
      end if
   end do
end do

end subroutine get_co_locate_points


subroutine combine_vertical_grid(mod_pres, obs_pres, obs_ak, psurf, pres_cb, &
     obs_pres_out,  obs_ak_out, nmod, olvl, np)
  
  implicit none
  
  integer, intent(in)::nmod, olvl, np
  real*8, intent(in):: mod_pres(np,nmod), obs_pres(np,olvl), & 
       obs_ak(np,olvl), psurf(np)
  
  real*8, intent(out)::pres_cb(np,nmod + olvl), & 
       obs_pres_out(np, olvl), obs_ak_out(np, olvl)
  
  integer::ilvl, next_mod, lm, cur_mod, cur_ol, ip
  real*8:: cur_obs_pres, next_obs_pres
  
  obs_pres_out=obs_pres
  obs_ak_out=obs_ak
  
  pres_cb=-999.0
  
  do ip=1, np
     ilvl=1
     cur_ol=1
     next_mod=0
     pres_cb(ip,ilvl)=obs_pres(ip,cur_ol)
     cur_obs_pres=obs_pres(ip,cur_ol)
     
     do cur_ol=2, olvl
        ! data will not be used
        next_obs_pres=obs_pres(ip,cur_ol)
        if (next_obs_pres>=psurf(ip)) exit
        
        if ((obs_ak(ip, cur_ol)<-990.0).or.(next_obs_pres<-990.0)) exit
        
        cur_mod=next_mod+1
        do lm=cur_mod, nmod
           if ((mod_pres(ip, lm)>cur_obs_pres).and.(mod_pres(ip,lm)<next_obs_pres)) then 
              ilvl=ilvl+1
              pres_cb(ip,ilvl)=mod_pres(ip,lm)
              next_mod=lm
           else if (mod_pres(ip,lm)>=next_obs_pres) then 
              exit
           end if
        end do
        
        ilvl=ilvl+1
        pres_cb(ip,ilvl)=next_obs_pres
        cur_obs_pres=next_obs_pres
     end do
     
     ! replace the bottom one with psurf
     if (cur_ol<olvl) then 
        obs_pres_out(ip, cur_ol)=psurf(ip)
        if (obs_ak(ip, cur_ol)<0.01) then
           obs_ak_out(ip, cur_ol)=obs_ak(ip, cur_ol-1)
        else
           obs_ak_out(ip, cur_ol)=obs_ak(ip, cur_ol)
        end if
        obs_pres_out(ip, cur_ol+1:)=-999.0
        obs_ak_out(ip, cur_ol+1:)=-999.0
        
        pres_cb(ip,ilvl)=psurf(ip)
        
     else 
        obs_pres_out(ip, olvl)=psurf(ip)
        obs_ak_out(ip, olvl)=obs_ak(ip, olvl-1)
        pres_cb(ip,ilvl)=psurf(ip)
     
     end if
     
     
     
  end do
  
end subroutine combine_vertical_grid

subroutine array_divide_array(x, y, add_val, yout, nx)
implicit none
integer, intent(in)::nx
real*8, intent(in)::x(nx), y(nx)
real*8, intent(in)::add_val
real*8, intent(out)::yout(nx)
integer::i
do i=1, nx
   if ((x(i)>-990.0).and.(y(i)>-990.0)) then
      yout(i)=x(i)/(add_val+y(i))
   else
      yout(i)=-999.0
   end if
end do
end  subroutine array_divide_array

subroutine array_divide_array_2d(x, y, add_val, yout, ny, nx)
implicit none
integer, intent(in)::ny, nx
real*8, intent(in)::x(ny, nx), y(ny, nx)
real*8, intent(in)::add_val
real*8, intent(out)::yout(ny,nx)
integer::i,j
do j=1,ny
   do i=1, nx
      if ((x(j, i)>-990.0).and.(y(j,i)>-990.0)) then
         yout(j,i)=x(j,i)/(add_val+y(j,i))
      else
         yout(j,i)=-999.0
      end if
   end do
end do
end subroutine array_divide_array_2d

subroutine array_divide_array_2d_1d(x, y, add_val, yout, ny, nx)
implicit none
integer, intent(in)::ny, nx
real*8, intent(in)::x(ny, nx), y(nx)
real*8, intent(in)::add_val
real*8, intent(out)::yout(ny,nx)
integer::i,j
do j=1,ny
   do i=1, nx
      if ((x(j, i)>-990.0).and.(y(i)>-990.0)) then
         yout(j,i)=x(j,i)/(add_val+y(i))
      else
         yout(j,i)=-999.0
      end if
   end do
end do
end subroutine array_divide_array_2d_1d
