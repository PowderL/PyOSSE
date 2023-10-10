! functions for calculate column values from profiles of tracer mixing ratios. 
!
! Fucntions: 
!
!  1. get_col_wgt_edge: 
!  Airmass weighting factor for calculating total columns from mixing ratio profiles,
!  when pressure is given at level edges.  
!
!  2. get_col_wgt_edge_mask: same as  get_col_wgt_edge, but for masked pressure array. 
!
!  3. get_col_wgt:  
!  Airmass weighting factor for calculating total columns from mixing ratio profiles,
!  when pressure is given at level centers. 
! 
!  4. get_col_wgt_mask: same as  get_col_wgt but for masked pressure array
! 
!  5. get_col_wgt_0d:  version for calculating column weight for single pressure profiles  
!
!  6. get_col_wgt_1d:  calculating column weights for pressure profiles along a track. 
!  
!  7. get_col_wgt_2d:  calculating column weights for pressure profiles over a 2D region.
!
!  8. col_int_0d:      calculating column for single profile 
!
!  9. col_int_1d:      calculating columns for profiles along a track
!
!  10. col_int_2d:      calculating columns for profiles over a 2D regions. 
!
!  11. col_int_0d_em:   As col_int_0d, but for multiple tracers
!
!  12. col_int_1d_em:   As col_int_1d, but for multiple tracers
!
!  13. col_int_2d_em:   As col_int_2d, but for multiple tracers
!
!  14. ak_col_int_0d:   calculating column for single profile weigted by instrument averaging kernels   
!
!  15. ak_col_int_1d:    calculating columns weighted by instrument kernels for multiple profiles!  16. ak_col_int_2d:   calculating columns weighted by instrument kernels for multiple profiles
!
!  17. ak_col_int_0d_em: ensemble version of ak_col_int_0d
!
!  18  ak_col_int_1d_em: ensemble version of ak_col_int_1d
!
!  19. ak_col_int_2d_em: ensemble version of ak_col_int_2d

!  




 
 




subroutine get_col_wgt_edge(nz, pres_in, wgt_out)
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  ! 
  ! Inputs: 
  ! -----------------------------------
  ! 1. integer*4 ::nz #  number of levels 
  ! 2. real*8    :: pres_in(nz) # pressure in hPa or log10(pressure/hPa). It is assumed that 
  !  pres_in is given in ascending order. 
  ! 
  !
  ! Outputs
  ! -------------------------------------
  ! 1. real*8 :: wgt_out(nz-1) # mass weighting factor sum(wgt_out)=1.0
  ! 
  ! 
  ! Notes:
  ! -------------------------------------------------------
  ! 1. Profile are assumed to be given at level centre, and pressure at level edge, so the pressure level is larger than profile level by 1
  ! 
  ! 2. weighting factor is given dp/P(top)-P(bottom)
  ! 
  !  
  
  implicit none
  
  ! Input
  integer, intent(in)::nz
  real*8, intent(in) ::pres_in(nz)
 
  ! Output
  real*8, intent(out)::wgt_out(nz-1)
  
  
  ! local
  real*8, allocatable::pres(:), logp(:)
  
  
  integer::i
  
  wgt_out(:)=0.0
  allocate(pres(nz), logp(nz))
  
  
  if (maxval(pres_in)>10) then 
     pres=pres_in
     logp=log10(pres)
  else
     pres=10.0**(pres)
     logp=pres_in
  end if
  
  do i=1, nz-1
     wgt_out(i)=pres(i+1)-pres(i)
  end do
  
  wgt_out=wgt_out/(pres(nz)-pres(1))
  
  deallocate(pres, logp)
  
end subroutine get_col_wgt_edge


subroutine get_col_wgt_edge_mask(nz, pres_in, wgt_out, mask_val)
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  ! 
  ! Inputs: 
  ! -----------------------------------
  ! 1. integer*4 ::nz #  number of levels 
  ! 2. real*8    :: pres_in(nz) # pressure in hPa or log10(pressure/hPa). It is assumed that 
  !  pres_in is given in ascending order. 
  ! 
  !
  ! Outputs
  ! -------------------------------------
  ! 1. real*8 :: wgt_out(nz-1) # mass weighting factor sum(wgt_out)=1.0
  ! 
  ! 
  ! Notes:
  ! -------------------------------------------------------
  ! 1. Profile are assumed to be given at level centre, and pressure at level edge, so the pressure level is larger than profile level by 1
  ! 
  ! 2. weighting factor is given dp/P(top)-P(bottom)
  ! 
  !  
  
  implicit none
  
  ! Input
  
  integer, intent(in)::nz
  real*8, intent(in)::pres_in(nz)
  real*8, intent(in), optional::mask_val
  
  ! Output
  real*8, intent(out)::wgt_out(nz-1)
  
  
  ! local
  real*8, allocatable::pres(:), logp(:), pres_a(:)
  integer::allocatabpidx
  real*8::fmask
  
  integer::i, ist, iend
  
  allocate(pres(nz), logp(nz), pres_a(nz))

 
  wgt_out(:)=0.0
  
  
  if (maxval(pres_in)>10) then 
     pres=pres_in
     logp=log10(pres)
  else
     pres=10.0**(pres)
     logp=pres_in
  end if
  
  do i=1, nz-1
     wgt_out(i)=pres(i+1)-pres(i)
  end do
  
  wgt_out=wgt_out/(pres(nz)-pres(1))
  
  deallocate(pres, logp)
  
end subroutine get_col_wgt_edge_mask



subroutine get_col_wgt(nz, pres_in, wgt_out)
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  !
  ! Inputs: 
  ! -----------------------------------
  ! 1. integer*4 ::nz #  number of levels 
  ! 2. real*8    :: pres_in(nz) # pressure in hPa or log10(pressure/hPa). It is assumed that 
  !  pres_in is given in ascending order. 
  ! 
  !
  ! Outputs
  ! -------------------------------------
  ! 1. real*8 :: wgt_out(nz) # mass weighting factor sum(wgt_out(nz))=1.0
  ! 
  ! Notes:
  ! -------------------------------------------------------
  ! 1. Profile and pressure are assumed to be given at grid edge, and have the same number of  vertical levels 
  
  ! 2. Vertical integration sum(fdz)=-sum(fdp/g)=sum_over_i{f(i+1)[p(i+1)-(p(i+1)-p(i))/ln(p(i+1)/p(i1))]-f(i)[p(i)-(p(i+1)-p(i))/ln(p(i+1)/p(i1))]}/g
  ! 3. the first term f(i+1)p(i+1) of each level will be canceled out when summed over different layers,  except the bottom and top one 
  ! 
  ! 
  
  implicit none
  ! Input
  integer, intent(in)::nz
  real*8, intent(in) ::pres_in(nz)
 
  ! Output
  real*8, intent(out)::wgt_out(nz)
  
  ! local

  real*8, allocatable::pres_a(:), pres(:), fsum(:), logp(:)
  real*8, allocatable::wgt(:)
  
  real*8  :: a=log(10.0)
  real*8::w1, w2
  integer::i, nz_real, ist, iend
  


  wgt_out(:)=0.0
  
  nz_real=0
  ist=1
  iend=nz
  
  ! remove repeated elements at the bottom and top of the pressure array
  
  do i=1, nz
     if (pres_in(i)==pres_in(ist)) then 
        ! last one == pres_in(1)
        ist=i
     else if (pres_in(i)==pres_in(iend)) then
        ! first one ==pres_in(nz)
        iend=i
        exit
     end if
  end do
  
  nz_real=iend-ist+1
  
  
  allocate(pres_a(nz_real), logp(nz_real), pres(nz_real), wgt(nz_real))
  allocate(fsum(nz_real))
  wgt=0.0
  
  pres(1:nz_real)=pres_in(ist:iend)
! guess it is log10(pressure) or pressure in hPa

  if (maxval(pres)>10) then 
     pres_a=pres
     logp=log10(pres)
  else
     pres_a=10.0**(pres)
     logp=pres
  end if

  fsum=0.0
  ! log_10 to ln  
  a=log(10.0) 
  ! weighting factor for single level 
  
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
  
  ! correction for first level
  wgt(1)=-pres_a(1)+(pres_a(2)-pres_a(1))/(a*(logp(2)-logp(1)))
  
  ! correction for top level 
  i=nz_real
  wgt(i)=pres_a(i)-(pres_a(i-1)-pres_a(i))/(a*(logp(i-1)-logp(i)))
 
  !  Normalized with P(top)-P(bottom)
  
  wgt(1:nz_real)=wgt(1:nz_real)/(pres_a(nz_real)-pres_a(1))
  wgt_out(ist:iend)=wgt(1:nz_real)
  
  deallocate(fsum)
  deallocate(pres_a, logp, pres, wgt)

end subroutine get_col_wgt

subroutine get_col_wgt_mask(nz, pres, wgt, mask_val)
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles. 
  ! Only the levels with pressure<>mask_val will be used
  ! 
  ! Inputs: 
  ! 1. integer*4 ::nz #  number of levels 
  ! 
  ! 2. real*8    :: pres_in(nz) # pres in hPa or log10(pressure/hPa). It is assumed that 
  !  pres_in is given in ascending order. 
  ! 3. real*8,  mask_val, optional::value for bad or missing data
  
  !
  ! Outputs
  ! 1. real*8 :: wgt(nz)
  
  
  implicit none
  
  integer, intent(in)::nz
  real*8, intent(in) ::pres(nz)
  real*8, optional, intent(in) ::mask_val
  ! output 
  
  real*8, intent(out)::wgt(nz)
  
  ! local variables
  
  real*8, allocatable::pres_use(:) ! pressure for levels in use
  real*8, allocatable::wgt_use(:)  ! weighting for levels in use
  integer, allocatable::pres_idx(:) ! the index for levels in use
  
  integer::i, k, nz_use
  
  real*8::fmask
  
  if (present(mask_val)) then
     fmask=mask_val

  else
     ! set to the default value
     fmask=-999.0
  end if
  
  
   
  ! check number of levels in use
  nz_use=0
  allocate(pres_use(nz), pres_idx(nz))
  k=0
  do i=1, nz
     if (pres(i)/=fmask) then
        nz_use=nz_use+1
        k=k+1
        pres_use(k)=pres(i)
        pres_idx(k)=i
     end if
  end do
  
  allocate(wgt_use(nz_use))
  wgt_use=0.0
  call get_col_wgt(nz_use, pres_use(1:nz_use),  wgt_use)
  
  wgt=0.0
  ! assign to output values
  
  do i=1, nz_use
     k=pres_idx(i)
     wgt(k)=wgt_use(i)
  end do
  
  deallocate(pres_use, wgt_use, pres_idx)
  
end subroutine get_col_wgt_mask


subroutine get_col_wgt_0d(nz, pres, wgt, mask_val)
  !  wrapper for get_col_wgt_mask
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of points 
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: pres(nx, nz) # pres in hPa or log10(pressure/hPa). It is assumed that 
  !  presis given in ascending order. 
  
  !
  ! Outputs
  ! 4. real*8 :: wgt(nx, nz): mass weighting factor 
  


  implicit none
  ! in 
  integer, intent(in)::nz
  real*8, intent(in) ::pres(nz)
  ! out
  real*8, intent(out)::wgt(nz)
  real*8, intent(in), optional::mask_val
  
  ! local variables
  
  real*8::fmask
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  call get_col_wgt_mask(nz, pres, wgt, fmask)
  
end subroutine get_col_wgt_0d



subroutine get_col_wgt_1d(nx, nz, pres, wgt, mask_val)
  
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of horizontal points 
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: pres(nx, nz) # pres in hPa or log10(pressure/hPa). It is assumed that 
  !  presis given in ascending order. 
  
  !
  ! Outputs
  ! 4. real*8 :: wgt(nx, nz): mass weighting factor 
  


  implicit none
  ! in 
  integer, intent(in)::nz, nx
  real*8, intent(in) ::pres(nx, nz)
  ! out
  real*8, intent(out)::wgt(nx, nz)
  real*8, intent(in), optional::mask_val
  
  ! local variables
  real*8, allocatable::cur_pres(:), cur_wgt(:)
  integer::i
  real*8::fmask
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  allocate(cur_pres(nz), cur_wgt(nz))
  
  
  do i=1, nx
     cur_pres(:)=pres(i,:)
     call get_col_wgt_mask(nz, cur_pres, cur_wgt, fmask)
     wgt(i,:)=cur_wgt(:)
  end do
  
  deallocate(cur_pres, cur_wgt)

end subroutine get_col_wgt_1d


subroutine get_col_wgt_2d(nx, ny, nz, pres, wgt, mask_val)
  
  !  Airmass weighting factor for calculating total columns from mixing ratio profiles, 
  ! Inputs: 
  ! 1. integer*4 ::nx #  number of horizontal points 
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: pres(nx, nz) # pres in hPa or log10(pressure/hPa). It is assumed that 
  !  presis given in ascending order. 
  
  !
  ! Outputs
  ! 4. real*8 :: wgt(nx, nz): mass weighting factor 
  


  implicit none
  ! in 
  integer, intent(in)::nz, ny, nx
  real*8, intent(in) ::pres(nx, ny, nz)
  ! out
  real*8, intent(out)::wgt(nx, ny, nz)
  real*8, intent(in), optional::mask_val
  
  ! local variables
  real*8, allocatable::cur_pres(:), cur_wgt(:)
  integer::i, j
  real*8::fmask
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  allocate(cur_pres(nz), cur_wgt(nz))
  
  
  do i=1, nx
     do j=1, ny
        cur_pres(:)=pres(i,j,:)
        call get_col_wgt_mask(nz, cur_pres, cur_wgt, fmask)
        wgt(i,j,:)=cur_wgt(:)
     end do
  end do
  
  deallocate(cur_pres, cur_wgt)
  
end subroutine get_col_wgt_2d


subroutine col_int_0d(nz, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  integer, intent(in)::nz
  real*8, intent(in)::prof(nz)
  real*8, intent(in)::colwgt(nz)
  real*8, intent(out)::colval
  real*8, intent(in), optional::mask_val
  
  real*8::fmask
  
  ! local variable
  
  integer::i, k
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  colval=0
  do k=1, nz
     ! skip the masked value
     
     if (prof(k)/=fmask) & 
          colval=colval+colwgt(k)*prof(k)
  end do
  
end subroutine col_int_0d


subroutine col_int_1d(nx, nz, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  integer, intent(in)::nx, nz
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::colwgt(nx, nz)
  real*8, intent(out)::colval(nx)
  real*8, intent(in), optional::mask_val

  real*8::fmask
  
  ! local variable
  
  integer::i, k
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  do i=1, nx
     colval(i)=0
     do k=1, nz
        ! skip the masked value
        
        if (prof(i,k)/=fmask) & 
             colval(i)=colval(i)+colwgt(i,k)*prof(i, k)
     end do
  end do
  
end subroutine col_int_1d


subroutine col_int_2d(nx, ny, nz, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  integer, intent(in)::nx, ny ,nz
  real*8, intent(in)::prof(nx, ny, nz)
  real*8, intent(in)::colwgt(nx, ny, nz)
  real*8, intent(out)::colval(nx, ny)
  real*8, intent(in), optional::mask_val

  real*8::fmask
  
  ! local variable
  
  integer::i, j , k
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  do i=1, nx
     do j=1, ny
        colval(i, j)=0
        do k=1, nz
           ! skip the masked value
           if (prof(i,j, k)/=fmask) & 
                colval(i, j)=colval(i,j)+colwgt(i,j,k)*prof(i, j, k)
        end do
     end do
  end do
end subroutine col_int_2d


subroutine col_int_0d_em(nz, ne, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  integer, intent(in)::nz, ne
  real*8, intent(in)::prof(nz, ne)
  real*8, intent(in)::colwgt(nz)
  real*8, intent(out)::colval(ne)
  real*8, intent(in), optional::mask_val
  
  real*8::fmask
  
  ! local variable
  
  integer::i, k
  integer::ie
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  colval=0
  do ie=1, ne
     colval(ie)=0.0
     do k=1, nz
        ! skip the masked value
        if (prof(k,ie)/=fmask) & 
             colval(ie)=colval(ie)+colwgt(k)*prof(k, ie)
     end do
  end do

end subroutine col_int_0d_em

subroutine col_int_1d_em(nx, nz, ne, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  integer, intent(in)::nx, nz,ne
  real*8, intent(in)::prof(nx, nz, ne)
  real*8, intent(in)::colwgt(nx, nz)
  real*8, intent(out)::colval(nx, ne)
  real*8, intent(in), optional::mask_val
  
  real*8::fmask
  
  ! local variable
  
  integer::i, k, ie
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  do ie=1, ne
     do i=1, nx
        colval(i, ie)=0
        do k=1, nz
           ! skip the masked value
           
           if (prof(i,k,ie)/=fmask) & 
                colval(i, ie)=colval(i,ie)+colwgt(i,k)*prof(i, k,ie)
        end do
     end do
  end do
  
end subroutine col_int_1d_em


subroutine col_int_2d_em(nx, ny, nz, ne, prof, colwgt, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1 integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ration 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  
  
  implicit none
  
  integer, intent(in)::nx, ny, nz,ne
  real*8, intent(in)::prof(nx, ny, nz, ne)
  real*8, intent(in)::colwgt(nx, ny, nz)
  real*8, intent(out)::colval(nx, ny, ne)
  real*8, intent(in), optional::mask_val
  
  real*8::fmask
  
  ! local variable
  
  integer::ie, i,j , k
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  do ie=1, ne
     do i=1, nx
        do j=1, ny
        colval(i, j, ie)=0
        do k=1, nz
           ! skip the masked value
           if (prof(i,j, k, ie)/=fmask) & 
                colval(i, j,ie)=colval(i, j,ie)+colwgt(i,j, k)*prof(i, j, k, ie)
        end do
     end do
  end do
end do

end subroutine col_int_2d_em



subroutine ak_col_int_0d(nz, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ratios 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  ! 5. real*8    :: ak(nx, nz) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nz
  real*8, intent(in)::prof(nz)
  real*8, intent(in)::ak(nz)
  real*8, intent(in)::colwgt(nz)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::colval
  integer::i, j, k
  
  real*8::fmask

  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  colval=0
  do k=1, nz
     if ((ak(k)/=fmask).and.(prof(k)/=fmask)) & 
          colval=colval+colwgt(k)*ak(k)*prof(k)
  end do
  
end subroutine ak_col_int_0d


subroutine ak_col_int_1d(nx, nz, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ratios 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  ! 5. real*8    :: ak(nx, nz) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nx, nz
  real*8, intent(in)::prof(nx, nz)
  real*8, intent(in)::ak(nx, nz)
  real*8, intent(in)::colwgt(nx, nz)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::colval(nx)
  integer::i, j, k
  
  real*8::fmask

  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  do i=1, nx
     colval(i)=0
     do k=1, nz
        if ((ak(i,k)/=fmask).and.(prof(i,k)/=fmask)) & 
             colval(i)=colval(i)+colwgt(i,k)*ak(i,k)*prof(i, k)
     end do
  end do
end subroutine ak_col_int_1d




subroutine ak_col_int_2d(nx, ny, nz, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ratios 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  ! 5. real*8    :: ak(nx, nz) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nx, ny, nz
  real*8, intent(in)::prof(nx, ny, nz)
  real*8, intent(in)::ak(nx, ny, nz)
  real*8, intent(in)::colwgt(nx, ny,nz)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::colval(nx, ny)
  integer::i, j, k
  
  real*8::fmask

  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  do i=1, nx
     do j=1, ny
        colval(i, j)=0
        do k=1, nz
           if ((ak(i,j,k)/=fmask).and.(prof(i,j, k)/=fmask)) & 
                colval(i,j)=colval(i,j)+colwgt(i,j,k)*ak(i,j,k)*prof(i,j, k)
        end do
     end do
  end do
end subroutine ak_col_int_2d


subroutine ak_col_int_0d_em(nz, ne, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles of ne tracers at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. integer*4 ::ne #  number of tracers 
  ! 4. real*8    :: prof(nx, ny, ne) # mixing ratios 
  ! 5. real*8    :: colwgt(nx, ny) # normalized airmass weighting factor 
  ! 6. real*8    :: ak(nx, ny) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nz, ne
  real*8, intent(in)::prof(nz,ne)
  real*8, intent(in)::ak(nz)
  real*8, intent(in)::colwgt(nz)

  real*8, intent(out)::colval(ne)
  
  real*8, intent(in), optional::mask_val

  integer::i, ie, k

  real*8::fmask
  
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  
  do ie=1, ne
     colval(ie)=0
     do k=1, nz
        if ((ak(k)/=fmask).and.(prof(k,ie)/=fmask)) & 
             colval(ie)=colval(ie)+colwgt(k)*ak(k)*prof(k,ie)
     end do
  end do

  
end subroutine ak_col_int_0d_em


subroutine ak_col_int_1d_em(nx, nz, ne, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles of ne tracers at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. integer*4 ::ne #  number of tracers 
  ! 4. real*8    :: prof(nx, ny, ne) # mixing ratios 
  ! 5. real*8    :: colwgt(nx, ny) # normalized airmass weighting factor 
  ! 6. real*8    :: ak(nx, ny) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nx, nz, ne
  real*8, intent(in)::prof(nx, nz,ne)
  real*8, intent(in)::ak(nx, nz)
  real*8, intent(in)::colwgt(nx, nz)

  real*8, intent(out)::colval(nx,ne)
  
  real*8, intent(in), optional::mask_val

  integer::i, ie, k

  real*8::fmask
  
  
  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  

  do ie=1, ne
     do i=1, nx
        colval(i, ie)=0
        do k=1, nz
           if ((ak(i,k)/=fmask).and.(prof(i, k,ie)/=fmask)) & 
                colval(i,ie)=colval(i,ie)+colwgt(i,k)*ak(i,k)*prof(i, k,ie)
        end do
     end do
  end do
  
end subroutine ak_col_int_1d_em


subroutine ak_col_int_2d_em(nx, ny, nz, ne, prof, colwgt, ak, colval, mask_val)
  
  ! Calculate the column values for mixing ratio  profiles at nx points, using mass weighting factor
  ! Inputs:  
  ! 1. integer*4 ::nx #  number of points = nx
  ! 2. integer*4 ::nz #  number of levels 
  ! 3. real*8    :: prof(nx, nz) # mixing ratios 
  ! 4. real*8    :: colwgt(nx, nz) # normalized airmass weighting factor 
  ! 5. real*8    :: ak(nx, nz) # normalized averaging kernels
  
  
  !
  ! Outputs
  ! 1. real*8 :: colval(nx): column mixing ratio
  


  implicit none
  integer, intent(in)::nx, ny, nz, ne
  real*8, intent(in)::prof(nx, ny, nz, ne)
  real*8, intent(in)::ak(nx, ny, nz)
  real*8, intent(in)::colwgt(nx, ny,nz)
  real*8, intent(in), optional::mask_val
  
  real*8, intent(out)::colval(nx, ny, ne)
  integer::i, j, k, ie
  
  real*8::fmask

  if (present(mask_val)) then 
     fmask=mask_val
  else
     fmask=-999.0
  end if
  
  
  do i=1, nx
     do j=1, ny
        colval(i, j,ie)=0
        do k=1, nz
           do ie=1, ne
              if ((ak(i,j,k)/=fmask).and.(prof(i,j, k,ne)/=fmask)) & 
                   colval(i,j,ie)=colval(i,j,ie)+colwgt(i,j,k)*ak(i,j,k)*prof(i,j, k,ie)
           end do
        end do
     end do
  end do
end subroutine ak_col_int_2d_em
