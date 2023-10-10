subroutine get_bs_wgt(x, nx, x2, nx2)
  ! get the weighting factor  for flux interpolation from small box to large box 
  ! 
  !  
  
  integer::intent(in)::nx
  integer::intent(in)::nx2
  real*8, intent(in)::x(nx)
  real*8, intent(in)::x2(nx2)
  
  pst=1


! the left is the first one larger than 
! the right one >1

do ix=1, nlon2-1
   plon(ix)=1
   do ilon=pst, nlon-1
      rlon=lon1(ilon+1)
   
      if (rlon>=lon2(ix+1)) then
         plon(ix)=ilon
         hflon=0.5*(lon1(ilon+1)+lon1(ilon))
         
         wlon(ix)=(lon2(ix+1)-hflon)/(lon2(ix+1)-lon2(ix))
         
         if (wlon(ix)>1) then 
            wlon(ix)=0
         else if (wlon(ix)<0) then
            wlon(ix)=1
         else
            wlon(ix)=1.0-wlon(ix)
         end if

         pst=ilon
         exit
      end if
   end do
end do

plon(nlon2)=nlon-1
wlon(nlon2)=1.0

end subroutine get_bs_wgt


subroutine regrid_flux_bs(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, flux1, flux2)

  !  regrid the flux from larger grid boxes to smaller ones 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   edges of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, edges of target grid box
  !  flux1    in  array of (nlon, nlat)  origin surface flux,which is assumed to be at the center of the box 
  !  flux2    out array of (nlon2, nlat2) surface flux at target  grid box   
  !  notes:
  !  we assume the edge box of longitude will not be used. but the latitude one will be used
  ! 
  
implicit none
integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
real*8, intent(in)::flux1(nlon-1, nlat-1)


integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
real*8, intent(out)::flux2(nlon2-1, nlat2-1)


integer, allocatable:: plon(:), plat(:), used(:,:)

real*8 :: w1, w2, w3, w4, total_w
real*8 :: flux_11, flux_12, flux_21, flux_22
integer::pst, ix, iy, xst, yst, xend, yend
integer::latp1, latp2, lonp1, lonp2

real*8::fll, ful, flu, fuu, hflon, hflat



real*8, allocatable::wlon(:)
real*8, allocatable::wlat(:)




real*8::rlon, rlat, rlonL, rlonR, rlatL, rlatR


integer::ilat, ilon, llon, ulon, llat, ulat 
! the left ledge 

allocate(plon(nlon2))
allocate(wlon(nlon2))

allocate(plat(nlat2))
allocate(wlat(nlat2))


! search the lon1
 
flux2=0.0




do ix=1, nlat2-1
   plat(ix)=1
   do ilat=pst, nlat-1
      rlat=lat1(ilat+1)
      if (rlat>=lat2(ix+1)) then
         
         hflat=0.5*(lon1(ilat+1)+lon1(ilat))
         
         wlat(ix)=(lat2(ix+1)-hflat)/(lat2(ix+1)-lat2(ix))
         
         if (wlat(ix)>1) then 
            wlat(ix)=0
         else if (wlat(ix)<0) then
            wlat(ix)=1
         else
            wlat(ix)=1.0-wlat(ix)
         end if


         plat(ix)=ilat
         pst=ilat
         exit
      end if
   end do
end do

plat(nlat2)=nlat-1
wlat(nlat2)=1.0


! fill the grid box 

do ix=1, nlon2-1
   llon=plon(ix)
   ulon=llon+1
   if (ulon>nlon-1) ulon=nlon-1
   
   do iy=1, nlat2-1
      llat=plat(iy)
      ulat=llat+1
      if (ulat>nlat-1) ulat=nlat-1
      
      fll=wlon(ix)*wlat(iy)*flux1(llon, llat)
      ful=(1.0-wlon(ix))*wlat(iy)*flux1(ulon, llat)
      flu=wlon(ix)*(1.0-wlat(iy))*flux1(llon, ulat)
      fuu=(1.0-wlon(ix))*(1.0-wlat(iy))*flux1(ulon, ulat)
      flux2(ix, iy)=fll+ful+flu+fuu
   end do
end do


deallocate(plon)
deallocate(wlon)

deallocate(plat)
deallocate(wlat)

end subroutine regrid_flux_bs



subroutine regrid_map_bs(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, flux1, flux2)

  !  regrid the flux from larger grid boxes to smaller ones 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   edges of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, edges of target grid box
  !  flux1    in  array of (nlon, nlat)  origin surface flux,which is assumed to be at the center of the box 
  !  flux2    out array of (nlon2, nlat2) surface flux at target  grid box   
  !  notes:
  !  we assume the edge box of longitude will not be used. but the latitude one will be used
  ! 
  
implicit none
integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
integer, intent(in)::flux1(nlon-1, nlat-1)


integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
real*8, intent(out)::flux2(nlon2-1, nlat2-1)


integer, allocatable:: plon(:), plat(:), used(:,:)

real*8 :: w1, w2, w3, w4, total_w
real*8 :: flux_11, flux_12, flux_21, flux_22
integer::pst, ix, iy, xst, yst, xend, yend
integer::latp1, latp2, lonp1, lonp2

real*8::fll, ful, flu, fuu, hflon, hflat



real*8, allocatable::wlon(:)
real*8, allocatable::wlat(:)




real*8::rlon, rlat, rlonL, rlonR, rlatL, rlatR


integer::ilat, ilon, llon, ulon, llat, ulat 
! the left ledge 

allocate(plon(nlon2))
allocate(wlon(nlon2))

allocate(plat(nlat2))
allocate(wlat(nlat2))


! search the lon1
 
flux2=0.0

pst=1


! the left is the first one larger than 
! the right one >1

do ix=1, nlon2-1
   plon(ix)=1
   do ilon=pst, nlon-1
      rlon=lon1(ilon+1)
   
      if (rlon>=lon2(ix+1)) then
         plon(ix)=ilon
         hflon=0.5*(lon1(ilon+1)+lon1(ilon))
         
         wlon(ix)=(lon2(ix+1)-hflon)/(lon2(ix+1)-lon2(ix))
         
         if (wlon(ix)>1) then 
            wlon(ix)=0
         else if (wlon(ix)<0) then
            wlon(ix)=1
         else
            wlon(ix)=1.0-wlon(ix)
      
            if (wlon(ix)>0.5) then 
               wlon(ix)=1.0
            else
               wlon(ix)=0.0
            end if

         end if

         pst=ilon
         exit
      end if
   end do
end do

plon(nlon2)=nlon-1
wlon(nlon2)=1.0



do ix=1, nlat2-1
   plat(ix)=1
   do ilat=pst, nlat-1
      rlat=lat1(ilat+1)
      if (rlat>=lat2(ix+1)) then
         
         hflat=0.5*(lon1(ilat+1)+lon1(ilat))
         
         wlat(ix)=(lat2(ix+1)-hflat)/(lat2(ix+1)-lat2(ix))
         
         if (wlat(ix)>1) then 
            wlat(ix)=0
         else if (wlat(ix)<0) then
            wlat(ix)=1
         else
            wlat(ix)=1.0-wlat(ix)
            if (wlon(ix)>0.5) then 
               wlon(ix)=1.0
            else
               wlon(ix)=0.0
            end if
         end if
         

         plat(ix)=ilat
         pst=ilat
         exit
      end if
   end do
end do

plat(nlat2)=nlat-1
wlat(nlat2)=1.0


! fill the grid box 

do ix=1, nlon2-1
   llon=plon(ix)
   ulon=llon+1
   if (ulon>nlon-1) ulon=nlon-1
   
   do iy=1, nlat2-1
      llat=plat(iy)
      ulat=llat+1
      if (ulat>nlat-1) ulat=nlat-1
      
      fll=wlon(ix)*wlat(iy)*flux1(llon, llat)
      ful=(1.0-wlon(ix))*wlat(iy)*flux1(ulon, llat)
      flu=wlon(ix)*(1.0-wlat(iy))*flux1(llon, ulat)
      fuu=(1.0-wlon(ix))*(1.0-wlat(iy))*flux1(ulon, ulat)
      flux2(ix, iy)=fll+ful+flu+fuu
   end do
end do


deallocate(plon)
deallocate(wlon)

deallocate(plat)
deallocate(wlat)

end subroutine regrid_map_bs

subroutine regrid_flux(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, flux1, inmap, &
     flux2, outmap)

  !  regrid the flux from grid boxes of 1x1 to larger grid boxes 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   centre of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, centre of target grid box
  !  inmap    in  array of (nlon, nlat), surface type  of origin grid box
  !  outmap   in  array of (nlon2, nlat2), surface type  of target grid box 
  !  flux1    in  array of (nlon, nlat)  origin surface flux
  !  flux2    out array of (nlon2, nlat2) surface flux at target  grid box   
  
implicit none
integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
real*8, intent(in)::flux1(nlon, nlat)
integer, intent(in)::inmap(nlon, nlat)


integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
real*8, intent(out)::flux2(nlon2, nlat2)
integer, intent(in)::outmap(nlon2, nlat2)


integer, allocatable:: plon(:), plat(:), used(:,:)
real*8 :: w1, w2, w3, w4, total_w
real*8 :: flux_11, flux_12, flux_21, flux_22
integer::pst, ix, iy, xst, yst, xend, yend
integer::latp1, latp2, lonp1, lonp2
real*8::lon_width, lat_width


integer::lon_grd(3), lat_grd(3)


real*8::rlon, rlat


integer::ilat, ilon 

allocate(plon(nlon))
allocate(plat(nlat))
allocate(used(nlon, nlat))


! search the lon1
 
used=0
flux2=0.0

pst=1



do ilon=1, nlon
   rlon=lon1(ilon)
   plon(ilon)=-999
   
   do ix=pst, nlon2
      
      if (rlon<=lon2(ix)) then
         if (ix==1) then 
            plon(ilon)=1
            pst=1
         else  
            plon(ilon)=ix-1
            pst=ix-1
         end if
         
         exit
         
      end if
   end do
   if (plon(ilon)<0) then 
      plon(ilon)=nlon2
   end if

end do



! search the closest lat 

pst=1

do ilat=1, nlat
   rlat=lat1(ilat)
   plat(ilat)=-999
   
   do iy=pst, nlat2
      
      if (rlat<=lat2(iy)) then
         if (iy==1) then 
            plat(ilat)=1
            pst=1
         else  
            plat(ilat)=iy-1
            pst=iy-1
         end if
         
         exit
         
      end if
   end do
   if (plat(ilat)<0) then 
      plat(ilat)=nlat2
   end if
end do

!  width of old grid box 
 
lon_width=lon1(3)-lon1(2)
lat_width=lat1(3)-lat1(2)





! assign flux to the new grid boxs  


do ilon=1, nlon
   lonp1=plon(ilon)
   lonp2=plon(ilon)+1
   if (lonp2>nlon2) then 
      lonp2=1
      w1=(0.5*(lon2(lonp1)+360.0+lon2(lonp2))-lon1(ilon))/lon_width
      lonp2=nlon2
      w1=(0.5*(lon2(lonp1)+lon2(lonp1))-lon1(ilon))/lon_width
   else
      w1=(0.5*(lon2(lonp1)+lon2(lonp2))-lon1(ilon))/lon_width
      
   end if
   
   if (w1>=1.0) then 
      w1=1.0
   else if (w1<=-1.0) then 
      w1=0.0
   else if (w1<0.0) then 
      w1=abs(w1)
   end if
   
   w2=1.0-w1
      
   
   do ilat=1, nlat
   
      
      latp1=plat(ilat)
      latp2=latp1+1
      
      if (latp2>nlat2) latp2=nlat2
      
      
      w3=(0.5*(lat2(latp1)+lat2(latp2))-lat1(ilat))/lat_width
      
      if (w3>=1.0) then 
         w3=1.0
      else if (w3<=-1.0) then 
         w3=0.0
      else if (w3<0.0) then 
         w3=abs(w3)
      end if
      
      w4=1.0-w3
      
      
      
      total_w=0.0
      
      flux_11=0.0
      flux_12=0.0
      flux_21=0.0
      flux_22=0.0
      
      if (inmap(ilon,ilat)==outmap(lonp1, latp1)) then 
         total_w=total_w+w1*w3
         flux_11=w1*w3*flux1(ilon, ilat)
      end if
      
      if (inmap(ilon,ilat)==outmap(lonp2, latp1)) then 
         total_w=total_w+w2*w3
         flux_21=w2*w3*flux1(ilon, ilat)
      end if
      
      if (inmap(ilon,ilat)==outmap(lonp1, latp2)) then 
         total_w=total_w+w1*w4
         flux_12=w1*w4*flux1(ilon, ilat)
      end if
      
      
      if (inmap(ilon,ilat)==outmap(lonp2, latp2)) then 
         total_w=total_w+w2*w4
         flux_22=w2*w4*flux1(ilon, ilat)
      end if
      ! added the flux to grid box 
      
      if (total_w>0) then 
      
         flux2(lonp1, latp1)=flux2(lonp1, latp1)+flux_11/total_w
         flux2(lonp2, latp1)=flux2(lonp2, latp1)+flux_21/total_w
         flux2(lonp1, latp2)=flux2(lonp1, latp2)+flux_12/total_w
         flux2(lonp2, latp2)=flux2(lonp2, latp2)+flux_22/total_w

         used(ilon,ilat)=1.0
         
      end if
      
      
      
      
      ! check whether it is fully contain 
      
   end do
   
end do





! handle the miss-match points grid points, and try to assign them to left and south one. 




do ilon=1, nlon
   do ilat=1, nlat
      if (used(ilon, ilat)==0) then 
         xst=1
         xend=3
         
         lon_grd(1)=plon(ilon)-1
         if (lon_grd(1)<1) then 
            lon_grd(1)=1
            xst=2
         end if
         
         lon_grd(2)=plon(ilon)
         lon_grd(3)=plon(ilon)+1
         
         if (lon_grd(3)>nlon2) then 
            ! goes back the first point 
            
            lon_grd(3)=1 
            
         end if

         yst=1
         yend=3
         
         lat_grd(1)=plat(ilat)-1
         
         if (lat_grd(1)<1) then 
            lat_grd(1)=1
            yst=2
         end if
         
         lat_grd(2)=plat(ilat)
         lat_grd(3)=plat(ilat)+1
         
         if (lat_grd(3)>nlat2) then 
            lat_grd(3)=nlat
            yend=2
         end if
         
         total_w=0.0
         do ix=xst, xend
            do iy=yst,yend
               if (inmap(ilon, ilat)==outmap(lon_grd(ix), lat_grd(iy))) then
                  total_w=total_w+abs(flux2(lon_grd(ix), lat_grd(iy)))
               end if
               
               
            end do
         end do
         
         
         
         if (total_w>0) then 
            do ix=xst, xend
               do iy=yst,yend
                  if (inmap(ilon, ilat)==outmap(lon_grd(ix), lat_grd(iy))) then
                     flux2(lon_grd(ix), lat_grd(iy))=flux2(lon_grd(ix), lat_grd(iy))+& 
                          flux1(ilon, ilat)* abs(flux2(lon_grd(ix), lat_grd(iy)))/total_w
                  end if
                  
               end do
            end do
         end if
         
         
      end if
      
   end do
   
end do
deallocate(plon)
deallocate(plat)
deallocate(used)


end subroutine regrid_flux

  
subroutine regrid_map(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, inmap, landbdr, landmap, outmap)
  
  !  regrid the flux from grid boxes of 1x1 to larger grid boxes 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   centre of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, centre of target grid box
  !  inmap    in  array of (nlon, nlat), surface type  of origin grid box
  !  outmap   in  array of (nlon2, nlat2), surface type  of target grid box 
  !  flux1    in  array of (nlon, nlat)  origin surface flux
  !  flux2    out array of (nlon2, nlat2) surface flux at target  grid box   
  
implicit none

integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
integer, intent(in)::inmap(nlon, nlat)
integer, intent(in)::landbdr

integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
integer, intent(in)::landmap(nlon2, nlat2)
integer, intent(out)::outmap(nlon2, nlat2)


integer, allocatable:: plon(:), plat(:), used(:,:)
real*8 :: w1, w2, w3, w4
real*8, allocatable::reg_cnt(:,:,:)

integer::pst, sel_ix, sel_iy, ix, iy, iz


integer::latp1, latp2, lonp1, lonp2
real*8::lon_width, lat_width
integer::maxreg, mk
integer::regix(1)


real*8::rlon, rlat


integer::ilat, ilon 

allocate(plon(nlon))
allocate(plat(nlat))
maxreg=maxval(inmap)
allocate(reg_cnt(nlon, nlat,maxreg+1))


! search the locations from  lon2 to encircle lon1
 

pst=1

do ilon=1, nlon
   rlon=lon1(ilon)
   plon(ilon)=-999
   
   do ix=pst, nlon2
      
      if (rlon<=lon2(ix)) then
         if (ix==1) then 
            plon(ilon)=1
            pst=1
         else  
            plon(ilon)=ix-1
            pst=ix-1
         end if
         
         exit
         
      end if
   end do
   if (plon(ilon)<0) then 
      plon(ilon)=nlon2
   end if

end do



! search the closest location and weight 

pst=1

do ilat=1, nlat
   rlat=lat1(ilat)
   plat(ilat)=-999
   
   do iy=pst, nlat2
      
      if (rlat<=lat2(iy)) then
         if (iy==1) then 
            plat(ilat)=1
            pst=1
         else  
            plat(ilat)=iy-1
            pst=iy-1
         end if
         
         exit
         
      end if
   end do
   if (plat(ilat)<0) then 
      plat(ilat)=nlat2
   end if
end do

!  width of old grid box 
 
lon_width=lon1(3)-lon1(2)
lat_width=lat1(3)-lat1(2)





! assign flux to grid box 


do ilon=1, nlon
   lonp1=plon(ilon)
   lonp2=plon(ilon)+1
   if (lonp2>nlon2) then 
      lonp2=1
      w1=(0.5*(lon2(lonp1)+360.0+lon2(lonp2))-lon1(ilon))/lon_width
   else
      w1=(0.5*(lon2(lonp1)+lon2(lonp2))-lon1(ilon))/lon_width
   end if
   
   if (w1>=1.0) then 
      w1=1.0
   else if (w1<=-1.0) then 
      w1=0.0
   else if (w1<0.0) then 
      w1=abs(w1)
   end if
   
   w2=1.0-w1
      
   
   do ilat=1, nlat
   
      
      latp1=plat(ilat)
      latp2=latp1+1
      
      if (latp2>nlat2) latp2=nlat2
      
      
      w3=(0.5*(lat2(latp1)+lat2(latp2))-lat1(ilat))/lat_width
      
      if (w3>=1.0) then 
         w3=1.0
      else if (w3<=-1.0) then 
         w3=0.0
      else if (w3<0.0) then 
         w3=abs(w3)
      end if
      
      w4=1.0-w3
      
      iz=inmap(ilon, ilat)+1
      reg_cnt(lonp1, latp1, iz)=reg_cnt(lonp1, latp1, iz)+w1*w3
      reg_cnt(lonp2, latp1, iz)=reg_cnt(lonp2, latp1, iz)+w2*w3
      reg_cnt(lonp1, latp2, iz)=reg_cnt(lonp1, latp2, iz)+w1*w4
      reg_cnt(lonp2, latp2, iz)=reg_cnt(lonp2, latp2, iz)+w2*w4
      
            
   end do
end do






! handle the miss-match points grid points, and try to assign them to left and south one. 

do ilon=1, nlon2
   do ilat=1, nlat2
      regix=maxloc(reg_cnt(ilon, ilat,:))
      outmap(ilon, ilat)=regix(1)-1
      
      do mk=1,maxreg 
         if (landmap(ilon ,ilat)==1) then 
            if (regix(1)<=landbdr+1) then 
               outmap(ilon, ilat)=regix(1)-1
               exit
            else
               reg_cnt(ilon, ilat,regix(1))=0.0
            end if
            
         else if (landmap(ilon ,ilat)==0) then
            
            if (regix(1)>landbdr+1) then 
               outmap(ilon, ilat)=regix(1)-1
               exit
            else
               reg_cnt(ilon, ilat,regix(1))=0.0
            end if
         end if
         regix=maxloc(reg_cnt(ilon, ilat,:))
      end do
      
   end do
end do

deallocate(reg_cnt)
deallocate(plon)
deallocate(plat)

end subroutine regrid_map


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
  ! real*8, intent(out)::outmap(nlon2, nlat2, ml)

  
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
  
  integer, intent(in)::nlon, nlat ! longitude and latitude size
  real*8, intent(in)::rlon1(nlon), rlat1(nlat)  ! original longitude and latitude grid
  integer, intent(in)::inmap(nlon, nlat) !  regional map
  integer, intent(in)::ml  ! ml=max(inmap) is the layer number of the map 
  real*8, intent(in)::influx(nlon, nlat) !  emissions from each grid box (not per area)
  
  
  integer, intent(in)::nlon2, nlat2
  real*8, intent(in)::rlon2(nlon2), rlat2(nlat2)
  
  ! outs
  real*8, intent(out)::outmap(nlon2, nlat2, ml)
  real*8, intent(out)::outflux(nlon2, nlat2, ml)
  
  
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






subroutine regrid_flux_ml(lon1, lat1, nlon, nlat,  ml, &
     lon2, lat2, nlon2, nlat2, flux1,flux2)
  
  !  regrid the flux from grid boxes of 1x1 to larger grid boxes. The original region will become mutiple layer. 
  ! 
  !  we need to check the map to make sure the surface type is co
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   centre of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, centre of target grid box
  !  inmap    in  array of (nlon, nlat), surface type  at origin grid box
  !  outmap   in  array of (nlon2, nlat2, ml), surface type  at target grid box 

  
implicit none
integer, intent(in)::nlon, nlat, ml
real*8, intent(in)::lon1(nlon), lat1(nlat)
real*8, intent(in)::flux1(nlon, nlat, ml)


integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
real*8, intent(out)::flux2(nlon2, nlat2, ml)


integer, allocatable:: plon(:), plat(:), used(:,:)

real*8 :: w1, w2, w3, w4, total_w
real*8 :: flux_11, flux_12, flux_21, flux_22
integer::pst, ix, iy, xst, yst, xend, yend
integer::latp1, latp2, lonp1, lonp2
real*8::lon_width, lat_width
integer::ilvl


integer::lon_grd(3), lat_grd(3)


real*8::rlon, rlat, rlonL, rlonR, rlatL, rlatR


integer::ilat, ilon, ilonL, iLonR, ilatL, ilatR 

allocate(plon(nlon))
allocate(plat(nlat))


! search the lon1
 

flux2=0.0

pst=1




do ilon=1, nlon
   rlon=lon1(ilon)
   plon(ilon)=-999
   
   do ix=pst, nlon2
      
      if (rlon<=lon2(ix)) then
         if (ix==1) then 
            plon(ilon)=1
            pst=1
         else  
            plon(ilon)=ix-1
            pst=ix-1
         end if
         
         exit
         
      end if
   end do
   
   if (plon(ilon)<0) then 
      ! print *, ix, rlon, lon2(nlon2)
      plon(ilon)=nlon2
   end if
   
end do



! search the closest lat 

pst=1

do ilat=1, nlat
   rlat=lat1(ilat)
   plat(ilat)=-999
   
   do iy=pst, nlat2
      
      if (rlat<=lat2(iy)) then
         if (iy==1) then 
            plat(ilat)=1
            pst=1
         else  
            plat(ilat)=iy-1
            pst=iy-1
         end if
         
         exit
         
      end if
   end do
   if (plat(ilat)<0) then 
      plat(ilat)=nlat2
   end if
end do

!  width of old grid box 
 
lon_width=0.5*(lon1(3)-lon1(2))

lat_width=0.5*(lat1(3)-lat1(2))






! assign flux to the new grid boxs  


do ilon=1, nlon
   ! choose the point
   
   rlon=lon1(ilon)
   lonp1=plon(ilon)
   ilonL=plon(ilon)-1
   
   if (ilonL<1) then
      ilonL=nlon2
      rlonL=0.5*(lon2(lonp1)+lon2(ilonL)-360.0)
   else
      rlonL=0.5*(lon2(lonp1)+lon2(ilonL))
   end if
   
   ilonR=plon(ilon)+1
   
   if (ilonR>nlon2) ilonR=1
   
   rlonR=rlonL+lon2(3)-lon2(2)
   
   
   
   if (rlonL>(rlon-lon_width)) then 
      lonp2=ilonL
      w1=1.0-0.5*((rlon-lon_width-rLonL)/lon_width)
      if (w1>1.0) w1=1.0
      if (w1<0.0) w1=0.0
      
   else if (rlonR<(rlon+lon_width)) then
      lonp2=ilonR
      w1=1.0-0.5*((rlon+lon_width-rlonR)/lon_width)
      if (w1>1.0) w1=1.0
      if (w1<0.0) w1=0.0
   else
      lonp2=ilonR
      w1=1.0
   end if
   w1=1.0
   
   
   w2=1.0-w1
   
   
   do ilat=1, nlat
      rlat=lat1(ilat)
      latp1=plat(ilat)
      ilatL=plat(ilat)-1
      if (ilatL<1) then
         ilatL=1
         rlatL=0.5*(lat2(latp1)+lat2(ilatL))
      else
         rlatL=0.5*(lat2(latp1)+lat2(ilatL))
      end if
   
      ilatR=plat(ilat)+1
      
      if (ilatR>nlat2) ilatR=nlat2
      
      rlatR=0.5*(lat2(latp1)+lat2(ilatR))
   
      
      if (rlatL>(rlat-lat_width)) then 
         latp2=ilatL
         w3=1.0-0.5*(rlat-lat_width-rLatL)/lat_width
         if (w3>1.0) w3=1.0
         if (w3<1.0) w3=0.0
      
      else if (rlatR<(rlat+lat_width)) then
         latp2=ilatR
         w3=1.0-0.5*(rlat+lat_width-rlatR)/lat_width
         if (w3>1.0) w3=1.0
         if (w3<0.0) w3=0.0
      
      else
         latp2=ilatR
         w3=1.0
      end if
      w3=1.0
      
      w4=1.0-w3
      
      total_w=0.0
      total_w=total_w+w1*w3
      total_w=total_w+w2*w3
      total_w=total_w+w1*w4
      total_w=total_w+w2*w4
      
      do ilvl=1, ml
         flux_11=w1*w3*flux1(ilon, ilat, ilvl)
         flux_21=w2*w3*flux1(ilon, ilat, ilvl)
         flux_12=w1*w4*flux1(ilon, ilat, ilvl)
         flux_22=w2*w4*flux1(ilon, ilat, ilvl)
         
         flux2(lonp1, latp1, ilvl)=flux2(lonp1, latp1, ilvl)+flux_11/total_w
         flux2(lonp2, latp1, ilvl)=flux2(lonp2, latp1, ilvl)+flux_21/total_w
         flux2(lonp1, latp2, ilvl)=flux2(lonp1, latp2, ilvl)+flux_12/total_w
         flux2(lonp2, latp2, ilvl)=flux2(lonp2, latp2, ilvl)+flux_22/total_w
         
      end do
   end do
end do

! handle the miss-match points grid points, and try to assign them to left and south one. 

deallocate(plon)
deallocate(plat)


end subroutine regrid_flux_ml


subroutine extract_region_flux_ml(lon1, lat1, nlon, nlat, &
     lon2, lat2, nlon2, nlat2, inmap, influx, ml, outmap, outflux)
  
  !  regrid the flux from grid boxes of 1x1 to larger grid boxes. The original region will become mutiple layer. 
  ! 
  !  we need to check the map to make sure the surface type is co
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   centre of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, centre of target grid box
  !  inmap    in  array of (nlon, nlat), surface type  at origin grid box
  !  outmap   in  array of (nlon2, nlat2, ml), surface type  at target grid box 

  
implicit none

integer, intent(in)::nlon, nlat
real*8, intent(in)::lon1(nlon), lat1(nlat)
integer, intent(in)::inmap(nlon, nlat), influx(nlon, nlat)

integer, intent(in)::ml

integer, intent(in)::nlon2, nlat2
real*8, intent(in)::lon2(nlon2), lat2(nlat2)
real*8, intent(out)::outmap(nlon2, nlat2, ml), outflux(nlon2, nlat2, ml)



integer, allocatable:: plon(:), plat(:)
real*8 :: wgt
real*8, allocatable::reg_cnt(:,:,:)

integer::pst, ix, iy, iz


integer::latp1, latp2, lonp1, lonp2
real*8::lon_width, lat_width
integer::mk



real*8::rlon, rlat


integer::ilat, ilon 

allocate(plon(nlon))
allocate(plat(nlat))
allocate(reg_cnt(nlon, nlat,ml))

reg_cnt=0.0
outmap=0.0

! search the locations from  lon2 to encircle lon1
 

pst=1

do ilon=1, nlon
   rlon=lon1(ilon)
   plon(ilon)=-999
   
   do ix=pst, nlon2
      
      if (rlon<=lon2(ix)) then
         if (ix==1) then 
            plon(ilon)=1
            pst=1
         else  
            plon(ilon)=ix-1
            pst=ix-1
         end if
         
         exit
         
      end if
   end do
   if (plon(ilon)<0) then 
      plon(ilon)=nlon2
   end if

end do



! search the closest location and weight 

pst=1

do ilat=1, nlat
   rlat=lat1(ilat)
   plat(ilat)=-999
   
   do iy=pst, nlat2
      
      if (rlat<=lat2(iy)) then
         if (iy==1) then 
            plat(ilat)=1
            pst=1
         else  
            plat(ilat)=iy-1
            pst=iy-1
         end if
         
         exit
         
      end if
   end do
   if (plat(ilat)<0) then 
      plat(ilat)=nlat2
   end if
end do

!  width of old grid box 
 
lon_width=lon1(3)-lon1(2)
lat_width=lat1(3)-lat1(2)

! assign flux to grid box 

do ilon=1, nlon
   lonp1=plon(ilon)
   do ilat=1, nlat
      latp1=plat(ilat)
      iz=inmap(ilon, ilat)+1
      reg_cnt(lonp1, latp1, iz)=reg_cnt(lonp1, latp1, iz)+1
      outflux(lonp1, latp1, iz)= outflux(lonp1, latp1, iz)+influx(lonp1, latp1)
   end do
end do

! handle the miss-match points grid points, and try to assign them to left and south one. 

do ilon=1, nlon2
   do ilat=1, nlat2
      wgt=sum(reg_cnt(ilon, ilat, :))
      do iz=1, ml
         if (reg_cnt(ilon, ilat, iz)>0) then
            outflux(ilon, ilat, iz)=outflux(ilon, ilat, iz)/reg_cnt(ilon, ilat, iz)
            outmap(ilon, ilat,iz)=reg_cnt(ilon, ilat, iz)/wgt
         end if
      end do
   end do
end do

deallocate(reg_cnt)
deallocate(plon)
deallocate(plat)

end subroutine extract_region_flux_ml

subroutine get_upper_index(nx0,x0, nx1, x1, pos)

  ! find the position of the upper boundary in axis \
  ! x0 for each point of axis x1 
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


subroutine get_closest_point(nx0,x0, nx1, x1, pos)

  ! find the position in axis x0 closest to point of axis x1 
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
  ! 
  ! here 

  
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








subroutine refine_map(lon1, lat1, nlon, nlat, & 
     lon2, lat2, nlon2, nlat2, inmap, outmap)
  
  !  regrid map for coarse grid to a finer grid 
  !  
  !  Inputs:
  !-----------------------------------------------------------------
  !  integer:: nlon, nlat # size of original longitude and latitude 
  !  real*8:: lon1, lat1  # origin longitude and latitude 
  !  integer:: nlon2, nlat2 # size of target longitude and latitude 
  !  real*8,  lon2 & lat2   # target longitudes and latitudes 
  !  real*8:: inmap(nlon, nlat) #  original region map
  
  !  Outputs
  !--------------------------------------------------------------------
  !  outmap   in  array of (nlon2, nlat2), surface type  of target grid box 
  !  flux1    in  array of (nlon, nlat)  origin surface flux
  !  flux2    out array of (nlon2, nlat2) surface flux at target  grid box   
  
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




subroutine get_area_regular_grid(nlon, nlat, area)

  implicit none
#include "CMN_GCTM"  ! Physical constants
  
  integer, intent(in)::nlon, nlat
  real*8, intent(out)::area(nlon, nlat)
  real*8, allocatable::rlon_edge(:), rlat_edge(:)
  real*8 ::            D_2_R, factor
  integer::ix, iy
  real*8::dlon, dlat

  allocate(rlon_edge(nlon+1))
  allocate(rlat_edge(nlat+1))
  

  rlon_edge(1)=-180.0
  dlon=360.0/nlon
  do ix=2, nlon+1
     rlon_edge(ix)=rlon_edge(ix-1)+dlon
  end do
  
  
  dlat=180.0/nlat
  rlat_edge(1)=-90.0
  
  do iy=2, nlat+1
     rlat_edge(iy)=rlat_edge(iy-1)+dlat
  end do
  
  D_2_R=2.0*PI/360.0
  rlon_edge=D_2_R*rlon_edge
  rlat_edge=D_2_R*rlat_edge
  rlat_edge=sin(rlat_edge)
  
  do iy=1, nlat
     factor=Re*Re*(rlat_edge(iy+1)-rlat_edge(iy))
     do ix=1, nlon
        area(ix, iy)=factor*(rlon_edge(ix+1)-rlon_edge(ix))
     end do
  end do
  
end subroutine get_area_regular_grid

subroutine get_area_regular_edge(nlon, nlat, area)

  implicit none
#include "CMN_GCTM"  ! Physical constants
  
  integer, intent(in)::nlon, nlat
  real*8, intent(out)::area(nlon, nlat)
  real*8, allocatable::rlon_edge(:), rlat_edge(:)
  real*8 ::            D_2_R, factor
  integer::ix, iy
  real*8::dlon, dlat

  allocate(rlon_edge(nlon+1))
  allocate(rlat_edge(nlat))
  
  area=0.0
  
  rlon_edge(1)=-180.0
  dlon=360.0/nlon
  do ix=2, nlon+1
     rlon_edge(ix)=rlon_edge(ix-1)+dlon
  end do
  
  
  dlat=180.0/(nlat-1)
  rlat_edge(1)=-90.0
  
  do iy=2, nlat
     rlat_edge(iy)=rlat_edge(iy-1)+dlat
  end do

  
  D_2_R=2.0*PI/360.0
  rlon_edge=D_2_R*rlon_edge
  rlat_edge=D_2_R*rlat_edge
  rlat_edge=sin(rlat_edge)
  
  do iy=1, nlat-1
     factor=Re*Re*(rlat_edge(iy+1)-rlat_edge(iy))
     do ix=1, nlon
        area(ix, iy)=factor*(rlon_edge(ix+1)-rlon_edge(ix))
     end do
  end do
  area(:, nlat)=area(:, 1)
  
  deallocate(rlon_edge)
  deallocate(rlat_edge)
  
end subroutine get_area_regular_edge




subroutine do_flux_by_coef(nx, ny, nz, flux, nm, coef, flux_out)
implicit none 
! calculate matrix multiplication: flux*coef

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



subroutine regrid_flux_smp(lon1, lat1, nlon, nlat, &

     lon2, lat2, nlon2, nlat2, influx, outflux)

  !  regrid the flux from a smaller grid boxes to larger grid boxes 
  !  we need to check the map to make sure the surface type is consistent with each other 
  !  
  !  parameters 
  !  lon1 & lat1  in  array of nlon or nlat   centre of origin grid box 
  !  lon2 & lat2  in  array  of nlon2, nlat2, centre of target grid box
  !  inmap    in  array of (nlon, nlat), surface type  of origin grid box
  !  outmap   in  array of (nlon2, nlat2), surface type  of target grid box 
  !  influx    in  array of (nlon, nlat)  origin surface flux
  !  outflux    out array of (nlon2, nlat2) surface flux at target  grid box   
  
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
 
used=0
flux2=0.0

pst=1


! search for closest lon2 location for origin lon
call get_closest_point(nlon2, rlon2, nlon, rlon1, plon)

! search for closest lat2 location for origin lat

call get_closest_point(nlat2, rlat2,nlat, rlat1, plat)
! add the flux from small boxes to larger boxes 

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

