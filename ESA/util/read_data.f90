!***********************************************************************
! This program reads the TransCom 3 input variables out of the binary
! file, "input.dat".
!************************************************************************
!

subroutine read_transcom_map(flnm, ff90, ff95, nep, ocean, sf6, &
     landunit, oceanunit)
  
  character(len=*), intent(in)::flnm
  Real*8, intent(out) :: ff90(720,360),ff95(720,360)
  Real*8, intent(out) :: sf6(720,360,11),landunit(720,360,11)
  Real*8, intent(out) :: nep(720,360,12),ocean(720,360,12)
  Real*8, intent(out) :: oceanunit(720,360,11,12)
  
  print *, trim(flnm)
  
  Open(unit=10,file=trim(flnm),form= 'unformatted')
  Read(10) ff90
  Read(10) ff95
  Read(10) nep
  Read(10) ocean
  Read(10) landunit
  
  print *, 'sum', sum(ff90)
  
  Read(10) oceanunit
  Read(10) sf6
  
  Close(10)
  
end subroutine read_transcom_map


      
! ************************************************************************
