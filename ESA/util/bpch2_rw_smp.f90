SUBROUTINE READ_BPCH2_HEAD( FILENAME, TMPFLNM, STAT, TITLE, TRACER_ID, & 
     LONRES_OUT, LATRES_OUT, IX,  JX,  &
     LX, &
     IFIRST,  JFIRST,    LFIRST, &
     HALFPOLAR, CENTER180, &
     TAU0_OUT, TAU1_OUT, NTRACERS,MAXTRACER)
  
  !
  IMPLICIT NONE
  ! References to F90 modules
  
  ! Arguments OUT
  
  CHARACTER(LEN=*),  INTENT(IN)   :: FILENAME, TMPFLNM
  INTEGER*4,         INTENT(OUT)  :: IX(MAXTRACER), JX(MAXTRACER), LX(MAXTRACER), TRACER_ID(MAXTRACER)
  INTEGER*4,         INTENT(OUT) ::  IFIRST(MAXTRACER),    JFIRST(MAXTRACER),    LFIRST(MAXTRACER)
  INTEGER*4,         INTENT(OUT) ::  HALFPOLAR(MAXTRACER), CENTER180(MAXTRACER)
  
  INTEGER*4,         INTENT(OUT)  :: NTRACERS, STAT
  INTEGER*4,         INTENT(IN)   :: MAXTRACER
  REAL*4, INTENT(OUT)             :: TAU0_OUT(MAXTRACER), TAU1_OUT(MAXTRACER)
  REAL*4, INTENT(OUT)             :: LONRES_OUT(MAXTRACER), LATRES_OUT(MAXTRACER)
  
  CHARACTER(LEN=80),  INTENT(OUT) :: TITLE
  
  
  ! Local variables
  LOGICAL            :: FOUND, TMP_QUIET
  INTEGER            :: I,  J,  L,  N,  IOS, M
  INTEGER            :: I1, I2, J1, J2, L1,  L2
  CHARACTER(LEN=255) :: MSG
  CHARACTER(LEN=40)                        :: FTI
  CHARACTER(LEN=80)                        :: TMP_TITLE
  
  INTEGER            :: IU_FILE
  
  
  ! For binary punch file, version 2.0
  INTEGER            :: NTRACER,   NSKIP
  INTEGER            :: HALFPOLAR_R, CENTER180_R
  INTEGER            :: NI,        NJ,        NL, IUNIT
  INTEGER            :: IFIRST_R,    JFIRST_R,    LFIRST_R
  REAL*4             :: LONRES,    LATRES
  REAL*8             :: ZTAU0,     ZTAU1
  CHARACTER(LEN=20)  :: MODELNAME
  CHARACTER(LEN=40)  :: CATEGORY
  CHARACTER(LEN=40)  :: UNIT     
  CHARACTER(LEN=40)  :: RESERVED
  
  INTEGER::NESTED_CH, NESTED_NA

  REAL*4             :: TEMPARRAY(540,361,73)   ! for high resolution 

  
  !=================================================================
  ! READ_BPCH2 begins here!
  !  
  ! Initialize some variables
  !=================================================================
  
  
  IU_FILE=39
  STAT=0
  ! print *, TRIM(TMPFLNM)
  OPEN(UNIT=36, FILE=TRIM(TMPFLNM), STATUS='UNKNOWN')
  
  IUNIT=IU_FILE
  OPEN(IUNIT,      FILE=TRIM( FILENAME), STATUS='OLD', &
       IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL')
  
  READ( IUNIT, IOSTAT=IOS ) FTI
  
  IF (IOS.NE.0) THEN 
     STAT=IOS
     print *, STAT
     RETURN
  END IF
  
   
  READ( IUNIT, IOSTAT=IOS ) TMP_TITLE
  
  IF (IOS.NE.0) THEN 
     STAT=IOS
     print *, STAT
     RETURN
  END IF
  
  
  TITLE=TRIM(TMP_TITLE)
  ! Error check
  
  NTRACERS=0
  DO
     READ( IU_FILE, IOSTAT=IOS ) &
          MODELNAME, LONRES, LATRES, HALFPOLAR_R, CENTER180_R
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        EXIT 
     END IF
     
     READ( IU_FILE, IOSTAT=IOS )  &
          CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
          NI,       NJ,       NL,   IFIRST_R, JFIRST_R, LFIRST_R,   &
          NSKIP
     IF (IOS.NE.0) THEN 
        STAT=IOS
        EXIT 
     END IF
     
     NTRACERS=NTRACERS+1
     
     IX(NTRACERS)=NI
     JX(NTRACERS)=NJ
     LX(NTRACERS)=NL
     
     IFIRST(NTRACERS)=IFIRST_R
     JFIRST(NTRACERS)=JFIRST_R
     LFIRST(NTRACERS)=LFIRST_R
     
     HALFPOLAR(NTRACERS)=HALFPOLAR_R
     CENTER180(NTRACERS)=CENTER180_R
     
     
     LONRES_OUT(NTRACERS)=LONRES
     LATRES_OUT(NTRACERS)=LATRES
     TAU0_OUT(NTRACERS)=ZTAU0
     TAU1_OUT(NTRACERS)=ZTAU1
     TRACER_ID(NTRACERS)=NTRACER
     !  save to file  
     
     WRITE(36, *)  TRIM(CATEGORY),' ,',  TRIM(MODELNAME), &
          ' ,', TRIM(UNIT), ' , ', TRIM(RESERVED)
     ! read data to memory 
     
     READ( IU_FILE, IOSTAT=IOS ) &
          ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        EXIT 
     END IF
     

     IF (NTRACERS==MAXTRACER) EXIT
         
     
  ENDDO
  
  CLOSE(IU_FILE)
  CLOSE(36)
  
  
END SUBROUTINE READ_BPCH2_HEAD



!------------------------------------------------------------------------------

SUBROUTINE WRITE_BPCH2_HDR (IUNIT, TITLE, STAT)
  ! write bpch2 file record 
  
  INTEGER*4       ,    INTENT(IN) :: IUNIT
  CHARACTER(LEN=*),  INTENT(IN) :: TITLE
  INTEGER*4      ,   INTENT(OUT):: STAT    
  ! Local variable
  INTEGER                       :: IOS
  CHARACTER(LEN=40)             :: FTI = 'CTM bin 02'
  
  IOS=0
  
  WRITE ( IUNIT, IOSTAT=IOS ) FTI
  IF (IOS/=0) THEN
     STAT=IOS
     RETURN
  ENDIF

  WRITE ( IUNIT, IOSTAT=IOS ) TITLE
  
  IF (IOS/=0) THEN
     STAT=IOS
     RETURN
  ENDIF
  
END SUBROUTINE WRITE_BPCH2_HDR


SUBROUTINE OPEN_BPCH2_FOR_WRITE( IUNIT, FILENAME, TITLE, STAT )
  !
  !******************************************************************************
  !  Subroutine OPEN_BPCH2_FOR_WRITE opens a binary punch file (version 2.0)
  !  for writing. (bmy, 7/30/02)
  !   
  !  Arguments as Input:
  !  ============================================================================
  !  (1 ) IUNIT    (INTEGER )  : Logical unit number of the file to be opened
  !  (2 ) FILENAME (CHARACTER) : Name of the file to be opened
  !  (3 ) TITLE    (CHARACTER) : Optional: title for top of file
  !  
  !  NOTES:
  !******************************************************************************
  !  Feng changed it to be called from python
  ! 
  ! References to F90 modules
  
  ! Arguments
  INTEGER*4,           INTENT(IN)           :: IUNIT
  INTEGER*4,           INTENT(OUT)           :: STAT
  CHARACTER(LEN=*),  INTENT(IN)           :: FILENAME
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TITLE
  INTEGER::FUNIT
  ! Local variables
  INTEGER                                 :: IOS
  CHARACTER(LEN=80)                       :: TMP_TITLE
  
  !=================================================================
  ! OPEN_BPCH2_FOR_WRITE begins here!
  !=================================================================

  ! If TITLE is not passed, create a default title string
  IF ( PRESENT( TITLE ) ) THEN
     TMP_TITLE = TITLE
  ELSE
     TMP_TITLE = 'GEOS-CHEM binary punch file v. 2.0'
  ENDIF
  FUNIT=IUNIT
  
  ! Open file for output
  OPEN( FUNIT,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN', & 
       IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

  
  IF ( IOS /= 0 ) THEN
     STAT=IOS
     RETURN
  END IF
   
  CALL WRITE_BPCH2_HDR( IUNIT, TMP_TITLE, STAT )
  
  ! Return to calling program
END SUBROUTINE OPEN_BPCH2_FOR_WRITE


SUBROUTINE OPEN_BPCH2_FOR_READ(FUNIT, FILENAME, FTI_OUT, TITLE_OUT, STAT)

  !******************************************************************************
  !  Subroutine OPEN_BPCH2_FOR_READ opens a binary punch file (version 2.0)
  !  for reading (lfeng)
  !
  !   
  !  Inputs:
  !  ============================================================================
  !  (1 ) FUNIT    (INTEGER )  : Logical unit number of the file to be opened
  !  (2 ) FILENAME (CHARACTER) : Name of the file to be opened
  !   
  !  Outputs
  !  ==============================================
  !  (1)  FTI_OUT (CHARACTER)      :  file information 
  !  (2)  TITLE_OUT (CHARACTER)    :  title 
  !  (3)  IOS (INTEGER)            :  IO status 
  !
  
  IMPLICIT NONE
  
  INTEGER*4,          INTENT(IN)         :: FUNIT
  CHARACTER(LEN=*),  INTENT(IN)           :: FILENAME
  CHARACTER(LEN=40), INTENT(OUT)           :: FTI_OUT
  CHARACTER(LEN=80), INTENT(OUT)           :: TITLE_OUT
  INTEGER*4,           INTENT(OUT)        :: STAT
  
  ! local variable 
  CHARACTER(LEN=40)                        :: FTI
  CHARACTER(LEN=80)                        :: TITLE
  INTEGER                                  :: IOS  
  
  OPEN( FUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD', &
       IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL')
  FTI_OUT=""
  TITLE_OUT=""
  
  IF (IOS/=0) THEN
     STAT=IOS
     RETURN 
  END IF
  
  READ( FUNIT, IOSTAT=IOS ) FTI
  IF (IOS/=0) THEN
     STAT=IOS
     RETURN 
  END IF
  
  FTI_OUT=TRIM(FTI)
  
  READ( FUNIT, IOSTAT=IOS ) TITLE
  
  IF (IOS/=0) THEN
     STAT=IOS
     RETURN 
  END IF
  
  TITLE_OUT=TRIM(TITLE)
  RETURN 

END SUBROUTINE OPEN_BPCH2_FOR_READ


SUBROUTINE READ_BPCH2_RECORD(FUNIT, & 
     NTRACER_OUT, HALFPOLAR_OUT, CENTER180_OUT, &
     NI_OUT,        NJ_OUT,        NL_OUT, &
     IFIRST_OUT,    JFIRST_OUT,    LFIRST_OUT, &
     LONRES_OUT,    LATRES_OUT, &
     TAU0_OUT,     TAU1_OUT, &
     MODELNAME_OUT, CATEGORY_OUT, UNIT_OUT, RESERVED_OUT, & 
     OUT_ARRAY, &
     STAT)
  

  ! Read CURRENT records out of BPCH2 file 
  ! Inputs
  ! =========================================================
  !  (1) FUNIT (INTEGER) : file unit 
  
  ! Outputs
  ! ========================================================
  !  (1)  NTRACER_OUT(INTEGER)   : TRACER ID 
  !  (2)  HALFPOLAR_OUT(INTEGER) : POLAR SETTING 
  !  (3)  CENTER180_OUT(INTEGER) : CENTER OF LONGITUDE GRID 
  !  (4)  NI_OUT       (INTEGER) : LONGITUDE SIZE 
  !  (5)  NJ_OUT       (INTEGER) : LATITUDE  SIZE 
  !  (6)  NL_OUT       (INTEGER) : VERTICAL   SIZE 
  !  (7)  IFIRST_OUT   (INTEGER) : SHIFT IN LONGITUDE 
  !  (8)  JFIRST_OUT   (INTEGER) : SHIFT IN LATITUDE  
  !  (9)  LFIRST_OUT   (INTEGER) : SHIFT IN VERTICAL
  !  (10) lONRES_OUT   (FLOAT)   : LONGITUDE RESOLUTION 
  !  (11) LATRES_OUT   (FLOAT)   : LATITUDE RESOLUTION
  !  (12) TAU0_OUT     (FLOAT)   : STARTING TIME
  !  (13) TAU1_OUT     (FLOAT)   : END TIME
  !  (14) MODELNAME_OUT (CHARACTER): MODEL NNAME
  !  (15) CATEGORY_OUT (CHARACTER): TRACER CATEGORY
  !  (16) UNIT_OUT (CHARACTER)    : UNIT 
  !  (17) RESERVED_OUT (CHARACTER): RESERVED WORDS
  !  (18) OUT_ARRAY    (ARRAY)    : DATA
  !  (19) STAT         (INTEGER)  :  IO status 
  
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: FUNIT
 
  INTEGER*4, INTENT(OUT)    :: NTRACER_OUT
  INTEGER*4, INTENT(OUT)           :: HALFPOLAR_OUT, CENTER180_OUT
  INTEGER*4, INTENT(OUT)          :: NI_OUT,        NJ_OUT,        NL_OUT
  INTEGER*4, INTENT(OUT)           :: IFIRST_OUT,    JFIRST_OUT,    LFIRST_OUT
  REAL*4,     INTENT(OUT)        :: LONRES_OUT,    LATRES_OUT
  REAL*4,     INTENT(OUT)        :: TAU0_OUT,     TAU1_OUT
  CHARACTER(LEN=20), INTENT(OUT) :: MODELNAME_OUT
  CHARACTER(LEN=40),  INTENT(OUT):: CATEGORY_OUT
  CHARACTER(LEN=40),  INTENT(OUT) :: UNIT_OUT
  CHARACTER(LEN=40),  INTENT(OUT) :: RESERVED_OUT

  REAL*4             ::  TEMPARRAY(540,361,73)   
  
  
  REAL*4, INTENT(OUT)            :: OUT_ARRAY(540,361,73)

  
  INTEGER*4,         INTENT(OUT) :: STAT
  
  ! Local variables
  LOGICAL            :: FOUND, TMP_QUIET
  INTEGER            :: I,  J,  L,  N,  IOS, M
  INTEGER            :: I1, I2, J1, J2, L1,  L2
  
  ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
  INTEGER            :: IU_FILE
  
  
  ! For binary punch file, version 2.0
  INTEGER            :: NTRACER,   NSKIP
  INTEGER            :: HALFPOLAR, CENTER180
  INTEGER            :: NI,        NJ,        NL, IUNIT
  INTEGER            :: IFIRST,    JFIRST,    LFIRST
  REAL*4             :: LONRES,    LATRES
  REAL*8             :: ZTAU0,     ZTAU1
  CHARACTER(LEN=20)  :: MODELNAME
  CHARACTER(LEN=40)  :: CATEGORY
  CHARACTER(LEN=40)  :: CATEGORY_TMP
  CHARACTER(LEN=40)  :: UNIT     
  CHARACTER(LEN=40)  :: RESERVED
  
  INTEGER::NESTED_CH, NESTED_NA

  OUT_ARRAY(:,:,:) = 0e0
  
  
  ! Define a temporary variable for QUIET
  
  ! CONVERSION OF CHARACTER NEEDED WHEN CALLED BY PYTHON
  ! Error check
  
  STAT=0
  READ( FUNIT, IOSTAT=IOS ) &
       MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
  !     print *, MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
  IF (IOS.NE.0) THEN 
     STAT=IOS
     RETURN
  END IF
  MODELNAME_OUT=MODELNAME
  LONRES_OUT=LONRES
  LATRES_OUT=LATRES
  HALFPOLAR_OUT=HALFPOLAR
  CENTER180_OUT=CENTER180
     !    
  READ( FUNIT, IOSTAT=IOS )  &
       CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
       NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,   &
       NSKIP                                               
  
  IF (IOS.NE.0) THEN 
     STAT=IOS
     RETURN
  END IF
  
  CATEGORY_OUT=CATEGORY
  NTRACER_OUT=NTRACER
  UNIT_OUT=UNIT
  TAU0_OUT=ZTAU0
  TAU1_OUT=ZTAU1
  RESERVED_OUT=RESERVED
  NI_OUT=NI
  NJ_OUT=NJ
  NL_OUT=NL
  IFIRST_OUT=IFIRST
  JFIRST_OUT=JFIRST
  LFIRST_OUT=LFIRST
     
  
  READ( FUNIT, IOSTAT=IOS ) &
       ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
  
!  PRINT *, 'READ RECORD==>NTRACER, NI, NJ, NL, MAXVAL', NTRACER, NI, NJ, NL, MAXVAL(OUT_ARRAY(1:NI,1:NJ,1:NL))
  
  
  
  IF (IOS.NE.0) THEN 
     STAT=IOS
     RETURN
  ELSE

     I1 = IFIRST
     J1 = JFIRST
     L1 = LFIRST
     
     I2 = NI + I1 - 1
     J2 = NJ + J1 - 1
     L2 = NL + L1 - 1
     
     OUT_ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
     
  END IF
  
     
END SUBROUTINE READ_BPCH2_RECORD


SUBROUTINE WRITE_BPCH2_DATA( IUNIT,     MODELNAME_IN, CATEGORY_IN,  RESERVED_IN, & 
     LONRES,   LATRES,&
     HALFPOLAR, CENTER180,  NTRACER, & 
     UNIT_IN,      TAU0_IN,      TAU1_IN,   &  
     NI,        NJ,        NL,       IFIRST,     & 
     JFIRST,    LFIRST,    ARRAY_IN, STAT)
  !
  !******************************************************************************
  !  Subroutine BPCH2 writes binary punch file (version 2.0) to disk.
  !  Information about the model grid is also stored with each data block.
  !  (bmy, 5/27/99, 7/30/02)
  !   
  !  Arguments as input:
  !  ============================================================================
  !  (1    ) IUNIT      : INTEGER  - logical unit number of the file 
  !  (2    ) MODELNAME  : CHAR*40  - Name of model used to create output
  !  (3    ) LONRES     : REAL*4   - Longitude resolution of grid, in degrees
  !  (4    ) LATRES     : REAL*4   - Latitude resolution of grid, in degrees
  !  (4    ) HALFPOLAR  : INTEGER  - flag, =1 if model has half-polar boxes
  !  (5    ) CENTER180  : INTEGER  - flag, =1 if model has lon center on 180 deg
  !  (6    ) CATEGORY   : CHAR*40  - diagnostic category name
  !  (7    ) NTRACER    : INTEGER  - number of tracer
  !  (8    ) UNIT       : CHAR*40  - units of data
  !  (9    ) TAU0       : REAL*8   - TAU at start of diagnostic interval
  !  (10   ) TAU1       : REAL*8   - TAU at end   of diagnostic interval
  !  (11   ) RESERVED   : CHAR*40  - Reserved for future use
  !  (12-14) NI,NJ,NL   : INTEGER  - dimensions of ARRAY
  !  (15   ) IFIRST     : INTEGER  - I-index of the first grid box
  !  (16   ) JFIRST     : INTEGER  - J-index of the first grid box
  !  (17   ) LFIRST     : INTEGER  - L-index of the first grid box
  !  (18   ) ARRAY      : REAL*4   - data block to be written to the file
  !
  !  NOTES:
  !  (1 ) Added indices to IOERROR calls (e.g. "bpch2:1", "bpch2:2", etc.) 
  !        (bmy, 10/4/99)
  !  (2 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
  !  (3 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
  !  (4 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
  !  (5)  LF CHANGE IT TO PYTHON INTERFACE
  !******************************************************************************
  !  
  ! References to F90 modules
  ! Arguments
  INTEGER*4,           INTENT(IN) :: IUNIT
  INTEGER*4,           INTENT(IN) :: NTRACER 
  INTEGER*4,           INTENT(IN) :: NI, NJ, NL 
  INTEGER*4,           INTENT(IN) :: IFIRST, JFIRST, LFIRST
  INTEGER*4,           INTENT(IN) :: HALFPOLAR, CENTER180
  REAL*4,            INTENT(IN) :: ARRAY_IN( NI, NJ, NL )
  REAL*4,            INTENT(IN) :: LONRES, LATRES
  REAL*4,            INTENT(IN) :: TAU0_IN,   TAU1_IN
  INTEGER*4,         INTENT(OUT):: STAT
  REAL*4,            ALLOCATABLE :: ARRAY(:, :, : )
  REAL*8                         :: TAU0,   TAU1
  
  CHARACTER(LEN=20), INTENT(IN) :: MODELNAME_IN
  CHARACTER(LEN=40), INTENT(IN) :: CATEGORY_IN
  CHARACTER(LEN=40), INTENT(IN) :: RESERVED_IN
  CHARACTER(LEN=40), INTENT(IN) :: UNIT_IN
  
  CHARACTER(LEN=20)  :: MODELNAME
  CHARACTER(LEN=40)  :: CATEGORY
  CHARACTER(LEN=40) :: RESERVED
  CHARACTER(LEN=40) :: UNIT
  
      
  INTEGER                       :: I, J, L, NSKIP, IOS

  ! For computing NSKIP
  INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 4
  INTEGER, PARAMETER            :: END_OF_RECORD    = 8
  
  !  TYPE CONVERT NEEDED WHEN CALLED BY PYTHON  
  IF (ALLOCATED(ARRAY)) DEALLOCATE(ARRAY)
  ALLOCATE(ARRAY(NI,NJ, NL))
  ARRAY=ARRAY_IN
  READ(MODELNAME_IN, *) MODELNAME
  READ(CATEGORY_IN, *)  CATEGORY
  UNIT=TRIM(UNIT_IN)
  ! print *, trim(unit), trim(unit_in)
  ! READ(UNIT_IN, *)  UNIT
  READ(RESERVED_IN, *) RESERVED
  
  ! RESERVED=" "
  
  TAU0=TAU0_IN
  TAU1=TAU1_IN
  
  

  ! Local variables
  
  !=================================================================
  ! BPCH2 begins here!!  
  !
  ! Compute the number of bytes to skip between the end of one 
  ! data block and the beginning of the next data header line
  !=================================================================
  NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD
  
  !=================================================================
  ! Write data block to binary punch file
  ! Check for I/O errors
  !=================================================================
  WRITE( IUNIT, IOSTAT=IOS ) & 
       MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
  
  IF ( IOS /= 0 ) THEN
     STAT=IOS
     RETURN
  END IF
  
  WRITE( IUNIT, IOSTAT = IOS ) & 
       CATEGORY, NTRACER,  UNIT, TAU0,   TAU1,   RESERVED, & 
       NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, &  
       NSKIP
  
  IF ( IOS /= 0 ) THEN
     STAT=IOS
     RETURN
  END IF
  
  
  WRITE( IUNIT, IOSTAT=IOS ) &  
       ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
  
  ! print*,  size(ARRAY(:,1,1)), size(ARRAY(1,:,1)), size(ARRAY(1,1,:))
  ! print*,  NI, NJ, NL
  ! print*,  ARRAY(NI-5:NI,1,1)  
  IF ( IOS /= 0 ) THEN
     STAT=IOS
     RETURN
  END IF
  
  
  !=================================================================
  ! Return to calling program      
  !=================================================================
END SUBROUTINE WRITE_BPCH2_DATA


SUBROUTINE CLOSE_BPCH2_FILE( IUNIT, STAT )
  INTEGER*4, INTENT(IN)::IUNIT
  INTEGER*4, INTENT(OUT)::STAT
  
  ! LOCAL VARIABLE
  INTEGER IOSTAT
  CLOSE(IUNIT, IOSTAT=IOS)
  STAT=IOS
END SUBROUTINE CLOSE_BPCH2_FILE


