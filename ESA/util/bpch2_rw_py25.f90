SUBROUTINE READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, & 
     TAU0_IN,  IX,          JX,          &
     LX,       OUT_ARRAY, UNIT_OUT, STAT ) 
  
  IMPLICIT NONE
  ! References to F90 modules
  ! Arguments
#     include "define.h"
#     include "CMN_SIZE"   

  INTEGER*4,           INTENT(IN)  :: IX, JX, LX, TRACER_IN
  CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME, CATEGORY_IN 
  REAL*4,            INTENT(IN)  :: TAU0_IN
  REAL*4,            INTENT(OUT) :: OUT_ARRAY(IX, JX, LX)      
  CHARACTER(LEN=40), INTENT(OUT)                        :: UNIT_OUT
    
  INTEGER*4,         INTENT(OUT) :: STAT
  ! Local variables
  LOGICAL            :: FOUND, TMP_QUIET
  INTEGER            :: I,  J,  L,  N,  IOS, M
  INTEGER            :: I1, I2, J1, J2, L1,  L2
  CHARACTER(LEN=255) :: MSG
  CHARACTER(LEN=40)                        :: FTI
  CHARACTER(LEN=80)                        :: TMP_TITLE
  
  ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
  
  !REAL*4             :: TEMPARRAY(360,181,70)

#if   defined( GRID05x0666 ) 
  REAL*4             :: TEMPARRAY(540,361,150)   
#else
  REAL*4             :: TEMPARRAY(360,181,150)   
#endif
  
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

  !=================================================================
  ! READ_BPCH2 begins here!
  !  
  ! Initialize some variables
  !=================================================================
!  PRINT *, IX, JX, LX, TRACER_IN
!  ALLOCATE(OUT_ARRAY(IX,JX,LX))
  
  print *, 'tau0', TAU0_IN,size(OUT_ARRAY)
  FOUND            = .FALSE.
  
  OUT_ARRAY(:,:,:)     = 0e0
  TEMPARRAY(:,:,:) = 0e0
  
  
  ! Define a temporary variable for QUIET
  TMP_QUIET = .FALSE.
  IU_FILE=39
  
  IUNIT=IU_FILE
  ! CONVERSION OF CHARACTER NEEDED WHEN CALLED BY PYTHON
  READ(CATEGORY_IN, *) CATEGORY_TMP
  
  OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD', &
       IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL')
  READ( IUNIT, IOSTAT=IOS ) FTI
  READ( IUNIT, IOSTAT=IOS ) TMP_TITLE
  ! Error check
  
  STAT=-1
  DO
     READ( IU_FILE, IOSTAT=IOS ) &
          MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
  !   print *, MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
     
     READ( IU_FILE, IOSTAT=IOS )  &
          CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
          NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,   &
          NSKIP                                               
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
     
     
     READ( IU_FILE, IOSTAT=IOS ) &
          ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
         
     
     ! Test for a match
     !  IF ('DUSTSRCE' == TRIM( CATEGORY_TMP )) & 
     ! print *, 'search match', TRACER_IN, NTRACER, TRIM( CATEGORY_IN),' ',  TRIM( CATEGORY ), TAU0_IN, ZTAU0
     
     IF ( TRIM( CATEGORY_TMP ) == TRIM( CATEGORY ) .and.  TRACER_IN  == NTRACER .and. TAU0_IN==ZTAU0) THEN
        ! PRINT *, 'FOUND==============='
        FOUND = .TRUE.
        STAT=0
        UNIT_OUT=TRIM(UNIT)
        CLOSE( IUNIT )
        EXIT
     ENDIF
     
  ENDDO
  
  !=================================================================unit
  ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
  ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
  ! the data block.
  !=================================================================
  IF ( FOUND ) THEN 

#if   defined( GRID1x1 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
         IF ((IX.EQ.1).AND.(JX.EQ.1)) THEN
            I1=1
            J1=1
            
            
         END IF
         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1
         OUT_ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
         ! PRINT *, 'FOUND', TEMPARRAY( 1:NI, 1:NJ, 1:NL )
         ! PRINT *, I1, I2, J1, J2, L1,L2, IFIRST, JFIRST
         ! PRINT *, OUT_ARRAY( 1:NI, 1:NJ, 1:NL )
         
  END IF
 
  
END SUBROUTINE READ_BPCH2

SUBROUTINE READ_BPCH2_FIRSTMATCH( FILENAME, CATEGORY_IN, &
     TRACER_IN, & 
     TAU0_IN, &
     CATEGORY_OUT, &
     TRACER_OUT, & 
     UNIT_OUT, &
     TAU0_OUT, &
     IX,       JX,          &
     LX,       OUT_ARRAY,  STAT ) 

  IMPLICIT NONE
  ! References to F90 modules
#     include "define.h"
#     include "CMN_SIZE"   
  ! Arguments
  INTEGER*4,         INTENT(IN)  :: TRACER_IN
  CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME, CATEGORY_IN 
  REAL*4,            INTENT(IN)  :: TAU0_IN
  
  INTEGER*4,         INTENT(OUT)  :: TRACER_OUT
  CHARACTER(LEN=40), INTENT(OUT)  :: CATEGORY_OUT 
  REAL*4,            INTENT(OUT)  :: TAU0_OUT



  REAL*4,      INTENT(OUT)  :: OUT_ARRAY(540,361,150)   
  CHARACTER(LEN=40), INTENT(OUT) :: UNIT_OUT
  INTEGER*4,         INTENT(OUT) :: IX, JX, LX
    
  INTEGER*4,         INTENT(OUT) :: STAT

  
  ! Local variables
  LOGICAL            :: FOUND, TMP_QUIET
  INTEGER            :: I,  J,  L,  N,  IOS, M
  INTEGER            :: I1, I2, J1, J2, L1,  L2
  CHARACTER(LEN=255) :: MSG
  CHARACTER(LEN=40)                        :: FTI
  CHARACTER(LEN=80)                        :: TMP_TITLE
  
  ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)

#if   defined( GRID05x0666 ) 
  REAL*4             :: TEMPARRAY(540,361,150)   
#else
  REAL*4             :: TEMPARRAY(360,181,150)   
#endif
  
  ! REAL*4             :: TEMPARRAY(360,181,70)
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
  INTEGER::MATCH_COUNT
  INTEGER::NESTED_CH, NESTED_NA

  !=================================================================
  ! READ_BPCH2 begins here!
  !  
  ! Initialize some variables
  !=================================================================
!  PRINT *, IX, JX, LX, TRACER_IN
!  ALLOCATE(OUT_ARRAY(IX,JX,LX))
  
!  print *, 'tau0', TAU0_IN,size(OUT_ARRAY)
  FOUND            = .FALSE.
  
  OUT_ARRAY(:,:,:)     = 0e0
  ! TEMPARRAY(:,:,:) = 0e0
  
  
  ! Define a temporary variable for QUIET
  TMP_QUIET = .FALSE.
  IU_FILE=39
  
  IUNIT=IU_FILE
  ! CONVERSION OF CHARACTER NEEDED WHEN CALLED BY PYTHON
  READ(CATEGORY_IN, *) CATEGORY_TMP
  
  OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD', &
       IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL')
  READ( IUNIT, IOSTAT=IOS ) FTI
  READ( IUNIT, IOSTAT=IOS ) TMP_TITLE
  ! Error check
  
  print *,  'select in f90', CATEGORY_TMP, TRACER_IN, TAU0_IN
  
  STAT=-1
  DO
     MATCH_COUNT=0
     READ( IU_FILE, IOSTAT=IOS ) &
          MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
  !   print *, MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
     
     READ( IU_FILE, IOSTAT=IOS )  &
          CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
          NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,   &
          NSKIP                                               
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
     
     
     READ( IU_FILE, IOSTAT=IOS ) &
          ( ( (TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        CLOSE( IUNIT )
        EXIT 
     END IF
         
     
     ! Test for a match
     !  IF ('DUSTSRCE' == TRIM( CATEGORY_TMP )) & 
     ! print *, 'search match', TRACER_IN, NTRACER, TRIM( CATEGORY_IN),' ',  TRIM( CATEGORY ), TAU0_IN, ZTAU0
     
     IF ( TRIM( CATEGORY_TMP ) == TRIM( CATEGORY ) .OR. &
          TRIM( CATEGORY_TMP ) == "NONE") MATCH_COUNT=MATCH_COUNT+1
       
     IF ( TRACER_IN  == NTRACER .OR. TRACER_IN==-999) MATCH_COUNT=MATCH_COUNT+1
   
     
     IF (ZTAU0>=TAU0_IN) MATCH_COUNT=MATCH_COUNT+1
     
     IF (MATCH_COUNT==3) THEN
        ! PRINT *, 'FOUND==============='
        FOUND = .TRUE.
        CLOSE( IUNIT )
        EXIT
     ENDIF
     
  ENDDO
  
  !=================================================================unit
  ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
  ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
  ! the data block.
  !=================================================================
  IF ( FOUND ) THEN 

#if   defined( GRID1x1 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
         IF ((NI.EQ.1).AND.(NJ.EQ.1)) THEN
            I1=1
            J1=1
            
            
         END IF
         ! VALUE OF PARAMETER
         
         STAT=0
         UNIT_OUT=TRIM(UNIT)
        
         TAU0_OUT=ZTAU0
         TRACER_OUT  = NTRACER
         CATEGORY_OUT =TRIM(CATEGORY )
        
         ! SIZE 
         
         IX=NI
         JX=NJ
         LX=NL

         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         OUT_ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
         
         
         ! PRINT *, 'FOUND', TEMPARRAY( 10, 10, 1)
         ! PRINT *, I1, I2, J1, J2, L1,L2, IFIRST, JFIRST
         ! PRINT *, OUT_ARRAY( 1:NI, 1:NJ, 1:NL )
         
  END IF
 
  
END SUBROUTINE READ_BPCH2_FIRSTMATCH



SUBROUTINE READ_BPCH2_HEAD( FILENAME, TMPFLNM, STAT, TITLE, TRACER_ID, & 
     LONRES_OUT, LATRES_OUT, IX,  JX,  &
     LX, &
     IFIRST,  JFIRST,    LFIRST, &
     HALFPOLAR, CENTER180, &
     TAU0_OUT, TAU1_OUT, NTRACERS,MAXTRACER)
  
  !
  IMPLICIT NONE
  ! References to F90 modules
  
  ! Arguments
  
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
  
  ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
  
  ! REAL*4             :: TEMPARRAY(360,181,190)
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

#if   defined( GRID05x0666 ) 
  REAL*4             :: TEMPARRAY(540,361,150)   
#else
  REAL*4             :: TEMPARRAY(360,181,150)   
#endif
  

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
     
     WRITE(36, *)  TRIM(CATEGORY),' ,',  TRIM(MODELNAME), &
          ' ,', TRIM(UNIT), ' , ', TRIM(RESERVED)
     
     
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
  
  IMPLICIT NONE
  ! References to F90 modules
#     include "define.h"
#     include "CMN_SIZE"   
  ! Arguments
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

#if   defined( GRID05x0666 ) 
  REAL*4             ::  TEMPARRAY(540,361,150)   
#else
  REAL*4             ::  TEMPARRAY(360,181,150)   
#endif
  
  
  REAL*4, INTENT(OUT)             :: OUT_ARRAY(540,361,150)

  
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

#if   defined( GRID1x1 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
         IF ((IFIRST.EQ.1).AND.(JFIRST.EQ.1)) THEN
            I1=1
            J1=1
            
            
         END IF
         
         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1
         ! PRINT *, 'INDEX RANGE', I1, I2, J1, J2, L1, L2
         ! PRINT *, 'NI, NJ, NL',  NI, NJ, NL
         ! PRINT *, 'MAX-MIN',  MAXVAL(TEMPARRAY), MINVAL(TEMPARRAY)
         
         OUT_ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
     
  END IF
  
     
END SUBROUTINE READ_BPCH2_RECORD


SUBROUTINE SEL_BPCH2_RECORD(FUNIT, CATEGORY_IN, & 
     TRACER_IN, TAU_IN, &
     NTRACER_OUT, HALFPOLAR_OUT, CENTER180_OUT, &
     NI_OUT,        NJ_OUT,        NL_OUT, &
     IFIRST_OUT,    JFIRST_OUT,    LFIRST_OUT, &
     LONRES_OUT,    LATRES_OUT, &
     TAU0_OUT,     TAU1_OUT, &
     MODELNAME_OUT, CATEGORY_OUT, UNIT_OUT, RESERVED_OUT, & 
     OUT_ARRAY, &
     STAT)
  
  IMPLICIT NONE
  ! References to F90 modules
#     include "define.h"
#     include "CMN_SIZE"   
  ! Arguments
  INTEGER, INTENT(IN) :: FUNIT
  CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN 
  INTEGER, INTENT(IN)            :: TRACER_IN
  REAL*8, INTENT(IN)             :: TAU_IN
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
  REAL*4, INTENT(OUT)             :: OUT_ARRAY(540,361,150)
  
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

#if   defined( GRID05x0666 ) 
  REAL*4             ::  TEMPARRAY(540,361,150)   
#else
  REAL*4             ::  TEMPARRAY(360,181,150)   
#endif
  

!  REAL*4             :: TEMPARRAY(360,181,70)

  INTEGER::NESTED_CH, NESTED_NA
  ! LOGICAL:: FOUND
  OUT_ARRAY(:,:,:) = 0e0
  TEMPARRAY(:,:,:) = 0e0

  
  READ(CATEGORY_IN, *) CATEGORY_TMP
  
  ! Define a temporary variable for QUIET
  
  ! CONVERSION OF CHARACTER NEEDED WHEN CALLED BY PYTHON
  ! Error check
  
  STAT=0
  DO 
     FOUND=.False.
     READ( FUNIT, IOSTAT=IOS ) &
          MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
     !     print *, MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        RETURN
     END IF
     
     !    
     READ( FUNIT, IOSTAT=IOS )  &
       CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
       NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,   &
       NSKIP                                               
     
     IF (IOS.NE.0) THEN 
        STAT=IOS
        RETURN
  
     END IF
  
     READ( FUNIT, IOSTAT=IOS ) &
       ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
  
     IF (IOS.NE.0) THEN 
        STAT=IOS
        RETURN
     END IF
     
     FOUND=(TRIM( CATEGORY_TMP ) == TRIM( CATEGORY ))
     IF (TAU_IN>0) THEN 
        FOUND=FOUND.AND.(ZTAU1>TAU_IN).AND.(ZTAU0<=TAU_IN)
     END IF
     
     IF (TRACER_IN>0) THEN 
        FOUND=FOUND.AND.TRACER_IN==NTRACER
     END IF
!     PRINT *, TRIM( CATEGORY_TMP ), TRACER_IN, TAU_IN
!     PRINT *, TRIM( CATEGORY ), NTRACER, ZTAU0, ZTAU1
!     PRINT *, 'FOUND', FOUND
     
     IF (FOUND) EXIT
     
  END DO
  
  IF (FOUND) THEN 
     MODELNAME_OUT=MODELNAME
     LONRES_OUT=LONRES
     LATRES_OUT=LATRES
     HALFPOLAR_OUT=HALFPOLAR
     CENTER180_OUT=CENTER180
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

#if   defined( GRID1x1 )     
#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
         IF ((IFIRST.EQ.1).AND.(JFIRST.EQ.1)) THEN
            I1=1
            J1=1
            
            
         END IF
         
         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1
         PRINT *, 'INDEX RANGE', I1, I2, J1, J2, L1, L2
         PRINT *, 'NI, NJ, NL',  NI, NJ, NL
         PRINT *, 'MAX-MIN',  MAXVAL(TEMPARRAY), MINVAL(TEMPARRAY)
         
         OUT_ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
!     OUT_ARRAY(1:NI, 1:NJ, 1:NL)=TEMPARRAY(1:NI,1:NL,1:NJ)
     STAT=0
     RETURN 
  
  END IF

  
  
     
END SUBROUTINE SEL_BPCH2_RECORD

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


!------------------------------------------------------------------------------

    
!------------------------------------------------------------------------------

      FUNCTION GET_MODELNAME() RESULT( MODELNAME )
!
!******************************************************************************
!  Function GET_MODELNAME returns the proper value of MODELNAME for GEOS-1,
!  GEOS-STRAT, GEOS-2, or GEOS-3 data.  MODELNAME is written to the binary
!  punch file and is used by the GAMAP package. (bmy, 6/22/00, 5/24/05)
!
!  NOTES:
!  (1 ) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (2 ) Added modelname for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (3 ) Added "GEOS4_30L" for reduced GEOS-4 grid.  Also now use C-preprocessor
!        switch "GRID30LEV" instead of IF statements. (bmy, 11/3/03)
!  (4 ) Updated for GCAP and GEOS-5 met fields.  Rearranged coding for
!        simplicity. (swu, bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"
#     include "CMN_SIZE"

      ! MODELNAME holds the return value for the function
      CHARACTER(LEN=20)   :: MODELNAME

      !=================================================================
      ! GET_MODELNAME begins here!
      !=================================================================

#if   defined( GEOS_1 ) 
      MODELNAME = 'GEOS1'
     
#elif defined( GEOS_STRAT ) 
      MODELNAME = 'GEOS_STRAT'

#elif defined( GEOS_3 ) && defined( GRID30LEV )
      MODELNAME = 'GEOS3_30L'

#elif defined( GEOS_3 )
      MODELNAME = 'GEOS3'

#elif defined( GEOS_4 ) && defined( GRID30LEV )
      MODELNAME = 'GEOS4_30L'

#elif defined( GEOS_4 )
      MODELNAME = 'GEOS4'

#elif defined( GEOS_5 ) && defined( GRID30LEV )
      MODELNAME = 'GEOS5_30L'
      
#elif defined( GEOS_5 ) 
      MODELNAME = 'GEOS5'

#elif defined( GCAP )
      MODELNAME = 'GCAP'

#endif

      ! Return to calling program
      END FUNCTION GET_MODELNAME

!------------------------------------------------------------------------------

      FUNCTION GET_NAME_EXT() RESULT( NAME_EXT )
!
!******************************************************************************
!  Function GET_NAME_EXT returns the proper filename extension for CTM
!  model name (i.e. "geos1", "geoss", "geos3", "geos4", "geos5", or "gcap").  
!  (bmy, 6/28/00, 5/24/05)
!  
!  NOTES:
!  (1 ) Added name string for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (2 ) Remove obsolete "geos2" model name strning (bmy, 11/3/03)
!  (3 ) Modified for GCAP and GEOS-5 met fields (bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"

#if   defined( GEOS_1 ) 
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos1'
     
#elif defined( GEOS_STRAT ) 
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geoss'

#elif defined( GEOS_3 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos3'

#elif defined( GEOS_4 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos4'

#elif defined( GEOS_5 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos5'

#elif defined( GCAP )
      CHARACTER(LEN=4) :: NAME_EXT
      NAME_EXT = 'gcap'

#endif
      
      ! Return to calling program
      END FUNCTION GET_NAME_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_NAME_EXT_2D() RESULT( NAME_EXT_2D )
!
!******************************************************************************
!  Function GET_NAME_EXT_2D returns the proper filename extension for CTM
!  model name for files which do not contain any vertical information
!  (i.e. "geos" or "gcap").  (bmy, 8/16/05)
!
!  NOTES: 
!******************************************************************************
!
      ! Local variables
      CHARACTER(LEN=4) :: NAME_EXT_2D
      CHARACTER(LEN=5) :: TEMP_NAME
      CHARACTER(LEN=5), EXTERNAL :: GET_NAME_EXT
      !=================================================================
      ! GET_NAME_EXT_2D begins here!
      !=================================================================

      ! Get the name extension
      TEMP_NAME   = GET_NAME_EXT()

      ! Take the 1st 4 characters ("geos" or "gcap") and return
      NAME_EXT_2D = TEMP_NAME(1:4)

      ! Return to calling program
      END FUNCTION GET_NAME_EXT_2D 

!------------------------------------------------------------------------------

      FUNCTION GET_RES_EXT() RESULT( RES_EXT )
!
!******************************************************************************
!  Function GET_RES_EXT returns the proper filename extension for
!  CTM grid resolution (i.e. "1x1", "2x25", "4x5").  (bmy, 6/28/00, 12/1/04)
! 
!  NOTES:
!  (1 ) Added extension for 1 x 1.25 grid (bmy, 12/1/04)
!******************************************************************************
!
#     include "define.h"
        CHARACTER(LEN=5):: RES_EXT

#if   defined( GRID4x5 )
!      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '4x5'
     
#elif defined( GRID2x25 ) 
!      CHARACTER(LEN=4) :: RES_EXT
      RES_EXT = '2x25'

#elif defined( GRID1x125 )
!      CHARACTER(LEN=5) :: RES_EXT
      RES_EXT = '1x125'

#elif defined( GRID1x1 ) 
!      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '1x1'

#endif
!       print *, RES_EXT

      END FUNCTION GET_RES_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_HALFPOLAR() RESULT( HALFPOLAR )
!
!******************************************************************************
!  Function GET_HALFPOLAR returns 1 if the current grid has half-sized polar
!  boxes (e.g. GEOS), or zero otherwise (e.g. GCAP).  (swu, bmy, 6/28/05)
!
!  NOTES: 
!******************************************************************************
!
#     include "define.h"

      ! Local variables
      INTEGER :: HALFPOLAR

      !=================================================================
      ! GET_HALFPOLAR begins here!
      !=================================================================

#if   defined( GCAP ) 

      ! GCAP grid does not have half-sized polar boxes
      HALFPOLAR = 0

#else

      ! All GEOS grids have half-sized polar boxes
      HALFPOLAR = 1

#endif

      ! Return to calling program
      END FUNCTION GET_HALFPOLAR

      SUBROUTINE READ_BPCH2_CO2_FLUX( FILENAME, STAT,&
           FLUX_ST,       TRACER_ST, &
           TRACER_END,  &
           IX, JX, LX, &
           ARRAY) 
        
      IMPLICIT NONE
#     include "define.h" 

        ! Arguments
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME 
      INTEGER, INTENT(OUT)           :: FLUX_ST, TRACER_ST, TRACER_END
      INTEGER, INTENT(OUT)           :: STAT
      REAL*8,            INTENT(OUT) :: ARRAY(360,181,500)      
      INTEGER, INTENT(OUT)            :: IX, JX, LX
      
      ! Local variables
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG
      
      ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
      REAL*4             :: TEMPARRAY(360,181,500)
      

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED
      INTEGER            :: IU_FILE  
      CHARACTER(LEN=40)  :: FTI
      CHARACTER(LEN=80)  :: TITLE
      INTEGER*4          :: STAT_OPEN
      
      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      ARRAY(:,:,:)     = 0e0
      TEMPARRAY(:,:,:) = 0e0
      IU_FILE=18
      
      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      !  PRINT *, TRIM(FILENAME)
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME, FTI, TITLE, STAT_OPEN )
      ! PRINT *, TRIM(FILENAME)
      
      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      READ( IU_FILE, IOSTAT=IOS )  &
           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
         
      
      IF ( IOS /= 0 ) THEN 
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN 
      END IF
      
      READ( IU_FILE, IOSTAT=IOS ) &
           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, &
           NSKIP
         
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      READ( IU_FILE, IOSTAT=IOS ) &
           ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
      
      
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      TRACER_ST=NTRACER
      FLUX_ST=2
      IF (TRACER_ST==1) FLUX_ST=1
      TRACER_END=TRACER_ST+NL-FLUX_ST
      ARRAY( 1:NI, 1:NJ, 1:NL ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )
      IX=NI
      JX=NJ
      LX=NL
      
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
    END SUBROUTINE READ_BPCH2_CO2_FLUX

    SUBROUTINE READ_BPCH2_CO2_STV( FILENAME, STAT,&
         FLUX_ST,       TRACER_ST, &
         TRACER_END,  &
         IX, JX, LX, &
         ARRAY) 
        
      IMPLICIT NONE
#     include "define.h" 

        ! Arguments
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME 
      INTEGER, INTENT(OUT)           :: FLUX_ST, TRACER_ST, TRACER_END
      INTEGER, INTENT(OUT)           :: STAT
      REAL*8,            INTENT(OUT) :: ARRAY(2000,500,2)      
      INTEGER, INTENT(OUT)            :: IX, JX, LX
      
      ! Local variables
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG
      
      ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
      REAL*4             :: TEMPARRAY(360,181,500)
      REAL*4             :: TEMPARRAY_2(2000,500,2)
      

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED
      INTEGER            :: IU_FILE  
      CHARACTER(LEN=40)  :: FTI
      CHARACTER(LEN=80)  :: TITLE
      INTEGER*4          :: STAT_OPEN
      
      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      ARRAY(:,:,:)     = 0e0
      TEMPARRAY(:,:,:) = 0e0
      TEMPARRAY_2(:,:,:) = 0e0
      
      IU_FILE=18
      
      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      !  PRINT *, TRIM(FILENAME)
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME, FTI, TITLE, STAT_OPEN )
      ! PRINT *, TRIM(FILENAME)
      
      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      READ( IU_FILE, IOSTAT=IOS )  &
           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
         
      
      IF ( IOS /= 0 ) THEN 
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN 
      END IF
      
      READ( IU_FILE, IOSTAT=IOS ) &
           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, &
           NSKIP
         
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      READ( IU_FILE, IOSTAT=IOS ) &
           ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
      
      
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      ! now into the section for state vector
      
      READ( IU_FILE, IOSTAT=IOS )  &
           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
      
      
      IF ( IOS /= 0 ) THEN 
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN 
      END IF
      
      
      READ( IU_FILE, IOSTAT=IOS ) &
           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, &
           NSKIP
         
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      READ( IU_FILE, IOSTAT=IOS ) &
           ( ( ( TEMPARRAY_2(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
      
      
      IF ( IOS /= 0 ) THEN
         STAT=IOS
         PRINT *, 'ERROR IN REAING', FILENAME
         RETURN
      END IF
      
      
      TRACER_ST=NTRACER
      FLUX_ST=2
      IF (TRACER_ST==1) FLUX_ST=1
      ! the state vector is stored as NST_Variable X NST_ENSEMBLE X 1 
      TRACER_END=TRACER_ST+NJ-FLUX_ST 
      ! PRINT*, NI, NJ, NL
      ! PRINT*, SIZE(ARRAY(:,1,1)), SIZE(ARRAY(1,:,1))
      ! PRINT*, SIZE(TEMPARRAY_2(:,1,1)), SIZE(TEMPARRAY_2(1,:,1))
      
      ARRAY( 1:NI, 1:NJ, 1:NL ) = TEMPARRAY_2( 1:NI, 1:NJ, 1:NL )
      IX=NI
      JX=NJ
      LX=NL
      
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
    END SUBROUTINE READ_BPCH2_CO2_STV


