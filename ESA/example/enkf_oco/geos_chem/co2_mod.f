! $Id: co2_mod.f,v 1.4 2006/11/07 19:01:56 bmy Exp $
      MODULE CO2_MOD
!
!******************************************************************************
!  Module CO2_MOD contains variables and routines used for the CO2 simulation.
!  (pns, bmy, 8/16/05, 9/27/06) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) EMFOSSCO2    : Array for fossil fuel CO2 emissions (annual mean)
!  (2 ) EMBIOCO2     : Balanced biosphere CO2 (CASA) emissions (daily)   
!  (3 ) EMOCCO2      : Ocean CO2 emissions (annual mean)  
!  (4 ) EMBIOBRNCO2  : Biomass burning emissions  
!  (5 ) EMBIOFUELCO2 : Biofuel emissions 
!  (6 ) EMBIONETCO2  : Net terrestrial CO2 emissions
!  (7 ) FMOL_CO2     : Molecular weight of CO2 
!  (8 ) LFOSSCO2     : Flag for switching on/off FOSSIL FUEL CO2 emissions
!  (9 ) LBIOCO2      : Flag for switching on/off BALANCED BIOSPHERE CO2 emiss. 
!  (10) LOCCO2       : Flag for switching on/off OCEAN CO2 emissions     
!  (11) LBIOBRNCO2   : Flag for switching on/off BIOMASS BURNING CO2 emissions 
!  (12) LBIOFUELCO2  : Flag for switching on/off BIOFUEL CO2 emissions 
!  (13) LBIONETCO2   : Flag for switching on/off NET BIOSPHERE EXCHANGE of CO2
!  (14) LUSECASANEP  : Flag for reading daily CASA NEP w/ diurnal cycle
!  (15) LPOST        : FLAG FOR READING POSTERIOR FLUXES

!  (16) XNUMOL_CO2   : molec CO2 / kg CO2 
!  (17) DO_BK_RUN    : Mode choice for tagged CO2 simulation 
!                        3 -- simulation with 1 OR MORE  tracers: the first one is the total, and the second one 
!                           is background. 
!                        4 -- ensemble run: 
!                           the first one is the total (prior fluxes) 
!                           the others are from perturbations. 
!                        5 -- run ONE OR MORE tracers: the first one is FORCED BY POSTERIOR FLUXES 
!  (17) EMS_PATH    : DIRECTORY FOR EMISSION DATA


!
!  Module Procedures:
!  ============================================================================
!  (1 ) EMISSCO2               : Emits CO2 into individual tracers
!  (2 ) READ_ANNUAL_FOSSILCO2  : Reads annual mean emission fields for CO2
!  (3 ) READ_ANNUAL_OCEANCO2   : Reads annual mean CO2 ocean emissions
!  (4 ) READ_BBIO_DAILYAVERAGE : Reads daily mean CASA Bal Bio CO2 (no diurnal)
!  (5 ) READ_BBIO_DIURNALCYCLE : Reads CASA NEP fluxes w/ imposed diurnal cycle
!  (6 ) READ_MONTH_BIOBRN_CO2  : Read monthly biomass burning emissions
!  (7 ) READ_ANNUAL_BIOFUELCO2 : Read annual mean biofuel emissions
!  (8 ) READ_ANNUAL_BIONET_CO2 : Read annual net terrestrial exchange
!  (9 ) INIT_CO2               : Allocates and initializes module arrays
!  (10) CLEANUP_CO2            : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "co2_mod.f"
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (3 ) diag04_mod.f           : Module w/ routines for CO2 diagnostics
!  (4 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (6 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (7 ) grid_mod.f             : Module w/ horizontal grid information
!  (8 ) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (9 ) time_mod.f             : Module w/ routines for computing time & date
!  (10) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) transfer_mod.f         : Module w/ routines to cast & resize arrays 
!
!  CO2 tracers:
!  ============================================================================
!  (1 ) Total CO2
!  (2 ) CO2 from oceans
!  (3 ) CO2 from fossil fuel 
!  (4 ) CO2 from balanced biosphere
!  (5 ) CO2 from biomass burning
!  (6 ) CO2 from biofuels
!
!  References:
!  ============================================================================
!  (1 ) Andres, R.J, G. Marland, I. Fung, and E. Matthews, "A 1x1 distribution
!        of carbon dioxide emissions from fossil fuel consumption and
!        cement manufacture", Glob. Biogeochem. Cycles, Vol 10, 419-429, 1996.
!  (2 ) Randerson, J.T, M.V. Thompson, T.J.Conway, I.Y. Fung, and C.B. Field,
!        "The contribution of terrestrial sources and sinks to trends in the
!        seasonal cycle of atmospheric carbon dioxide", Glob. Biogeochem. 
!        Cycles, Vol 11, 535-560, 1997.
!  (3 ) Takahashi, T, R. Feely, R. Weiss, R. Wanninkof, D. Chipman, 
!        S. Sutherland, and T. Takahashi, "Global air-sea flux of CO2: An
!        estimate based on measurements of sea-air pCO2 difference", 
!        Proceedings of the National Academy of Sciences, 94, 8292-8299,
!        1997
!  (4 ) Yevich, R. and J. A. Logan, "An assesment of biofuel use and burning 
!        of agricultural waste in the developing world", Glob. Biogeochem. 
!        Cycles, Vol 17, 1095, doi:10.1029/2002GB001952, 2003
!
!  NOTES: 
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now references biomass_mod.f (bmy, 9/27/06)
!******************************************************************************
!
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "co2_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 

      ! ... except these routines
      PUBLIC :: CLEANUP_CO2
      PUBLIC :: EMISSCO2, INIT_CO2_ENKF
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Logical switches
      ! lf for test
      LOGICAL, PARAMETER   :: LFOSSCO2     = .TRUE.
      LOGICAL, PARAMETER   :: LBIOCO2      = .TRUE.
      LOGICAL, PARAMETER   :: LOCCO2       = .TRUE.
      LOGICAL, PARAMETER   :: LBIOBRNCO2   = .TRUE.
      LOGICAL, PARAMETER   :: LBIOFUELCO2  = .FALSE.
      LOGICAL, PARAMETER   :: LUSECASANEP  = .TRUE.
      LOGICAL, PARAMETER   :: LPOST     = .TRUE.
      
!     LF 2009.10.6
!     INCLUDE FOSSIL AND BIO FUEL OR NOT IN THE TAGGED REGIONS      
      LOGICAL, PARAMETER   :: TAG_FUEL    = .FALSE.
      
      
! lf test end

      LOGICAL, PARAMETER   :: LBIONETCO2   = .FALSE. 
      
! LF   THE USE OF INVERSION RESULTS
      LOGICAL, PARAMETER   ::  DO_FF_SMOOTH =.FALSE.
      
      
      ! Arrays
      REAL*8,  ALLOCATABLE :: EMFOSSCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMFOSSCO2_NEXT(:,:)
      
      REAL*8,  ALLOCATABLE :: EMOCCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOBRNCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOFUELCO2(:,:)
      
      REAL*8,  ALLOCATABLE :: EMBIONETCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOMASSCO2(:,:)
     
      
      
!     scaling applied to  biosphere emissions

      REAL*8,  PARAMETER   :: PB_BIO_FACTOR=1.0

!     scaling applied to  fossil fuel  emissions

      REAL*8,  PARAMETER   :: PB_FF_FACTOR=1.0

!     scalin apply to  all emissions
      
      REAL*8,  PARAMETER   :: FLUX_FACTOR=1.0
      
!     location of emissions
      CHARACTER(LEN=180), PARAMETER::EMS_PATH='./surface_flux/'
      

     
!     added by lf for  ENKF simulation
      
      REAL*8,  ALLOCATABLE :: POST_FLUX(:,:,:)   

!     Basis functions for flux perturbations 
      
      REAL*8,  ALLOCATABLE :: BF_FLUX(:,:,:)   
      INTEGER              :: NBF ! number of basis functions, (read in from bf flux files)
      
      
      ! 
      ! THE STARTING DATE OF CO2 EMISSION PERTURBATIONS, 
      ! USED WHEN THE EMISSIONS ARE READ FROM FILES FROM EMISSION_PATH
      ! 
      INTEGER, ALLOCATABLE :: ENR_YYYY(:), ENR_DOY(:), ENR_PBUSE(:)
      
      REAL*8,  ALLOCATABLE  :: ENR_DD(:)
      ! CURRENT TIME 
      INTEGER              :: CUR_YYYY=0
      INTEGER              :: CUR_DOY=0
      INTEGER              :: CUR_PBUSE=0
      INTEGER, PARAMETER   :: MAX_PBFLUX=150
      
      !                      
      INTEGER              :: ENR_NTIME
     
      ! ENSEMBLE FIRST AND LAST MEMBER 
      
      
      INTEGER              :: ENR_ST
      INTEGER              :: ENR_END
      !                       THE POSITION OF ENR_ST AT FLUX BF FILE (i.e ensemble files)
      
      INTEGER              :: ENR_BFST
      
      ! THE NAME OF THE OUTSIDE EMISSION FILES
      CHARACTER(LEN=255)   :: ENR_FILE
      
      ! THE FIRST AND THE LAST MEMBER OF OUSIDE EMISSION FILES 
      INTEGER              :: TRACER_ST, TRACER_END, FLUX_ST, DO_BK_RUN
      
      ! RESERVED NUMBER OF TRACERS
      
      INTEGER               :: NTR_RESVD
      !end ENKF extension by lf
      
                                ! 
      ! FMOL_CO2     - kg CO2 / mole CO2 
      REAL*8,  PARAMETER   :: FMOL_CO2   = 44d-3

      ! XNUMOL_CO2   - molecules CO2 / kg CO2 
      REAL*8,  PARAMETER   :: XNUMOL_CO2 = 6.022d+23 / FMOL_CO2
      REAL*8,  PARAMETER   :: EMFACTCO2CO = 12.068

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
      

!-----------------------------------------------------------------------------
   
      SUBROUTINE EMISSCO2
!
!******************************************************************************
!  Subroutine EMISSCO2 is the driver routine for CO2 emissions. 
!  (pns, bmy, 8/16/05, 9/27/06)
!
!  The initial condition for CO2 has to be at least 50 ppm or higher or else
!  the balanced biosphere fluxes will make STT negative. (pns, bmy, 8/16/05)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) We now get CO2 biomass emissions from biomass_mod.f.  This allows us 
!        to use either GFED2 or default Duncan et al biomass emissions. 
!        (bmy, 9/27/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS,       IDBCO2
      USE DIAG04_MOD,    ONLY : AD04,          ND04
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_DAY,       GET_DAY_OF_YEAR
      USE TIME_MOD,      ONLY : GET_HOUR,      GET_MONTH
      USE TIME_MOD,      ONLY : GET_YEAR,      GET_TS_CHEM 
      USE TIME_MOD,      ONLY : ITS_A_NEW_DAY, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT
      
#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      LOGICAL, SAVE          :: FIRST  = .TRUE.

      ! Local variables
      INTEGER                :: I,     IJLOOP, J,    L,     N
      INTEGER                :: DAY,   DOY,    HOUR, MONTH, YEAR   
      REAL*8                 :: A_CM2, DTSRCE, E_CO2
      
      REAL*8                 :: RDD
      INTEGER                :: TRACER_ID
      
      !=================================================================
      ! EMISSCO2 begins here!
      !=================================================================
      
      print *, '================== in co2 mod -------------' 
      
      ! First-time initialization
      IF ( FIRST ) THEN 
         print *, 'do co2 init' 
         
         ! Allocate arrays and read annual-mean data
         CALL INIT_CO2
         
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read in monthly and daily emissions fields
      !=================================================================      

      ! Emission timestep (assume it is called at chemistry step as well)
      
      DTSRCE = 60d0 * GET_TS_CHEM()
      
      ! Time variables
      DAY    = GET_DAY()
      DOY    = GET_DAY_OF_YEAR()
      HOUR   = GET_HOUR()
      MONTH  = GET_MONTH()
      YEAR   = GET_YEAR()
      ! rdd is used to read emission bf 
      
      RDD=(YEAR-1985.)*365.24+DOY
      
      
      !Check if Balanced Biosphere emissions are required  
      
      IF ( LBIOCO2 ) THEN  
         IF ( LUSECASANEP ) THEN
            
            IF ( MOD( HOUR, 3 ) == 0 ) THEN
               print *, 'CALL CASA_BIO'
               CALL READ_BBIO_DIURNALCYCLE( YEAR, MONTH, DAY, HOUR )
               ! CALL READ_CASA_BIO( YEAR, MONTH, DAY, HOUR )
               
            ENDIF
            
         ELSE
            
            ! ... otherwise use constant daily emissions of NEP for Bal Bio
            
            IF ( ITS_A_NEW_DAY() ) THEN
               
              CALL READ_BBIO_DAILYAVERAGE( YEAR, MONTH, DAY, DOY ) 
            
           ENDIF

         ENDIF
         
         
         
      ENDIF
      
      IF ( LOCCO2 ) THEN
         IF ( ITS_A_NEW_DAY() ) THEN
            CALL READ_MONTHLY_OCEANCO2( YEAR, MONTH, DAY) 
         END IF   
        
      END IF
      
      
      IF (LBIOBRNCO2) THEN
         IF ( ITS_A_NEW_DAY() ) THEN
            
            IF (LUSECASANEP) THEN
               CALL READ_DAILY_BIOMASSCO2( YEAR, MONTH, DAY)
            ELSE
               CALL READ_MONTHLY_BIOMASSCO2( YEAR, MONTH, DAY)
            END IF 
         ENDIF
      END IF
      
      
      
      IF ( MOD( HOUR, 3 ) == 0 ) THEN
         IF (DO_BK_RUN==5) THEN
            PRINT *, 'CALL READ POST FLUX FROM DO_BK_RUN==5'
            CALL READ_POST_FLUX(YEAR, MONTH, DAY)
         ELSE IF (DO_BK_RUN==4) THEN 
            PRINT *, 'CALL READ PERTURB FLUX FROM DO_BK_RUN==4'
            CALL READ_BASIS_FUNCTION(YEAR, MONTH, DAY, HOUR)
         END IF
      END IF
      
      
      
      !=================================================================
      ! use emissions to force model simulations  
      !=================================================================

      
      
      IF (DO_BK_RUN==4) THEN 
      
         print *, 'CALL RUN_ENKF (forced by flux perturbations)' 
         CALL RUN_ENKF()
         print *, 'AFTER CALL RUN_ENKF' 
         
      ELSE IF (DO_BK_RUN==5) THEN

!     run with posterior regional flux map      
         PRINT *,'CALL RUN_BK_5'
         
         print *, 'The first one tracer is forced by prior emissions' 
         print *,  'and posterior flux corrections' 
         
         CALL RUN_BK_5()

         print *, 'AFTER CALL RUN_ENKF_5' 
         
      ELSE
!     run with 1 or more tracers. Only the first one is forced by emissions

         PRINT *,'CALL RUN_BK_3' 
         print *, 'Only first one tracer is forced by emissions' 
         CALL RUN_BK_3()
         print *, 'AFTER RUN_BK_3' 
      END IF
      
      ! Return to calling program
      END SUBROUTINE EMISSCO2
      

      SUBROUTINE RUN_BK_3()
      
!     References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS,       IDBCO2
      USE DIAG04_MOD,    ONLY : AD04,          ND04
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT
      USE TIME_MOD,      ONLY : GET_TS_CHEM,  GET_DAY_OF_YEAR,
     &     GET_MONTH
#     include "CMN_SIZE"        ! Size parameters

      INTEGER                :: I,  DOY, J,    L,     N
      INTEGER                :: MONTH
      REAL*8                 :: A_CM2, DTSRCE, E_CO2
      
      
      DTSRCE = 60d0 * GET_TS_CHEM()
      DOY    = GET_DAY_OF_YEAR()
      MONTH   = GET_MONTH()
      
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, A_CM2, E_CO2 )
      
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

!-------------------------------------------
! #2: CO2 from fossil fuel emissions 
!-------------------------------------------
            IF ( LFOSSCO2 ) THEN
! Fossil fuel emissions of CO2 [molec/cm2/s]

               E_CO2          = EMFOSSCO2(I,J)
               IF (DO_FF_SMOOTH) E_CO2  = E_CO2+ 
     &              EMFOSSCO2_NEXT(I,J)*(MONTH-1)
               
! ND04 diag: Fossil Fuel CO2 [molec/cm2/s] 
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,1) = AD04(I,J,1) + FLUX_FACTOR*E_CO2
               ENDIF

!     Convert from [molec/cm2/s] to [kg]
               
               E_CO2  = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
!     Add to Tracer #1: Total CO2 [kg]
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
               
               
               
            ENDIF
         
!-------------------------------------------
! #3: CO2 from ocean emissions
!-------------------------------------------
         
            IF ( LOCCO2 ) THEN
               
! Ocean CO2 emissions in [molec/cm2/s]
         
               E_CO2          = EMOCCO2(I,J)
               
! ND04 diag: Ocean CO2 [molec/cm2/s]
               
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,2) = AD04(I,J,2) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert from [molec/cm2/s] to [kg]
               
               E_CO2  = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
!     Add to Tracer #1: Total CO2 [kg]
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
               
               
            ENDIF
         
!-------------------------------------------
!     #4: CO2 from balanced biosphere emissions
!-------------------------------------------
         
            IF ( LBIOCO2 ) THEN
!     Balanced biosphere CO2 [molec/cm2/s]
               E_CO2         = EMBIOCO2(I,J)
               
! ND04 diag: Bal Bio CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,3) = AD04(I,J,3) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert from [molec/cm2/s] to [kg CO2]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
               
            ENDIF

!-------------------------------------------
! #5: CO2 from biomass burning emissions
!-------------------------------------------
            IF ( LBIOBRNCO2 ) THEN
               
               E_CO2   = EMBIOMASSCO2(I,J)
            
               
! ND04 diag: Biomass burning CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,4) = AD04(I,J,4) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert from [molec/cm2/s] to [kg]
               E_CO2  = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
            
            ENDIF
         
!-------------------------------------------
!     #6: CO2 from biofuel emissions
!-------------------------------------------
            IF ( LBIOFUELCO2 ) THEN
               
! Biofuel CO2 emissions [molec/cm2/s]
               E_CO2          = EMBIOFUELCO2(I,J)
! ND04 diag: terrial CO2 [molec/cm2/s] 
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,5) = AD04(I,J,5) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
            
            ENDIF
            
         
         ENDDO
      ENDDO
      
!$OMP END PARALLEL DO
      print *, 'EXIT RUN_BK_3'
      
      END SUBROUTINE
      
      SUBROUTINE RUN_BK_5()
      
                                ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS,       IDBCO2
      USE DIAG04_MOD,    ONLY : AD04,          ND04
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT
      USE TIME_MOD,      ONLY : GET_TS_CHEM,  GET_DAY_OF_YEAR, 
     &     GET_MONTH
      
#     include "CMN_SIZE"      ! Size parameters

      INTEGER                :: I,  DOY, J,    L,     N
      REAL*8                 :: A_CM2, DTSRCE, E_CO2
      
      INTEGER                :: MONTH
      DTSRCE = 60d0 * GET_TS_CHEM()
      DOY    = GET_DAY_OF_YEAR()
      MONTH  = GET_MONTH()
      
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, A_CM2, E_CO2 )
      
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR
            
!-------------------------------------------
! #1: Total CO2
!     #2: CO2 from fossil fuel emissions 
!-------------------------------------------
            IF ( LFOSSCO2 ) THEN

! Fossil fuel emissions of CO2 [molec/cm2/s]
               E_CO2          = EMFOSSCO2(I,J)
               IF (DO_FF_SMOOTH) E_CO2  = E_CO2+ 
     &              EMFOSSCO2_NEXT(I,J)*(MONTH-1)

! ND04 diag: Fossil Fuel CO2 [molec/cm2/s] 
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,1) = AD04(I,J,1) + FLUX_FACTOR*E_CO2
               ENDIF

! Convert from [molec/cm2/s] to [kg]
               
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
            
            ENDIF
         
!-------------------------------------------
! #3: CO2 from ocean emissions
!-------------------------------------------
            IF ( LOCCO2 ) THEN
               
!     Ocean CO2 emissions in [molec/cm2/s]
               E_CO2          = EMOCCO2(I,J)
               
! ND04 diag: Ocean CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,2) = AD04(I,J,2) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert from [molec/cm2/s] to [kg]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
! Add to Tracer #1: Total CO2 [kg]
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
            ENDIF
            
!-------------------------------------------
! #4: CO2 from balanced biosphere emissions
!-------------------------------------------
            IF ( LBIOCO2 ) THEN
               
! Balanced biosphere CO2 [molec/cm2/s]
               E_CO2         = EMBIOCO2(I,J)
               
! ND04 diag: Bal Bio CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,3) = AD04(I,J,3) + FLUX_FACTOR*E_CO2
               ENDIF
               
! Convert from [molec/cm2/s] to [kg CO2]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
               
            ENDIF

!-------------------------------------------
! #5: CO2 from biomass burning emissions
!-------------------------------------------
            IF ( LBIOBRNCO2 ) THEN
! Biomass burning emissions [molec/cm2/s]
               E_CO2          = EMBIOMASSCO2(I,J)
               
! ND04 diag: Biomass burning CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,4) = AD04(I,J,4) + FLUX_FACTOR*E_CO2
               ENDIF
               
!     Convert from [molec/cm2/s] to [kg]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
            
            ENDIF
         
!-------------------------------------------
! #6: CO2 from biofuel emissions
!-------------------------------------------
            IF ( LBIOFUELCO2 ) THEN
               
!     Biofuel CO2 emissions [molec/cm2/s]
               E_CO2          = EMBIOFUELCO2(I,J)
               
! ND04 diag: Biofuel CO2 [molec/cm2/s] 
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,5) = AD04(I,J,5) + E_CO2
               ENDIF
               
! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
               STT(I,J,1,1)   = STT(I,J,1,1) + FLUX_FACTOR*E_CO2
               
            ENDIF
         
         
!-------------------------------------------
! #7: CO2 from posterior net terrestrial exchange
!-------------------------------------------
            
            IF ( LPOST ) THEN
               E_CO2          = POST_FLUX(I,J,1)
               
               IF ( ND04 > 0 ) THEN
!     save purterbation to bio-sphere IN MOLC CO2/CM2/S
                  AD04(I,J,3) = AD04(I,J,3) + 
     &                 E_CO2
               ENDIF
! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
               
!     ND04 diag: Bal Bio CO2 [molec/cm2/s]
               
!     Convert from [molec/cm2/s] to [kg]
               STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
               
            ENDIF
            
         ENDDO
      ENDDO

      
!$OMP END PARALLEL DO
      print *, 'EXIT RUN_BK_5'
      
      
      END SUBROUTINE



! THE CODE IS WRITTEN TO TAG ENKF TRACERS. THE DIFFERENCE FROM PREVIOUS VERSIONS 
! IS THAT WE ARE NOW ABLE TO USE HOURLY FLUXES 
      
      SUBROUTINE RUN_ENKF()

!     References to F90 modules
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT
      USE GRID_MOD,      ONLY : GET_AREA_CM2
#     include "CMN_SIZE"        ! Size parameters
      INTEGER                :: TRACER_ID
      INTEGER                :: I, J,  N
      REAL*8                 :: A_CM2, DTSRCE
      
      DTSRCE = 60d0 * GET_TS_CHEM()
      DTSRCE=DTSRCE / XNUMOL_CO2 
      
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N,TRACER_ID, A_CM2)
      
      DO N=1, N_TRACERS
         
         
!     THE FIRST TRACE CORRESPONDING TO NO-EMISSION FLUXES

         IF (N==1) THEN
            PRINT *, 'THE FIRST TRACER--PASS'
         ELSE
            TRACER_ID=ENR_BFST+N-1

! print *, 'TRACER_ID',  TRACER_ID
            
            DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT(I,J,1,N) = STT(I,J,1,N) +
     $                 DTSRCE*A_CM2*BF_FLUX(I, J,TRACER_ID)
               END DO
            ENDDO
         END IF
      END DO                          
      
!$OMP END PARALLEL DO
      print *,'EXIT ENKF RUN'
      END SUBROUTINE
      
      
!-----------------------------------------------------------------------------

      SUBROUTINE READ_ANNUAL_FOSSILCO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_FOSSILCO2 reads in annual mean fossil CO2 emissions 
!  from a binary punch file. (pns, bmy, 8/16/05, 10/3/05)
!
!  References:
!  ============================================================================
!  (1 ) CDIAC gridded (1x1) dataset for 1995 (Andres et al.)
!
!  NOTES:
!  (1 ) Emissions read in from directory : DATA_DIR/CO2XXX
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,   ONLY : GET_YEAR
      
#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      ! INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY

      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      INTEGER                :: YEAR, YEAR2
      character(len=4)       ::SYEAR
      character(len=2)       ::SMONTH
      character(len=2)       ::SDAY
      
      YEAR   = GET_YEAR()
      
      IF (YEAR>2010) THEN 
         WRITE( SYEAR, '(i4.4)' ) 2010 
      ELSE
         WRITE( SYEAR, '(i4.4)' ) YEAR
      END IF
!     WRITE( SMONTH, '(I2.2)') MONTH
      ! WRITE( SMONTH, '(I2.2)') MONTH
      ! SYEAR='2006'
      
      SMONTH='01'
      SDAY='01'
      
      
      FILENAME = TRIM(RUN_DIR )                      //
     & TRIM(EMS_PATH)// 
     &           'fossil.' //
     &           SYEAR //"."     //
     &           GET_NAME_EXT_2D() // '.'            //
     &           GET_RES_EXT()
     
      ! Echo info
       !=================================================================
      ! READ_ANNUAL_FOSSILCO2 begins here!
      !=================================================================

      ! File contaning fossil fuel CO2 data

      
      
     
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_FOSSILCO2: Reading ', a ) 

      ! Get TAU value corresponding to DOY in year 2000
      TAU = GET_TAU0(1, 1, 1985)
      
      ! Read fossil fuel CO2 [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), EMFOSSCO2 )
      ! ADD PURTERBATION
      EMFOSSCO2=PB_FF_FACTOR*EMFOSSCO2
      
      
      IF (DO_FF_SMOOTH) THEN
         
         YEAR2=YEAR+1
         IF (YEAR2>2006) THEN 
            WRITE( SYEAR, '(i4.4)' ) 2006 
         ELSE
            WRITE( SYEAR, '(i4.4)' ) YEAR2
         END IF
!     
         FILENAME = TRIM(RUN_DIR )                      //
     &        TRIM(EMS_PATH)// 
     &        'fossil.' //
     &        SYEAR //"."     //
     &        GET_NAME_EXT_2D() // '.'            //
     &        GET_RES_EXT()
         
         
         
         
         WRITE( 6, 100 ) TRIM( FILENAME )
         
         CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &        TAU,       IIPAR,     JJPAR,     
     &        1,         ARRAY,     QUIET=.TRUE. )
         
         
         CALL TRANSFER_2D( ARRAY(:,:,1), EMFOSSCO2_NEXT )
         EMFOSSCO2_NEXT =EMFOSSCO2_NEXT - EMFOSSCO2
         EMFOSSCO2_NEXT=EMFOSSCO2_NEXT/12.0 ! MONTHLY CHANGE RATE 
      END IF
      

      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_FOSSILCO2

!-----------------------------------------------------------------------------

      SUBROUTINE READ_MONTHLY_OCEANCO2( YEAR, MONTH, DAY)
!     
!******************************************************************************
!     Subroutine READ_ANNUAL_OCEANCO2 reads in annual mean oceanic CO2 exchange  
!     from a binary punch file. (pns, bmy, 8/16/05, 10/3/05)
!     
!     References:
!     ============================================================================
!     (1 ) Takahashi et al. (1997)
!     
!     NOTES:
!     (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!     
                                ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      
#     include "CMN_SIZE"        ! Size parameters
      
                                ! Local variables
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      
      character(len=4)       ::SYEAR
      character(len=2)       ::SMONTH
      character(len=2)       ::SDAY
      
      IF (YEAR>2010) THEN
         WRITE( SYEAR, '(i4.4)' ) 2010
      ELSE
         WRITE( SYEAR, '(i4.4)' ) YEAR 
      END IF 

      WRITE( SMONTH, '(I2.2)') MONTH
      WRITE( SMONTH, '(I2.2)') MONTH
      
      FILENAME = TRIM(RUN_DIR ) //
     &        TRIM(EMS_PATH)// 
     &     'ocean_co2.' //
     &     SYEAR//SMONTH//"." //
     &     GET_NAME_EXT_2D() // '.'  //
     &     GET_RES_EXT()      
      
!=================================================================
! READ_ANNUAL_OCEANCO2 begins here!
!=================================================================

      
! Echo info

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_OCEANCO2: Reading ', a ) 
      
! TAU value for start of "generic" year 1985 
      TAU = GET_TAU0( MONTH, 1, YEAR )
      print *, 'TAU', TAU
! Read ocean CO2 data [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &     TAU,       IIPAR,     JJPAR,     
     &     1,         ARRAY,     QUIET=.TRUE. )
      
 ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMOCCO2 )
      
                                ! Return to calling program
      END SUBROUTINE READ_MONTHLY_OCEANCO2
      

!------------------------------------------------------------------------------

      SUBROUTINE  READ_DAILY_BIOMASSCO2(YEAR, MONTH, DAY)
!
!******************************************************************************
!  Subroutine READ_BBIO_DIURNALCYCLE reads CASA daily Net Ecosystem Production
!  (NEP) fluxes but with a diurnal cycle imposed.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month  (1-31)
!  (3 ) HOUR  (INTEGER) : Current hour of day   (0-23)
!
!  References:
!  ============================================================================
!  (1 ) Olsen et al. [2004]. Fluxes are Net Ecosystem Production (NEP) 
!        from CASA model
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D 
      
      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY 

      ! Local variables
      INTEGER                :: DOY, TMP_HH
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3 )      :: SDOY
      CHARACTER(LEN=255)     :: FILENAME
      
      character(len=4)       ::SYEAR
      character(len=2)       ::SMONTH
      character(len=2)       ::SDAY
      
      IF (YEAR>2010) THEN 
         if (mod(year, 4)==0) then
            WRITE( SYEAR, '(i4.4)' ) 2008
         else
            WRITE( SYEAR, '(i4.4)' ) 2010
         end if
      ELSE
         WRITE( SYEAR, '(i4.4)' ) YEAR 
      END IF
      
      WRITE( SMONTH, '(I2.2)') MONTH
      WRITE( SDAY, '(I2.2)') DAY
      
      FILENAME = TRIM(RUN_DIR )   //
     &        TRIM(EMS_PATH)// 
     &           'casa_fire.' //
     &           SYEAR//SMONTH //"."     //
     &           GET_NAME_EXT_2D() // '.'            //
     &           GET_RES_EXT()
     
      !=================================================================
      ! READ_BBIO_DIURNALCYCLE begins here!
      !=================================================================
      
      ! Create string for day of year
      
      ! TMP_HH=HOUR/3
      TMP_HH=0
      
      
      

          
         
      
      
      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DIURNALCYCLE: Reading ', a )
      
      ! Get TAU of this month & day in "generic" year 1985
      IF (YEAR>2010) THEN 
         if (mod(year, 4)==0) then 
            TAU = GET_TAU0( MONTH, DAY, 2008, TMP_HH)
         else
            TAU = GET_TAU0( MONTH, DAY, 2010, TMP_HH)
         end if 
      ELSE
         TAU = GET_TAU0( MONTH, DAY, YEAR, TMP_HH) 
      END IF
      
         
      ! Read Net Ecosytem Productivity [molec CO/cm2/s] from disk
      ! The CASA fluxes use atmospheric convention: 
                                ! positive = into atm; negative = into biosphere
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
       
      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOMASSCO2)
      
      ! Return to calling program
      END SUBROUTINE READ_DAILY_BIOMASSCO2
      


!-----------------------------------------------------------------------------

      SUBROUTINE READ_MONTHLY_BIOMASSCO2(YEAR, MONTH, DAY)
!     
!******************************************************************************
!     Subroutine READ_ANNUAL_OCEANCO2 reads in annual mean oceanic CO2 exchange  
!     from a binary punch file. (pns, bmy, 8/16/05, 10/3/05)
!     
!     References:
!     ============================================================================
!     (1 ) Takahashi et al. (1997)
!     
!     NOTES:
!     (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!     
                                ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      
#     include "CMN_SIZE"        ! Size parameters
      
                                ! Local variables
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      
      character(len=4)       ::SYEAR
      character(len=2)       ::SMONTH
      character(len=2)       ::SDAY
      
      IF (YEAR>2006) THEN 
         WRITE( SYEAR, '(i4.4)' ) 2006 
      ELSE
         WRITE( SYEAR, '(i4.4)' ) YEAR 
      END IF
      
      WRITE( SMONTH, '(I2.2)') MONTH
      WRITE( SMONTH, '(I2.2)') MONTH
      
      FILENAME = TRIM(RUN_DIR ) //
     &     TRIM(EMS_PATH)// 
     &     SYEAR//SMONTH//"." //
     &     GET_NAME_EXT_2D() // '.' //
     &     GET_RES_EXT()
      
      
    !=================================================================
    ! READ_ANNUAL_OCEANCO2 begins here!
    !=================================================================

      
                                ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_MONTHLY_BIOMASSCO2: Reading ', a ) 
      
                                ! TAU value for start of "generic" year 1985 
      IF (YEAR>2006) THEN 
         TAU = GET_TAU0( MONTH, 1, 2006 )
      ELSE
         TAU = GET_TAU0( MONTH, 1, YEAR )
      END IF
      
                                ! Read ocean CO2 data [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'BIOSRCE', 4, 
     &     TAU,       IIPAR,     JJPAR,     
     &     1,         ARRAY,     QUIET=.TRUE. )
      
                                ! Cast to REAL*8 and resize if necessary
      print *,'AFTER READ BIO'

      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOMASSCO2 )
      
      EMBIOMASSCO2=EMFACTCO2CO*EMBIOMASSCO2
      
                                ! Return to calling program
      END SUBROUTINE READ_MONTHLY_BIOMASSCO2
      


      SUBROUTINE READ_ANNUAL_BIOFUELCO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_BIOFUELCO2 reads in annual mean biofuel CO2 
!  emissions from a binary punch file (pns, bmy, 8/16/05, 10/3/05)
!
!  References:
!  ============================================================================
!  (1 ) Yevich and Logan 2001 gridded (1x1) dataset in combination with 
!        emission factors for CO2 per kg drymatter burned
!  (2 ) See routines in /users/trop/pns/GEOSCHEM/EMISSIONS/BIOFUEL
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_BIOFUELCO2 begins here!
      !=================================================================

      ! File contaning biofuel CO2 data
      FILENAME = TRIM( DATA_DIR )          // 
     &           'CO2_200508/biofuel_CO2.' // GET_NAME_EXT_2D() //
     &           '.'                       // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_BIOFUELCO2: Reading ', a ) 
      TAU = GET_TAU0(1, 1, 1985)
      
      ! Read biofuel CO2 emissions [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 5, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOFUELCO2 )

      
      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_BIOFUELCO2


      SUBROUTINE READ_BBIO_DAILYAVERAGE( YEAR, MONTH, DAY, DOY ) 
!
!******************************************************************************
!  Subroutine READ_DAILY_BBIO_CO2 reads in daily values for balanced 
!  biospheric exchange from a binary punch file.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month (1-31)
!  (3 ) DOY   (INTEGER) : Current day of year (0-366)
!
!  Data Sources:
!  ============================================================================
!  (1 ) CASA gridded (1x1) dataset for from M. Thompson
!        Monthly values interpolated to daily values : 365 daily files 
!        NB : These files DO NOT have the diurnal cycle in daily emissions
!        See routine ' ' to read in files with diurnal cycle imposed
!
!  References
!  ============================================================================
!  (1 ) Randerson et al. [1997]
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      
#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY, DOY

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3  )     :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      INTEGER                :: TMP_DOY
      
      !=================================================================
      ! READ_BBIO_DAILYAVERAGE begins here!
      !=================================================================

      ! Make a string from DOY
      
      TMP_DOY=DOY
      
      IF (MOD(YEAR, 4).NE.0) THEN
         
         IF (MONTH>2) TMP_DOY=TMP_DOY+1
      END IF
      

          
         
      
      WRITE( SDOY, '(i3.3)' ) TMP_DOY 
      
      ! Name of file with Balanced Bio CO2 data
      FILENAME = TRIM( DATA_DIR )                      //
     &           'CO2_200508/BBIO_DAILYAVG/CO2.daily.' //
     &           GET_NAME_EXT_2D() // '.'              //
     &           GET_RES_EXT()     // '.'              // SDOY

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DAILYAVERAGE: Reading ', a ) 

      ! Get TAU value corresponding to DOY in year 2000
      TAU = GET_TAU0( MONTH, DAY, 2000 )
      print *, 'TAU', TAU
      
     

      ! Read balanced biosphere CO2 [molec/cm2/s] from disk
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 3, 
     &                 TAU,       IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )


      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )

      ! Return to calling program
      END SUBROUTINE READ_BBIO_DAILYAVERAGE


      SUBROUTINE READ_CASA_BIO( YEAR, MONTH, DAY, HOUR )
!
!******************************************************************************
!  Subroutine READ_BBIO_DIURNALCYCLE reads CASA daily Net Ecosystem Production
!  (NEP) fluxes but with a diurnal cycle imposed.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month  (1-31)
!  (3 ) HOUR  (INTEGER) : Current hour of day   (0-23)
!
!  References:
!  ============================================================================
!  (1 ) Olsen et al. [2004]. Fluxes are Net Ecosystem Production (NEP) 
!        from CASA model
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D 
      
      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY, HOUR 

      ! Local variables
      INTEGER                :: DOY, TMP_HH
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3 )      :: SDOY
      CHARACTER(LEN=255)     :: FILENAME
      
      character(len=4)       ::SYEAR
      character(len=2)       ::SMONTH
      character(len=2)       ::SDAY
      
      IF (YEAR>2010) THEN 
         if (mod(year, 4)==0) then
            WRITE( SYEAR, '(i4.4)' ) 2008
         else
            WRITE( SYEAR, '(i4.4)' ) 2010
         end if
      ELSE
         WRITE( SYEAR, '(i4.4)' ) YEAR 
      END IF
      
      WRITE( SMONTH, '(I2.2)') MONTH
      WRITE( SDAY, '(I2.2)') DAY
      
      FILENAME = TRIM(RUN_DIR )   //
     &     TRIM(EMS_PATH)// 
     &     'casa_co2.'//
     &     SYEAR//SMONTH //"." //
     &     GET_NAME_EXT_2D() // '.' //
     &     GET_RES_EXT()
      
     
      !=================================================================
      ! READ_BBIO_DIURNALCYCLE begins here!
      !=================================================================
      
      ! Create string for day of year
      
      TMP_HH=HOUR/3
      TMP_HH=3*TMP_HH
      
      

          
         
      
      
      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DIURNALCYCLE: Reading ', a )

      ! Get TAU of this month & day in "generic" year 1985
      IF (YEAR>2010) THEN 
         if (mod(year, 4)==0) then 
            TAU = GET_TAU0( MONTH, DAY, 2008, TMP_HH)
         else
            TAU = GET_TAU0( MONTH, DAY, 2010, TMP_HH)
         end if 
      ELSE
         TAU = GET_TAU0( MONTH, DAY, YEAR, TMP_HH) 
      END IF
      
         
      ! Read Net Ecosytem Productivity [molec CO/cm2/s] from disk
      ! The CASA fluxes use atmospheric convention: 
                                ! positive = into atm; negative = into biosphere
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      
      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )
         
      ! Return to calling program
      END SUBROUTINE READ_CASA_BIO
      
      
!------------------------------------------------------------------------------

      SUBROUTINE READ_BBIO_DIURNALCYCLE( YEAR, MONTH, DAY, HOUR )
!
!******************************************************************************
!  Subroutine READ_BBIO_DIURNALCYCLE reads CASA daily Net Ecosystem Production
!  (NEP) fluxes but with a diurnal cycle imposed.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month  (1-31)
!  (3 ) HOUR  (INTEGER) : Current hour of day   (0-23)
!
!  References:
!  ============================================================================
!  (1 ) Olsen et al. [2004]. Fluxes are Net Ecosystem Production (NEP) 
!        from CASA model
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      USE TIME_MOD,      ONLY : GET_DAY,       GET_DAY_OF_YEAR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D 

      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: YEAR, MONTH, DAY, HOUR 

      ! Local variables
      INTEGER                :: DOY, TMP_DOY
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3 )      :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_BBIO_DIURNALCYCLE begins here!
      !=================================================================
      
      ! Create string for day of year
      DOY    = GET_DAY_OF_YEAR()
      TMP_DOY=DOY
      
      IF (MOD(YEAR, 4).Eq.0) THEN
         
         IF (MONTH>2) TMP_DOY=TMP_DOY-1
      
      END IF
      
      
          
         
      
      WRITE( SDOY, '(i3.3)' ) TMP_DOY 
      
      
! if the 4x5 is used please switch the following on 
      
      FILENAME = TRIM(RUN_DIR ) //
     &     TRIM(EMS_PATH)//     
     &     'nep.' // 
     &     GET_NAME_EXT_2D() // '.' //
     &     GET_RES_EXT()     // '.' // SDOY
      
      
                                ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DIURNALCYCLE: Reading ', a )

      ! Get TAU of this month & day in "generic" year 1985
      TAU = GET_TAU0( MONTH, DAY, 1985, HOUR ) 
!     TAU = GET_TAU0( MONTH, DAY, YEAR, HOUR ) 
      
      ! Read Net Ecosytem Productivity [molec CO/cm2/s] from disk
      ! The CASA fluxes use atmospheric convention: 
      ! positive = into atm; negative = into biosphere
      ! print *, 'TAU', TAU
      
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP', 2, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
       
      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )
         
      ! Return to calling program
      END SUBROUTINE READ_BBIO_DIURNALCYCLE



      
      
      SUBROUTINE INIT_CO2 
!
!******************************************************************************
!  Subroutine INIT_CO2 allocates memory to module arrays and reads in annual
!  mean emissions. (pns, bmy, 8/16/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS 

      !=================================================================
      ! INIT_CO2 begins here!
      !=================================================================
           
      ! Exit if we have already intialized 
      IF ( IS_INIT ) RETURN

      ! callocate BB_REGION
      
      
      ! Array for Fossil fuel CO2
      ALLOCATE( EMFOSSCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMFOSSCO2' )         
      EMFOSSCO2 = 0d0 

      ! Array for Fossil fuel CO2
      ALLOCATE( EMFOSSCO2_NEXT( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMFOSSCO2' )         
      EMFOSSCO2_NEXT = 0d0 
      
      ! Array for CO2 from ocean exchange
      ALLOCATE( EMOCCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMOCCO2' )         
      EMOCCO2 = 0d0 

      ! Array for Balanced Bio CO2
      ALLOCATE( EMBIOCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOCO2' )         
      EMBIOCO2 = 0d0 

      ! Array for Biomass burning CO2
      ALLOCATE( EMBIOBRNCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOBRNCO2' )
      EMBIOBRNCO2 = 0d0

      ! Array for Biofuel CO2
      ALLOCATE( EMBIOFUELCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOFUELCO2' )         
      EMBIOFUELCO2 = 0d0 

      ! Array for NET BIO CO2
      ALLOCATE( EMBIONETCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIONETCO2' )
      EMBIONETCO2  = 0d0


      ! Array for NET BIO CO2
      ALLOCATE( EMBIOMASSCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIONETCO2' )
      EMBIOMASSCO2  = 0d0
      
      
      
      
      
      ALLOCATE(BF_FLUX( IIPAR, JJPAR,MAX_PBFLUX), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIONETCO2' )
      BF_FLUX  = 0d0
      

      
      

      
      !=================================================================
      ! Read in annual mean emissions
      !=================================================================

      IF ( LFOSSCO2    ) CALL READ_ANNUAL_FOSSILCO2
      
      ! Fossil fuel emissions
      
      ! Oceanic exchange
      ! DO IT MONTHLY 
!     
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_CO2 
      
      SUBROUTINE INIT_CO2_ENKF(YYYYS, DOYS, PBUSE, NTIME, ENST, ENEND, 
     &     PBSTART,   
     &     CO2FILE, DOBKRUN)
      ! INITIALIZE  TAGGED CO2 ENKF 
      ! NTIME           IN  THE NUMBER OF ENSEMBLE 
      ! YYYYS(NITIME)   IN  THE YEAR OF THE STARTING DATES FOR SURFACE EMISSIONS
      ! DOYS(NTIME)     IN  THE DOY OF THE STARTING DATES FOR SURFACE EMISSIONS
      ! ENST            IN  THE FIRST MEMBER OF ENSEMBLE 
      ! ENEND           IN  THE LAST MEMBER OF ENSEMBLE 
      ! DOBKRUN         IN  CHOOSE TO TAG THE BKGROUND OR TAGGED TRACER
      INTEGER, INTENT(IN)::NTIME  
      INTEGER, INTENT(IN)::YYYYS(NTIME)
      INTEGER, INTENT(IN)::DOYS(NTIME), PBUSE(NTIME)
      INTEGER, INTENT(IN)::ENST, ENEND, PBSTART, DOBKRUN
      CHARACTER(LEN=*), INTENT(IN)::CO2FILE
      LOGICAL, SAVE::IS_INIT=.FALSE.
      INTEGER::I
      
      IF (IS_INIT) RETURN
      IS_INIT=.TRUE.
      
      ENR_FILE=TRIM(CO2FILE)
      ALLOCATE(ENR_YYYY(NTIME))
      ALLOCATE(ENR_DOY(NTIME))
      ALLOCATE(ENR_DD(NTIME))
      ALLOCATE(ENR_PBUSE(NTIME))
      DO I=1, NTIME
         ENR_YYYY(I)=YYYYS(I)
         ENR_DOY(I)=DOYS(I)
         ENR_PBUSE(I)=PBUSE(I)
         ENR_DD(I)=(YYYYS(I)-1985.0)*365.24+DOYS(I)
      END DO
      ENR_NTIME=NTIME
      ENR_ST=ENST
      ENR_END=ENEND
      ENR_BFST= PBSTART
      DO_BK_RUN=DOBKRUN
      
      END SUBROUTINE INIT_CO2_ENKF

      

! -----------------------------------------------------------------------------
      
      SUBROUTINE READ_BASIS_FUNCTION(YEAR, MONTH, DAY, HOUR)
!     
!******************************************************************************
!  Subroutine emission ensemble

!  References:
!  ============================================================================
!  (1 ) (feng et al., 2009)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,  READ_BPCH2_BF
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      
#     include "CMN_SIZE"      ! Size parameters
      
      ! Local variables
      INTEGER,     INTENT(IN) :: YEAR, MONTH, DAY, HOUR
      INTEGER                :: I, J
      CHARACTER(LEN=255)     :: FILENAME
      INTEGER                :: TMP_HH 
      REAL*8                 :: TAU
      
      
      TMP_HH=HOUR/3
      TMP_HH=3*TMP_HH
      
      
      ! Echo file name to stdout

      ! Get TAU of this month & day in "generic" year 1985
      TAU = GET_TAU0( MONTH, DAY, YEAR, TMP_HH)
            
      
      
      ! ============================================================
      ! READ_CO2_BF begins here!
      !=================================================================

      ! File containing Basis functions 
      
      FILENAME = TRIM( RUN_DIR ) //
     &     TRIM(EMS_PATH)// 
     &     TRIM(ENR_FILE)//"."//
     $     GET_NAME_EXT_2D() //"."
     $     //GET_RES_EXT()
      
                                ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME)
 100  FORMAT( '     - READ BASIS FUNCTION- ', a ) 
      
      ! Read fossil fuel CO2 [molec/cm2/s]
      BF_FLUX=0.0
      
      CALL READ_BPCH2_BF(FILENAME, 
     &     IGLOB, 
     &     JGLOB, MAX_PBFLUX, TAU, BF_FLUX, NBF)
      
      PRINT *, 'PB FLUX, MAX, MIN, TOTAL (KG/S)'
      PRINT *,  'NBF', NBF, MAXVAL(BF_FLUX), MINVAL(BF_FLUX)
      
                                ! Return to calling program
      END SUBROUTINE READ_BASIS_FUNCTION
      


      SUBROUTINE READ_POST_FLUX(YEAR, MONTH, DAY)

!******************************************************************************
!  Subroutine to read posterior fluxes 

!  References:
!  ============================================================================
!  (1 ) (feng et al., 2009)
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : READ_BPCH2_REG_FLUX
      USE DIRECTORY_MOD, ONLY : DATA_DIR, RUN_DIR
      
#     include "CMN_SIZE"      ! Size parameters
      
      ! Local variables
      INTEGER,     INTENT(IN) :: YEAR
      INTEGER,     INTENT(IN) :: MONTH, DAY
      
      CHARACTER(LEN=8)       :: SDATE
      CHARACTER(LEN=4)      :: SYYYY
      CHARACTER(LEN=2)       ::SMONTH
      CHARACTER(LEN=2)       ::SDAY
      
      INTEGER                :: I, J
      CHARACTER(LEN=255)     :: FILENAME
      
      
      WRITE(SYYYY, '(i4.4)')  YEAR
      WRITE(SMONTH, '(i2.2)') MONTH
      WRITE(SDAY, '(i2.2)')   DAY
      SDATE=SYYYY//SMONTH//SDAY
      
      IF (.NOT.ALLOCATED(POST_FLUX)) 
     $     ALLOCATE(POST_FLUX(IGLOB, JGLOB, 1))
      
      
      !=================================================================
      ! READ_POST_FLUX begins here!
      !=================================================================
      
                                ! File contaning fossil fuel CO2 data
      FILENAME = TRIM( RUN_DIR ) //
     &        TRIM(EMS_PATH)// 
     & "post_flux."//SDATE//"."// 
     &     GET_NAME_EXT_2D() //"."
     $     //GET_RES_EXT()
     
                                ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - POSTER_CO2_FLUX2: Reading ', a ) 
      
      ! Read fossil fuel CO2 [molec/cm2/s]
      CALL READ_BPCH2_REG_FLUX(FILENAME, 
     &     IGLOB, 
     &     JGLOB, 1, POST_FLUX)
      
                                ! Return to calling program
      END SUBROUTINE READ_POST_FLUX
      
      SUBROUTINE GET_YYYY_DOY(DD, R_YYYY,R_DOY)
      ! FIND THE START DATE (R_YYYY & R_DOY) OF THE EMISSIONS
      ! DD 
      
      REAL*8, INTENT(IN) :: DD
      INTEGER, INTENT(OUT) :: R_DOY, R_YYYY
      INTEGER ::I
      
      DO I=2, ENR_NTIME
         IF (DD>=ENR_DD(I-1).AND.DD< ENR_DD(I)) EXIT
      END DO
      
      R_DOY=ENR_DOY(I-1)
      R_YYYY=ENR_YYYY(I-1)
      END SUBROUTINE GET_YYYY_DOY
      
!!!       
      
!------------------------------------------------------------------------------
  
      SUBROUTINE CLEANUP_CO2 
!
!******************************************************************************
!  Subroutine CLEANUP_CO2 deallocates all module arrays (pns, bmy, 8/16/05)
! 
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_CO2 begins here!
      !=================================================================
      IF ( ALLOCATED( EMFOSSCO2    ) ) DEALLOCATE( EMFOSSCO2    )
      IF ( ALLOCATED( EMFOSSCO2_NEXT    ) ) DEALLOCATE( EMFOSSCO2_NEXT)
      
      IF ( ALLOCATED( EMOCCO2      ) ) DEALLOCATE( EMOCCO2      )
      IF ( ALLOCATED( EMBIOCO2     ) ) DEALLOCATE( EMBIOCO2     )
      IF ( ALLOCATED( EMBIOBRNCO2  ) ) DEALLOCATE( EMBIOBRNCO2  )
      IF ( ALLOCATED( EMBIOFUELCO2 ) ) DEALLOCATE( EMBIOFUELCO2 )
      IF ( ALLOCATED( EMBIOMASSCO2 ) ) DEALLOCATE( EMBIOMASSCO2 )
      
      ! added by lf to remove the memory for ENKF simulation 
      IF (ALLOCATED(ENR_YYYY)) DEALLOCATE(ENR_YYYY)
      IF (ALLOCATED(ENR_DOY))  DEALLOCATE(ENR_DOY)
      IF (ALLOCATED(ENR_DD))  DEALLOCATE(ENR_DD)
      IF (ALLOCATED(POST_FLUX)) DEALLOCATE(POST_FLUX)
      IF (ALLOCATED(BF_FLUX)) DEALLOCATE(BF_FLUX)
      
      
      
      
      
      !     

      ! Return to calling program
      END SUBROUTINE CLEANUP_CO2

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CO2_MOD
