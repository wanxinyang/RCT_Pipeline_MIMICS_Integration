C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************     MIMICS Version 1.5a MODEL              ***************
C***********************************************************************
C*************           Master Program                  ***************
C***********************************************************************
C
      PROGRAM MIMICS_version_1_5a
        save


C
C***********************************************************************
C------------------- Development begun 10-5-88 -------------------------
C--------------------  revised 10-5-88  --------------------------------
C***********************************************************************
C   DECLARE THE PARAMETERS IN THE 'PARAMETERS.INCLUDE' FILE.
C----------------------------
c%include 'parameters.include'
         include 'parameters.include'
C----------------------------
C***********************************************************************
C-----------------------------------------------------------------------
C--------------------- VARIABLE DECLARATIONS ---------------------------
C-----------------------------------------------------------------------
C***********************************************************************
C
        LOGICAL COMPUTE, WRITE_HEADER
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
        REAL FREQ_HERTZ, WAVELENGTH, k0, THETA, CTHETA, STHETA
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
C
C***********************************************************************
C   DECLARATIONS FROM SUBROUTINE INITIALIZE
C***********************************************************************
C
        INTEGER LOOP_NUM(N_VARIABLES), LOOP_COUNT(N_VARIABLES)
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL CALL_SUB(N_CALLS,N_SUB_CALLS)
C
        LOGICAL LOG_CROWN
C
        COMMON /I_COUNT/ LOOP_NUM, LOOP_COUNT
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
C***********************************************************************
C   DECLARATIONS FROM SUBROUTINE READ_INPUT
C***********************************************************************
C
C-------- DECLARATIONS FROM SUBROUTINE READ_CONFIGURE ------------------
C
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
C-------- DECLARATIONS FROM SUBROUTINE READ_SENSOR ---------------------
C
        INTEGER FREQ_NUMBER, THETA_NUMBER
        REAL FREQ_DELTA, FREQ_VECT(I_VECT_LEN), FREQ_GHZ
        REAL THETA_DELTA, THETA_VECT(I_VECT_LEN), THETA_DEGREES
        LOGICAL LOG_FREQ_TABLE, LOG_THETA_TABLE
C
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
        COMMON /R_SENSOR_DELTA/ THETA_DELTA, FREQ_DELTA
        COMMON /I_SENSOR/ THETA_NUMBER, FREQ_NUMBER
        COMMON /R_SENSOR_VECT/ THETA_VECT, FREQ_VECT
        COMMON /L_SENSOR/ LOG_FREQ_TABLE,LOG_THETA_TABLE
C
C-------- DECLARATIONS FROM SUBROUTINE READ_ENVIRONMENT ----------------
C
        INTEGER T_SOIL_NUMBER, T_WATER_NUMBER, T_VEG_NUMBER
        REAL T_SOIL_DELTA, T_SOIL_VECT(I_VECT_LEN), T_SOIL
        REAL T_WATER_DELTA, T_WATER_VECT(I_VECT_LEN), T_WATER
        REAL T_VEG_VECT(I_VECT_LEN)
        REAL T_VEG_DELTA, T_VEG
        LOGICAL LOG_ENVIRONMENT(N_ENVIRONMENT)
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_ENVIRON_DELTA/ T_SOIL_DELTA,T_WATER_DELTA,T_VEG_DELTA
        COMMON /I_ENVIRON/ T_SOIL_NUMBER, T_WATER_NUMBER, T_VEG_NUMBER
        COMMON /R_ENVIRON_VECT/ T_SOIL_VECT, T_WATER_VECT, T_VEG_VECT
        COMMON /L_ENVIRON/ LOG_ENVIRONMENT
C
C-------- DECLARATIONS FROM SUBROUTINE READ_GROUND ---------------------
C
        INTEGER MV_SOIL_NUMBER, RMS_SOIL_NUMBER, LS_SOIL_NUMBER
        INTEGER SAND_NUMBER, CLAY_NUMBER, SALT_NUMBER
        REAL MV_SOIL_DELTA, RMS_SOIL_DELTA, LS_SOIL_DELTA
        REAL SAND_DELTA, CLAY_DELTA, SALT_DELTA
        REAL MV_SOIL_VECT(I_VECT_LEN), RMS_SOIL_VECT(I_VECT_LEN)
        REAL LS_SOIL_VECT (I_VECT_LEN)
        REAL SAND_VECT(I_VECT_LEN), CLAY_VECT(I_VECT_LEN)
        REAL SALT_VECT(I_VECT_LEN)
        REAL MV_SOIL, SAND, CLAY, RMS_SOIL, LS_SOIL
        REAL SALT
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
        LOGICAL LOG_MV_SOIL,LOG_ROUGH_SOIL(2),LOG_SOIL_TEXT,LOG_SALT
C
        COMMON /R_SURFACE/ RMS_SOIL, LS_SOIL
        COMMON /R_SOIL/ MV_SOIL, SAND, CLAY
        COMMON /R_WATER/ SALT
C
        COMMON /R_SURFACE_DELTA/ RMS_SOIL_DELTA, LS_SOIL_DELTA
        COMMON /R_SOIL_DELTA/ MV_SOIL_DELTA, SAND_DELTA, CLAY_DELTA
        COMMON /R_WATER_DELTA/ SALT_DELTA
C
        COMMON /I_SURFACE/ RMS_SOIL_NUMBER, LS_SOIL_NUMBER
        COMMON /I_SOIL/ MV_SOIL_NUMBER, SAND_NUMBER, CLAY_NUMBER
        COMMON /I_WATER/ SALT_NUMBER
C
        COMMON /R_SURFACE_VECT/ RMS_SOIL_VECT, LS_SOIL_VECT
        COMMON /R_SOIL_VECT/ MV_SOIL_VECT,SAND_VECT,CLAY_VECT
        COMMON /R_WATER_VECT/ SALT_VECT
C
        COMMON /L_SURFACE/ LOG_ROUGH_SOIL
        COMMON /L_SOIL/ LOG_MV_SOIL, LOG_SOIL_TEXT
        COMMON /L_WATER/ LOG_SALT
C
        COMMON /L_SURFACE_TYPE/ SOIL_SURFACE,WATER_SURFACE,
     &                          SNOW_SURFACE,ICE_SURFACE
C
C-------- DECLARATIONS FROM SUBROUTINE READ_TRUNK ----------------------
C
        INTEGER MG_TRUNK_NUM, RHO_TRUNK_NUM, DENSITY_NUM   
        INTEGER CROWN_HGHT_NUM, TRUNK_DIAM_NUM, TRUNK_HGHT_NUM
        REAL MG_TRUNK_DELTA, RHO_TRUNK_DELTA, DENSITY_DELTA
        REAL CROWN_HGHT_DELTA, TRUNK_DIAM_DELTA, TRUNK_HGHT_DELTA
        REAL MG_TRUNK_VECT(I_VECT_LEN), RHO_TRUNK_VECT(I_VECT_LEN)
        REAL DENSITY_VECT(I_VECT_LEN)
        REAL CROWN_HGHT_VECT(I_VECT_LEN), TRUNK_DIAM_VECT(I_VECT_LEN)
        REAL TRUNK_HGHT_VECT(I_VECT_LEN)
        REAL MG_TRUNK, RHO_TRUNK, DENSITY, CROWN_HGHT
        REAL TRUNK_DIAM,TRUNK_HGHT
        LOGICAL LOG_MG_TRUNK, LOG_RHO_TRUNK, LOG_DENSITY
        LOGICAL LOG_CROWN_HGHT, LOG_TRUNK_DIAM, LOG_TRUNK_HGHT
C
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_TRUNK_DELTA/ MG_TRUNK_DELTA, RHO_TRUNK_DELTA, 
     &                         TRUNK_DIAM_DELTA
        COMMON /I_TRUNK/ MG_TRUNK_NUM, RHO_TRUNK_NUM, TRUNK_DIAM_NUM
        COMMON /R_TRUNK_VECT/MG_TRUNK_VECT, RHO_TRUNK_VECT,
     &                         TRUNK_DIAM_VECT
        COMMON /L_TRUNK/ LOG_MG_TRUNK,LOG_RHO_TRUNK,LOG_TRUNK_DIAM
C
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_CANOPY_DELTA/ DENSITY_DELTA, CROWN_HGHT_DELTA,
     &                          TRUNK_HGHT_DELTA
        COMMON /I_CANOPY/ DENSITY_NUM, CROWN_HGHT_NUM, TRUNK_HGHT_NUM
        COMMON /R_CANOPY_VECT/DENSITY_VECT, CROWN_HGHT_VECT,
     &                        TRUNK_HGHT_VECT
        COMMON /L_CANOPY/ LOG_DENSITY, LOG_CROWN_HGHT, LOG_TRUNK_HGHT
C
C-------- DECLARATIONS FROM SUBROUTINE READ_LEAF -----------------------
C
        INTEGER MG_LEAF_NUM, RHO_LEAF_NUM, LEAF_DENS_NUM   
        INTEGER LEAF_DIAM_NUM, LEAF_TAU_NUM
        REAL MG_LEAF_DELTA, RHO_LEAF_DELTA, LEAF_DENS_DELTA
        REAL LEAF_DIAM_DELTA, LEAF_TAU_DELTA
        REAL MG_LEAF_VECT(I_VECT_LEN), RHO_LEAF_VECT(I_VECT_LEN)
        REAL LEAF_DENS_VECT(I_VECT_LEN)
        REAL LEAF_DIAM_VECT(I_VECT_LEN),LEAF_TAU_VECT(I_VECT_LEN)
        REAL MG_LEAF, RHO_LEAF, LEAF_DENS, LEAF_DIAM,LEAF_TAU
        LOGICAL LOG_MG_LEAF, LOG_RHO_LEAF, LOG_LEAF_DENS
        LOGICAL LOG_LEAF_DIAM, LOG_LEAF_TAU
C
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /R_LEAF_DELTA/ MG_LEAF_DELTA, RHO_LEAF_DELTA, 
     &                LEAF_DENS_DELTA, LEAF_DIAM_DELTA, LEAF_TAU_DELTA
        COMMON /I_LEAF/ MG_LEAF_NUM, RHO_LEAF_NUM, LEAF_DENS_NUM,
     &                  LEAF_DIAM_NUM, LEAF_TAU_NUM
        COMMON /R_LEAF_VECT/MG_LEAF_VECT, RHO_LEAF_VECT,
     &              LEAF_DENS_VECT,LEAF_DIAM_VECT, LEAF_TAU_VECT
        COMMON /L_LEAF/ LOG_MG_LEAF, LOG_RHO_LEAF,
     &              LOG_LEAF_DENS, LOG_LEAF_DIAM, LOG_LEAF_TAU
C
C-------- DECLARATIONS FROM SUBROUTINE READ_NEEDLE ---------------------
C
        INTEGER MG_NDL_NUM, RHO_NDL_NUM, NDL_DENS_NUM   
        INTEGER NDL_DIAM_NUM, NDL_LNG_NUM
        REAL MG_NDL_DELTA, RHO_NDL_DELTA, NDL_DENS_DELTA
        REAL NDL_DIAM_DELTA, NDL_LNG_DELTA
        REAL MG_NDL_VECT(I_VECT_LEN), RHO_NDL_VECT(I_VECT_LEN)
        REAL NDL_DENS_VECT(I_VECT_LEN)
        REAL NDL_DIAM_VECT(I_VECT_LEN), NDL_LNG_VECT(I_VECT_LEN)
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM,NDL_LNG
        LOGICAL LOG_MG_NDL, LOG_RHO_NDL, LOG_NDL_DENS
        LOGICAL LOG_NDL_DIAM, LOG_NDL_LNG
C
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /R_NDL_DELTA/ MG_NDL_DELTA,RHO_NDL_DELTA,NDL_DENS_DELTA,
     &                        NDL_DIAM_DELTA, NDL_LNG_DELTA
        COMMON /I_NDL/ MG_NDL_NUM, RHO_NDL_NUM, NDL_DENS_NUM,
     &                  NDL_DIAM_NUM, NDL_LNG_NUM
        COMMON /R_NDL_VECT/MG_NDL_VECT, RHO_NDL_VECT,NDL_DENS_VECT,
     &                      NDL_DIAM_VECT, NDL_LNG_VECT
        COMMON /L_NDL/ LOG_MG_NDL, LOG_RHO_NDL, LOG_NDL_DENS,
     &                      LOG_NDL_DIAM, LOG_NDL_LNG
C
C-------- DECLARATIONS FROM SUBROUTINE READ_PRIMARY_BRANCH -------------
C
        INTEGER MG_BR1_NUM, RHO_BR1_NUM, BR1_DENS_NUM   
        INTEGER BR1_DIAM_NUM, BR1_LNG_NUM
        REAL MG_BR1_DELTA, RHO_BR1_DELTA, BR1_DENS_DELTA
        REAL BR1_DIAM_DELTA, BR1_LNG_DELTA
        REAL MG_BR1_VECT(I_VECT_LEN), RHO_BR1_VECT(I_VECT_LEN)
        REAL BR1_DENS_VECT(I_VECT_LEN)
        REAL BR1_DIAM_VECT(I_VECT_LEN), BR1_LNG_VECT(I_VECT_LEN)
        REAL MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        LOGICAL LOG_MG_BR1, LOG_RHO_BR1, LOG_BR1_DENS
        LOGICAL LOG_BR1_DIAM, LOG_BR1_LNG
C
        COMMON /R_BR1/ MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        COMMON /R_BR1_DELTA/ MG_BR1_DELTA,RHO_BR1_DELTA,BR1_DENS_DELTA,
     &                        BR1_DIAM_DELTA, BR1_LNG_DELTA
        COMMON /I_BR1/ MG_BR1_NUM, RHO_BR1_NUM, BR1_DENS_NUM,
     &                  BR1_DIAM_NUM, BR1_LNG_NUM
        COMMON /R_BR1_VECT/MG_BR1_VECT, RHO_BR1_VECT, BR1_DENS_VECT,
     &                      BR1_DIAM_VECT, BR1_LNG_VECT
        COMMON /L_BR1/ LOG_MG_BR1, LOG_RHO_BR1, LOG_BR1_DENS,
     &                      LOG_BR1_DIAM, LOG_BR1_LNG
C
C-------- DECLARATIONS FROM SUBROUTINE READ_NESTING --------------------
C
        INTEGER PARAM_NUM(N_VARIABLES), INEST(N_VARIABLES)
        COMMON /I_NEST/ PARAM_NUM, INEST
C
        REAL GRND_REFLECT_MAT(4,4), GRND_BACK_MAT(4,4)
        COMMON /R_GROUND_MATS/ GRND_REFLECT_MAT, GRND_BACK_MAT
C
C-------- DECLARATIONS FROM SUBROUTINE TRUNK_LAYER ---------------------
C
        REAL EXP_KAPPA_T_p(4,4), EXP_KAPPA_T_m (4,4)
        COMMON /TRUNK_EXT/ EXP_KAPPA_T_p, EXP_KAPPA_T_m

        REAL TRUNK_PHASE_MAT_p(4,4,2), TRUNK_PHASE_MAT_m(4,4,2)
        COMMON /TRUNK_PHASE/ TRUNK_PHASE_MAT_p, TRUNK_PHASE_MAT_m

C-----------------------------------------------------------------------
C DECLARE THE CONSTANTS IN THE 'CONSTANTS.INCLUDE' FILE.
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C
C***********************************************************************
C***********************************************************************
C******************  MAIN DRIVER ROUTINE   *****************************
C***********************************************************************
C***********************************************************************
C
C-- Read input data, initalize fixed constants, set output files flag --
C
CC        print *,'calling read_input'
        CALL READ_INPUT
CC        print *,'calling init'
        CALL INITIALIZE(COMPUTE)
        i = 0
c
        DO WHILE(COMPUTE)
            i=i+1
c
c-----------   set constants based on input parmeters ------------------
c                these constants vary as a function of input parameters
c                i.e. wavelength and theta(radians)
C----  also check for looping which varies input parameters  -----------
c
CC            print *,'calling start_loop'
            CALL START_LOOP(WRITE_HEADER,COMPUTE,LOG_CROWN)

CC            DO i=1,N_CALLS   
CC               DO j=1,N_SUB_CALLS
CC                  print *,'call_sub(',i,',',j,')= ',call_sub(i,j)
CC               ENDDO
CC            ENDDO


c
c----------  compute dielectric constants  -----------------------------
c
            IF(CALL_SUB(1,1)) THEN
CC              print *,'calling dielecric'
              CALL DIELECTRIC(EPSILONR,EPSILONRC)
            ENDIF
C
c--------- compute ground specular reflectivity matrix and -------------
c---------         backscattering phase matrix             -------------
c                                                                       
            IF(CALL_SUB(2,1))  THEN
CC                 print *,'calling ground layer'
                 CALL GROUND_LAYER
            ENDIF
C
C  The numbers below are sigvv and sighh.
C
        PRINT*,' '
        PRINT*,'I = ',I
c
c----------  compute trunk layer phase and extinction matrices  --------
c
            IF(CALL_SUB(3,1)) THEN
CC               print *,'calling trunk layer'
               CALL TRUNK_LAYER(LOG_CONSTITUENT(1))
            ENDIF
c
c--------- compute crown layer phase and extinction matrices  ----------
c
            IF(CALL_SUB(4,1)) THEN
                print *,'calling crown_layer'
                CALL CROWN_LAYER(LOG_CROWN)
            ENDIF

c
c------------ compute stokes matrix of the canopy ----------------------
c
              print *,'calling solve_canopy'
              CALL SOLVE_CANOPY(LOG_CONSTITUENT(1),LOG_CROWN)
C
C----------- FORMAT DATA FOR DESIRED OUTPUT TYPE -----------------------
C
CC            print *,'calling format_output'
            CALL FORMAT_OUTPUT
C
C-----------         and output results          -----------------------
C
CC            print *,'calling write_data'
            CALL WRITE_DATA(WRITE_HEADER)
c
        ENDDO
c
c--------   clean up and close files  ----------------------------------
c
CC       print *,'calling finish'
       CALL FINISH
c
c-----------------------------------------------------------------------
c
        STOP
        END
