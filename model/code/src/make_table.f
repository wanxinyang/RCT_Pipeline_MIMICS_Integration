C make_table.f        10-30-92
C
C by Leland Pierce
C
C Prints out the order of requests to the parameter_value_table.input
C file during a read.
C Gives the variable name when the request is made.
C 
C This is a modification to source/read_input.f, the true MIMICS 
C reader subroutines file.
C
C usage: after having filled in every input file EXCEPT:
C        input_dir/parameter_value_table.input,
C        go to the directory "input_dir", then just run: make_table.
C        It will create the annotated file:
C          input_dir/parameter_value_table.input
C         with the parameter names, etc, and empty templates for
C         each.
C         You may then edit this file to specify these parameters.
C
C------------------------------------------------------------------------

      CALL READ_INPUT
 
      STOP
      END




C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to Read input data           ***************
C***    Calling routine:      main MIMICS program                    ***
C***    Called subroutines:   READ_CONFIGURE                         ***
C***                          READ_SENSOR                            ***
C***                          READ_TRUNK                             ***
C***                          READ_LEAF                              ***
C***                          READ_NEEDLE                            ***
C***                          READ_GROUND                            ***
C***                          READ_PRIMARY_BRANCH                    ***
C***                          READ_SECONDARY_BRANCH                  ***
C***                          READ_ENVIRONMENT                       ***
C***                          READ_NESTING                           ***
C***                          READ_EPS_TABLE                         *** 
C***                          READ_HISTOGRAM                         *** 
C*********************************************************************** 
C
        SUBROUTINE READ_INPUT
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C     THETA_DEGREES  =   radar look angle (degrees)
C     FREQ_GHZ       =   radar frequency  (GHz)
C----------------------------------------------------------------------
C-------------------- GROUND PARAMETERS -------------------------------
C----------------------------------------------------------------------
C     MV_SOIL      =  soil volumetric moisture
C     RMS_SOIL     =   surface RMS height  (centimeters)
C     LS_SOIL      =   surface correlation length (centimeters)
C     SAND         =   soil % sand
C     CLAY         =   soil % clay
C     SALT         =   standing water salt content (parts per thousand)
C----------------------------------------------------------------------
C-------------------- TRUNK PARAMETERS --------------------------------
C----------------------------------------------------------------------
C     MG_TRUNK       =   gravimetric trunk moisture
C     DENSITY        =   canopy tree density (trees per square meter)
C     RHO_TRUNK      =   dry density of trunk material
C     CROWN_HGHT     =   crown layer height (meters)
C     TRUNK_DIAM     =   trunk diameter  (centimeters)
C     TRUNK_HGHT     =   trunk layer height (meters)
C----------------------------------------------------------------------
C-------------------- LEAF PARAMETERS ---------------------------------
C----------------------------------------------------------------------
C     MG_LEAF        =   gravimetric leaf moisture
C     RHO_LEAF       =   dry density of leaf material
C     LEAF_DENS      =   leaf number density  (leaves per cubic meter)
C     LEAF_DIAM      =   leaf diameter  (centimeters)
C     LEAF_TAU       =   leaf thickness (centimeters)
C----------------------------------------------------------------------
C-------------------- NEEDLE PARAMETERS -------------------------------
C----------------------------------------------------------------------
C     MG_NDL         =  gravimetric needle moisture
C     RHO_NDL        =  dry density of needle material
C     NDL_DENS       =  needle number density (needles per cubic meter)
C     NDL_DIAM       =  needle diameter (centimeters)
C     NDL_LNG        =  needle length   (centimeters)
C----------------------------------------------------------------------
C------------------ PRIMARY BRANCH PARAMETERS -------------------------
C----------------------------------------------------------------------
C     MG_BR1         =  gravimetric branch moisture
C     RHO_BR1        =  dry density of branch material
C     BR1_DENS       =  branch number density (branches per cubic meter)
C     BR1_DIAM       =  branch diameter (centimeters)
C     BR1_LNG        =  branch length   (meters) 
C----------------------------------------------------------------------
C------------------ SECONDARY BRANCH PARAMETERS -----------------------
C----------------------------------------------------------------------
C     MG_BR2         =  gravimetric branch moisture
C     RHO_BR2        =  dry density of branch material
C     BR2_DENS       =  branch number density (branches per cubic meter)
C     BR2_DIAM       =  branch diameter (centimeters)
C     BR2_LNG        =  branch length   (meters) 
C----------------------------------------------------------------------
C---------------- ENVIRONMENT PARAMETERS ------------------------------
C----------------------------------------------------------------------
C     T_SOIL         =   temperature of soil  (degrees C)
C     T_VEG          =   temperature of vegetation (degrees C)
C     T_WATER        =   temperature of standing water (degrees C)
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER I
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
C
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
        LOGICAL LOG_MV_SOIL,LOG_ROUGH_SOIL(N_ROUGH_SOIL)
        LOGICAL LOG_SOIL_TEXT,LOG_SALT
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
        REAL BR1_DIAM_VECT(I_VECT_LEN)
        REAL BR1_LNG_VECT(I_VECT_LEN)
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
C-------- DECLARATIONS FROM SUBROUTINE READ_SECONDARY_BRANCH -----------
C
        INTEGER MG_BR2_NUM, RHO_BR2_NUM, BR2_DENS_NUM   
        INTEGER BR2_DIAM_NUM, BR2_LNG_NUM
        REAL MG_BR2_DELTA, RHO_BR2_DELTA, BR2_DENS_DELTA
        REAL BR2_DIAM_DELTA, BR2_LNG_DELTA
        REAL MG_BR2_VECT(I_VECT_LEN), RHO_BR2_VECT(I_VECT_LEN)
        REAL BR2_DENS_VECT(I_VECT_LEN)
        REAL BR2_DIAM_VECT(I_VECT_LEN)
        REAL BR2_LNG_VECT(I_VECT_LEN)
        REAL MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        LOGICAL LOG_MG_BR2, LOG_RHO_BR2, LOG_BR2_DENS
        LOGICAL LOG_BR2_DIAM, LOG_BR2_LNG
C
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        COMMON /R_BR2_DELTA/ MG_BR2_DELTA,RHO_BR2_DELTA,BR2_DENS_DELTA,
     &                        BR2_DIAM_DELTA, BR2_LNG_DELTA
        COMMON /I_BR2/ MG_BR2_NUM, RHO_BR2_NUM, BR2_DENS_NUM,
     &                  BR2_DIAM_NUM, BR2_LNG_NUM
        COMMON /R_BR2_VECT/MG_BR2_VECT, RHO_BR2_VECT, BR2_DENS_VECT,
     &                      BR2_DIAM_VECT, BR2_LNG_VECT
        COMMON /L_BR2/ LOG_MG_BR2, LOG_RHO_BR2, LOG_BR2_DENS,
     &                      LOG_BR2_DIAM, LOG_BR2_LNG
C
C-------- DECLARATIONS FROM SUBROUTINE READ_EPS_TABLE ------------------
C
        COMPLEX EPS_TR_VECT(I_VECT_LEN) ,EPS_TR_VECT_C(I_VECT_LEN)
        COMPLEX EPS_LF_VECT(I_VECT_LEN) ,EPS_LF_VECT_C(I_VECT_LEN)
        COMPLEX EPS_BR_1_VECT(I_VECT_LEN) ,EPS_BR_1_VECT_C(I_VECT_LEN)
        COMPLEX EPS_BR_2_VECT(I_VECT_LEN) ,EPS_BR_2_VECT_C(I_VECT_LEN)
        COMPLEX EPS_NDL_VECT(I_VECT_LEN) ,EPS_NDL_VECT_C(I_VECT_LEN)
        COMPLEX EPS_GRND_VECT(I_VECT_LEN) ,EPS_GRND_VECT_C(I_VECT_LEN)
        COMPLEX EPS_SNOW_VECT(I_VECT_LEN) ,EPS_SNOW_VECT_C(I_VECT_LEN) 
C
        COMMON /C_EPS_VECT/ EPS_TR_VECT, EPS_LF_VECT, EPS_BR_1_VECT,
     &       EPS_BR_2_VECT, EPS_NDL_VECT, EPS_GRND_VECT, EPS_SNOW_VECT
        COMMON /C_EPS_VECT_C/ EPS_TR_VECT_C, EPS_LF_VECT_C,
     &              EPS_BR_1_VECT_C, EPS_BR_2_VECT_C, EPS_NDL_VECT_C,
     &              EPS_GRND_VECT_C, EPS_SNOW_VECT_C
C
C-------- DECLARATIONS FROM SUBROUTINE READ_NESTING --------------------
C
        INTEGER PARAM_NUM(N_VARIABLES), INEST(N_VARIABLES)
C
        COMMON /I_NEST/ PARAM_NUM, INEST
C
C-------- DECLARATIONS FROM SUBROUTINE READ_TABLE ----------------------
C
        LOGICAL TABLE_OPEN
        COMMON /L_TABLE/ TABLE_OPEN
C
C***********************************************************************
C
        TABLE_OPEN = .FALSE.
C
        CALL READ_CONFIGURE
C
        CALL READ_SENSOR
C
        CALL READ_TRUNK
C
        IF(LOG_CONSTITUENT(4)) CALL READ_LEAF
C
        IF(LOG_CONSTITUENT(5)) CALL READ_NEEDLE
C
        IF(LOG_CONSTITUENT(2)) CALL READ_PRIMARY_BRANCH
        IF(LOG_CONSTITUENT(3)) CALL READ_SECONDARY_BRANCH
C
        CALL READ_GROUND
C
        CALL READ_ENVIRONMENT
C
        DO 100 I=1,N_CONSTITUENTS
            IF(LOG_EPS_TABLE(I)) CALL READ_EPS_TABLE(I)
100     CONTINUE
C
        IF(SNOW_SURFACE) CALL READ_EPS_TABLE(N_CONSTITUENTS+1)
C
        IF(LOG_HIST(1)) CALL READ_HISTOGRAM(1)
C
        CALL READ_NESTING
C
        IF(TABLE_OPEN)THEN
            CLOSE(3)
        ENDIF
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to Read Configuation.input   ***************
C***    Calling routine:      main MIMICS program                    ***
C***    Called subroutines:   SET_FLAG                               ***
C***                          SET_EPS_FLAG                           ***
C*********************************************************************** 
C
        SUBROUTINE READ_CONFIGURE
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C-----------------------------------------------------------------------
C LOGICAL VARIABLES:  .TRUE. = CREAT OUTPUT FILE
C                    .FLASE. = DO NOT CREAT OUTPUT FILE
C   L_POL = polarimetric.out
C   L_PHASE = forest_phase_mats.out
C   L_KAPPA = forest_kappa_mats.out
C   L_TRANS = forest_trans_mats.out
C-----------------------------------------------------------------------
C LOGICAL VARIABLES:  .TRUE. = CONSTITUENT PRESENT
C                    .FLASE. = CONSTITUENT NOT PRESENT
C   LOG_CONSTITUENT(1)  = trunks
C   LOG_CONSTITUENT(2)  = primary branches
C   LOG_CONSTITUENT(3)  = secondary branches
C   LOG_CONSTITUENT(4)  = leaves
C   LOG_CONSTITUENT(5)  = needles
C-----------------------------------------------------------------------
C LOGICAL VARIABLES:
C        .TRUE. = COMPUTE CONSTITUENT DIELECTRIC FROM TABLE VALUES
C       .FLASE. = COMPUTE CONSTITUENT DIELECTRIC FROM MOISTURE MODELS
C   LOG_EPS_TABLE(1)    = trunks
C   LOG_EPS_TABLE(2)    = primary branches
C   LOG_EPS_TABLE(3)    = secondary branches
C   LOG_EPS_TABLE(4)    = leaves
C   LOG_EPS_TABLE(5)    = needles
C-----------------------------------------------------------------------
C LOGICAL VARIABLES:
C       .TRUE. = COMPUTE DIELECTRIC USING A GIVEN DRY DENSITY
C      .FLASE. = COMPUTE DIELECTRIC ASSUMING DENSITY OF LEAFY VEGETATION
C   LOG_DRY_DENSITY(1)  = trunks
C   LOG_DRY_DENSITY(2)  = primary branches
C   LOG_DRY_DENSITY(3)  = secondary branches
C   LOG_DRY_DENSITY(4)  = leaves
C   LOG_DRY_DENSITY(5)  = needles
C-----------------------------------------------------------------------
C LOGICAL VARIABLES:
C   .TRUE. = P.D.F. COMPUTED USING SPECIFIED FUNCTION IN PDF SUBROUTINE
C  .FALSE. = P.D.F. COMPUTED USING DEFAULT FUNCTION READ IN INPUT FILES
C--- variation in size ---
C   LOG_PDF_TYPE(1,1)   = trunks
C   LOG_PDF_TYPE(2,1)   = primary branches
C   LOG_PDF_TYPE(3,1)   = secondary branches
C   LOG_PDF_TYPE(4,1)   = leaves
C   LOG_PDF_TYPE(5,1)   = needles
C
C LOGICAL VARIABLES:
C   .TRUE. = P.D.F. COMPUTED USING HISTOGRAM FILE
C   LOG_HIST(1)   = trunks
C   LOG_HIST(2)   = primary branches
C   LOG_HIST(3)   = secondary branches
C   LOG_HIST(4)   = leaves
C   LOG_HIST(5)   = needles
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER I, J
C
        INTEGER T, B_1, B_2, L, N
        INTEGER T_FORM, B_1_FORM, B_2_FORM, L_FORM, N_FORM, G_FORM
        INTEGER T_PDF_S, B_1_PDF_S, B_2_PDF_S, L_PDF_S, N_PDF_S
        integer pol, phase, KAPPA, TRANS
C
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
        LOGICAL L_POL, L_PHASE, L_KAPPA, L_TRANS
C
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                 LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
        COMMON /L_FILES/ L_POL, L_PHASE, L_KAPPA, L_TRANS
C
C***********************************************************************
C
        OPEN(UNIT=2,FILE='configuration.input',STATUS='UNKNOWN') 
C
        READ(2,5)
        READ(2,6)  POL, PHASE, KAPPA, TRANS
        READ(2,10)
        READ(2,20) T, B_1, B_2, L, N
        READ(2,30)
        READ(2,40) T_FORM, B_1_FORM, B_2_FORM, L_FORM, N_FORM, G_FORM
        READ(2,50)
        READ(2,40) T_PDF_S, B_1_PDF_S, B_2_PDF_S, L_PDF_S, N_PDF_S
C
        CLOSE(2)
C
C***********************************************************************
C------------------ SET LOGICAL VARIABLES ------------------------------
C***********************************************************************
C
        CALL SET_FLAG(POL,L_POL)
        CALL SET_FLAG(PHASE,L_PHASE)
        CALL SET_FLAG(KAPPA,L_KAPPA)
        CALL SET_FLAG(TRANS,L_TRANS)
C
        CALL SET_FLAG(T,LOG_CONSTITUENT(1))
        CALL SET_FLAG(B_1,LOG_CONSTITUENT(2))
        CALL SET_FLAG(B_2,LOG_CONSTITUENT(3))
        CALL SET_FLAG(L,LOG_CONSTITUENT(4))
        CALL SET_FLAG(N,LOG_CONSTITUENT(5))
C
        CALL SET_TWO_FLAGS(T_FORM,LOG_EPS_TABLE(1),LOG_DRY_DENSITY(1))
        CALL SET_TWO_FLAGS(B_1_FORM,LOG_EPS_TABLE(2),LOG_DRY_DENSITY(2))
        CALL SET_TWO_FLAGS(B_2_FORM,LOG_EPS_TABLE(3),LOG_DRY_DENSITY(3))
        CALL SET_TWO_FLAGS(L_FORM,LOG_EPS_TABLE(4),LOG_DRY_DENSITY(4))
        CALL SET_TWO_FLAGS(N_FORM,LOG_EPS_TABLE(5),LOG_DRY_DENSITY(5))
        CALL SET_TWO_FLAGS(G_FORM,LOG_EPS_TABLE(6),LOG_DRY_DENSITY(6))
C
C        CALL SET_FLAG(T_PDF_S, LOG_PDF_TYPE(1,1))
C        CALL SET_FLAG(B_1_PDF_S, LOG_PDF_TYPE(2,1))
C        CALL SET_FLAG(B_2_PDF_S, LOG_PDF_TYPE(3,1))
C        CALL SET_FLAG(L_PDF_S, LOG_PDF_TYPE(4,1))
C        CALL SET_FLAG(N_PDF_S, LOG_PDF_TYPE(5,1))
C
        CALL SET_TWO_FLAGS(T_PDF_S, LOG_HIST(1), LOG_PDF_TYPE(1,1))
        CALL SET_TWO_FLAGS(B_1_PDF_S, LOG_HIST(2), LOG_PDF_TYPE(2,1))
        CALL SET_TWO_FLAGS(B_2_PDF_S, LOG_HIST(3), LOG_PDF_TYPE(3,1))
        CALL SET_TWO_FLAGS(L_PDF_S, LOG_HIST(4), LOG_PDF_TYPE(4,1))
        CALL SET_TWO_FLAGS(N_PDF_S, LOG_HIST(5), LOG_PDF_TYPE(5,1))
C
C***********************************************************************
C---------------- CHECK FOR ALLOWED INPUT TYPES ------------------------
C***********************************************************************
C
        DO 101 J=1,1
          DO 100 I=2,5
            IF(LOG_PDF_TYPE(I,J).OR.LOG_HIST(I))THEN
                CALL WRITE_ERROR(3)
                STOP
            ENDIF
100       CONTINUE
101     CONTINUE
C
C***********************************************************************
C------------------ FORMAT STATEMENTS ----------------------------------
C***********************************************************************
C
5      FORMAT(11(/))
6      FORMAT(12X,I1,12X,I1,12X,I1,13x,I1)
10      FORMAT(10(/))
20      FORMAT(10X,I1,15X,I1,12X,I1,14X,I1,14X,I1)
30      FORMAT(11(/))
40      FORMAT(13X,I1,11X,I1,12X,I1,15X,I1,11X,I1,11X,I1)
50      FORMAT(10(/))
60      FORMAT(13X,I1,11X,I1,12X,I1,15X,I1,11X,I1)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to set two logical variables  **************
C***    Calling routine:      READ_CONFUGIRE                         ***
C***    Called subroutines:   WRITE_ERROR                            ***
C*********************************************************************** 
C
        SUBROUTINE SET_TWO_FLAGS(I, FLAG_1, FLAG_2)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   I      = 0,1 or 2 (integer)
C   FLAG_1 = logical variable
C          = .TRUE.  if I = 2
C          = .FALSE. if I = 0 or 1
C   FLAG_2 = logical variable
C          = .TRUE.  if I = 1
C          = .FALSE. if I = 0 or 2
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER I
        LOGICAL FLAG_1, FLAG_2
C
C***********************************************************************
C
        IF(I.EQ.0)THEN
          FLAG_1 = .FALSE.
          FLAG_2 = .FALSE.
        ELSE IF(I.EQ.1)THEN
          FLAG_1 = .FALSE.
          FLAG_2 = .TRUE.
        ELSE IF(I.EQ.2)THEN
          FLAG_1 = .TRUE.
          FLAG_2 = .FALSE.
        ELSE
          CALL WRITE_ERROR(10)
          STOP
        ENDIF
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to set a logical variable     **************
C***    Calling routine:      READ_CONFIGURE                         ***
C***                          READ_SENSOR                            ***
C***                          READ_ENVIRONMENT                       ***
C***                          READ_TRUNK                             ***
C***                          READ_GROUND                            ***
C***                          READ_LEAF                              ***
C***                          READ_NEEDLE                            ***
C***                          READ_PRIMARY_BRANCH                    ***
C***                          READ_SECONDARY_BRANCH                  ***
C***    Called subroutines:   WRITE_ERROR                            ***
C*********************************************************************** 
C
        SUBROUTINE SET_FLAG(I,FLAG)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   I    = 0 or 1 (integer)
C   FLAG = logical variable
C        = .TRUE.  if I = 1
C        = .FALSE. if I = 0
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER I
        LOGICAL FLAG
C
C***********************************************************************
C
        IF(I.EQ.0)THEN
            FLAG = .FALSE.
        ELSE IF(I.EQ.1)THEN
            FLAG = .TRUE.
        ELSE
            CALL WRITE_ERROR(15)
            STOP
        ENDIF
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C*** Function to compute number of loops through the supplied input  ***
C*** data based on the starting, stopping, and incremental value of  ***
C*** that parameter.                                                 ***
C***    Calling routine:      READ_SENSOR                            ***
C***    Called subroutines:   WRITE_ERROR                            ***
C*********************************************************************** 
C
        INTEGER FUNCTION LOOP(START,STOP,DELTA)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   START   =   starting value
C   STOP    =   stopping value
C   DELTA   =   increment value
C   LOOP    =   number of loops through the input data
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
        REAL START,STOP,DELTA
C
C***********************************************************************
C----- compute number of increments ------------------------------------
C***********************************************************************
C
        IF(DELTA.NE.0)THEN
            LOOP = INT((STOP-START)/DELTA) + 1
            IF ((START+(LOOP-1)*DELTA).LT.(STOP-0.1*DELTA))THEN
                LOOP = LOOP + 1
            ENDIF
        ELSE IF(START.NE.STOP)THEN
            CALL WRITE_ERROR(20)
            STOP
        ELSE
            LOOP = 1
        ENDIF
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read SENSOR.INPUT          **************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_SENSOR
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   FREQ_START  =   starting value of frequency (GHz)
C   FREQ_STOP   =   stopping value of frequency (GHz)
C   FREQ_DELTA  =   increment for stepping frequency (GHz)
C   THETA_START =   starting value of radar look angle (degrees)
C   THETA_STOP  =   stopping value of radar look angle (degrees)
C   THETA_DELTA =   increment for stepping radar look angle (degrees)
C
C   FREQ_NUMBER     = number of frequency steps
C   THETA_NUMBER    = number of radar look angle steps
C
C   LOG_FREQ_TABLE  = logical variable indicating look-up table data
C   LOG_THETA_TABLE = logical variable indicating look-up table data
C                   = .TRUE. use PARAMETER_VALUE_TABLE.INPUT
C                   = .FALSE. do not use PARAMETER_VALUE_TABLE.INPUT
C
C   THETA_VECT  =   VECTOR CONTAINING LOOK-UP TABLE VALUES FOR THETA
C   FREQ_VECT   =   VECTOR CONTAINING LOOK-UP TABLE VALUES FOR FREQUENCY
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER FREQ_NUMBER, FREQ_T, THETA_NUMBER, THETA_T
        INTEGER LOOP
C
        REAL FREQ_START, FREQ_STOP, FREQ_DELTA
        REAL FREQ_VECT(I_VECT_LEN)
        REAL THETA_START, THETA_STOP, THETA_DELTA
        REAL THETA_VECT(I_VECT_LEN)
        REAL THETA_DEGREES, FREQ_GHZ
C
        LOGICAL LOG_FREQ_TABLE, LOG_THETA_TABLE
C
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
        COMMON /R_SENSOR_DELTA/ THETA_DELTA, FREQ_DELTA
        COMMON /I_SENSOR/ THETA_NUMBER, FREQ_NUMBER
        COMMON /R_SENSOR_VECT/ THETA_VECT, FREQ_VECT
        COMMON /L_SENSOR/ LOG_FREQ_TABLE,LOG_THETA_TABLE
        COMMON /SENSOR_START/ FREQ_START,THETA_START
C
C***********************************************************************
C   READ INPUT DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='sensor.input',STATUS='UNKNOWN') 
C
        READ(2,10)
        READ(2,20) FREQ_START, FREQ_STOP, FREQ_DELTA, FREQ_T
        READ(2,30)
        READ(2,20) THETA_START, THETA_STOP, THETA_DELTA, THETA_T
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(FREQ_T,LOG_FREQ_TABLE)
        CALL SET_FLAG(THETA_T,LOG_THETA_TABLE)
C
        IF(LOG_FREQ_TABLE)THEN
            CALL READ_TABLE(FREQ_VECT,FREQ_NUMBER,1)
C            print *,'1: FREQ_GHZ'
            FREQ_GHZ = FREQ_VECT(1)
            FREQ_START = FREQ_VECT(1)
        ELSE
            FREQ_NUMBER = LOOP(FREQ_START, FREQ_STOP, FREQ_DELTA)
            FREQ_GHZ = FREQ_START
        ENDIF
C
        IF(LOG_THETA_TABLE)THEN
            CALL READ_TABLE(THETA_VECT,THETA_NUMBER,2)
C            print *,'2: THETA_DEGREES'
            THETA_DEGREES = THETA_VECT(1)
            THETA_START = THETA_VECT(1)
        ELSE
            THETA_NUMBER = LOOP(THETA_START, THETA_STOP, THETA_DELTA)
            THETA_DEGREES = THETA_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(16x,f5.2,21x,f5.2,19x,f5.2,8X,I1)
30      FORMAT(5(/))
40      FORMAT(16x,f4.1,22x,f4.1,20x,f4.1,9X,I1)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read ENVIRONMENT.INPUT     **************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_ENVIRONMENT
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   T_SOIL_START   =   starting value of soil temperature (degrees C)
C   T_SOIL_STOP    =   stopping value of soil temperature (degrees C)
C   T_SOIL_DELTA   =   increment for stepping soil temperature (deg. C)
C
C   T_WATER_START  =   starting value of water temperature (degrees C)
C   T_WATER_STOP   =   stopping value of water temperature (degrees C)
C   T_WATER_DELTA  =   increment for stepping water temperature (deg. C)
C
C   T_VEG_START    =   starting value of vegetation temperature (deg C)
C   T_VEG_STOP     =   stopping value of vegetation temperature (deg C)
C   T_VEG_DELTA    =   increment for stepping vegetation temp. (deg C)
C
C   T_SOIL_NUMBER  = number of soil temperature steps
C   T_WATER_NUMBER = number of standing water temperature steps
C   T_VEG_NUMBER   = number of vegetation temperature steps
C
C   LOG_ENVIRONMENT() = logical variable indicating look-up table data
C                     = .TRUE. use PARAMETER_VALUE_TABLE.INPUT
C                     = .FALSE. do not use PARAMETER_VALUE_TABLE.INPUT
C   LOG_ENVIRONMENT(1) = data for soil temperature
C   LOG_ENVIRONMENT(2) = data for standing water temperature
C   LOG_ENVIRONMENT(3) = data for vegetation temperature
C
C   SOIL_T_VECT  = vector containing look-up table values for soil temp
C   WATER_T_VECT = vector containing look-up table values for water temp
C   VEG_T_VECT   = vector containing look-up table values for veg. temp
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_SOIL_NUMBER, T_WATER_NUMBER, T_VEG_NUMBER
        INTEGER T_SOIL_T, T_WATER_T, T_VEG_T
        INTEGER LOOP
C
        REAL T_SOIL_START, T_SOIL_STOP, T_SOIL_DELTA
        REAL T_WATER_START, T_WATER_STOP, T_WATER_DELTA
        REAL T_VEG_START, T_VEG_STOP, T_VEG_DELTA
        REAL T_SOIL_VECT(I_VECT_LEN), T_WATER_VECT(I_VECT_LEN)
        REAL T_VEG_VECT(I_VECT_LEN)
        REAL T_SOIL,T_VEG,T_WATER
C
        LOGICAL LOG_ENVIRONMENT(N_ENVIRONMENT)
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_ENVIRON_DELTA/ T_SOIL_DELTA,T_WATER_DELTA,T_VEG_DELTA
        COMMON /I_ENVIRON/ T_SOIL_NUMBER, T_WATER_NUMBER, T_VEG_NUMBER
        COMMON /R_ENVIRON_VECT/ T_SOIL_VECT, T_WATER_VECT, T_VEG_VECT
        COMMON /L_ENVIRON/ LOG_ENVIRONMENT
        COMMON /ENV_START/ T_SOIL_START,T_WATER_START,T_VEG_START
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='environment.input',STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) T_SOIL_START, T_SOIL_STOP, T_SOIL_DELTA, T_SOIL_T
        READ(2,30)
        READ(2,20) T_WATER_START, T_WATER_STOP, T_WATER_DELTA, T_WATER_T
        READ(2,30)
        READ(2,20) T_VEG_START, T_VEG_STOP, T_VEG_DELTA, T_VEG_T
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_SOIL_T,LOG_ENVIRONMENT(1))
        CALL SET_FLAG(T_WATER_T,LOG_ENVIRONMENT(2))
        CALL SET_FLAG(T_VEG_T,LOG_ENVIRONMENT(3))
C
        IF(LOG_ENVIRONMENT(1))THEN
            CALL READ_TABLE(T_SOIL_VECT,T_SOIL_NUMBER,3)
C            print *,'3: T_SOIL'
            T_SOIL = T_SOIL_VECT(1)
            T_SOIL_START = T_SOIL_VECT(1)
        ELSE
            T_SOIL_NUMBER = LOOP(T_SOIL_START,T_SOIL_STOP,T_SOIL_DELTA)
            T_SOIL = T_SOIL_START
        ENDIF
C
        IF(LOG_ENVIRONMENT(2))THEN
         CALL READ_TABLE(T_WATER_VECT,T_WATER_NUMBER,4)
C         print *,'4: T_WATER'
         T_WATER = T_WATER_VECT(1)
         T_WATER_START = T_WATER_VECT(1)
        ELSE
         T_WATER_NUMBER = LOOP(T_WATER_START,T_WATER_STOP,T_WATER_DELTA)
         T_WATER = T_WATER_START
        ENDIF
C
        IF(LOG_ENVIRONMENT(3))THEN
            CALL READ_TABLE(T_VEG_VECT,T_VEG_NUMBER,5)
C            print *,'5: T_VEG'
            T_VEG = T_VEG_VECT(1)
            T_VEG_START = T_VEG_VECT(1)
        ELSE
            T_VEG_NUMBER = LOOP(T_VEG_START, T_VEG_STOP, T_VEG_DELTA)
            T_VEG = T_VEG_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(16X,F5.1,21X,F5.1,19X,F5.1,8X,I1)
30      FORMAT(5(/))
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read GROUND.INPUT          **************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   WRITE_ERROR                            ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_GROUND
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   SURFACE_TYPE   =   indicates type of surface
C                      1 = SOIL
C                      2 = STANDING WATER
C                      3 = ICE
C   SNOW_FLAG      = INDICATES PRESENCE OR ABSENCE OF A SNOW LAYER
C                      0 = NO SNOW
C                      1 = SNOW PRESENT
C---LOGICAL VARIABLES---
C                   = .TRUE.   -- analize for this type of surface
C                   = .FALSE.  -- surface is not of this type
C     SOIL_SURFACE  
C     WATER_SURFACE
C     SNOW_SURFACE
C     ICE_SURFACE
C-----------------------
C
C   MV_SOIL_START  =   starting value of soil volumetric moisture
C   MV_SOIL_STOP   =   stopping value of soil volumetric moisture
C   MV_SOIL_DELTA  =   increment for stepping soil volumetric moisture
C
C   RMS_SOIL_START =   starting value of soil RMS height (centimeters)
C   RMS_SOIL_STOP  =   stopping value of soil RMS height (centimeters)
C   RMS_SOIL_DELTA =   increment for stepping soil RMS height (cm)
C
C   LS_SOIL_START  =   starting value of soil correlation length (cm)
C   LS_SOIL_STOP   =   stopping value of soil correlation length (cm)
C   LS_SOIL_DELTA  =   increment for stepping soil corr. length (cm)
C
C   SAND_START     =   starting value of soil sand percent by weight
C   SAND_STOP      =   stopping value of soil sand percent by weight
C   SAND_DELTA     =   increment for stepping soil sand % by weight
C
C   CLAY_START     =   starting value of soil clay percent by weight
C   CLAY_STOP      =   stopping value of soil clay percent by weight
C   CLAY_DELTA     =   increment for stepping soil clay % by weight
C
C   SALT_START     =   starting value of standing water salt (ppt)
C   SALT_STOP      =   stopping value of standing water salt (ppt)
C   SALT_DELTA     =   increment for stepping standing water salt (ppt)
C
C   MV_SOIL_NUMBER  = number of soil moisture steps
C   RMS_SOIL_NUMBER = number of soil RMS height steps
C   LS_SOIL_NUMBER  = number of soil correlation length steps
C   SAND_NUMBER     = number of soil % sand steps
C   CLAY_NUMBER     = number of soil % clay steps
C   SALT_NUMBER     = number of standing water salt ppt steps
C
C--- logical variable indicating look-up table data --------------------
C                  = .TRUE. use PARAMETER_VALUE_TABLE.INPUT
C                  = .FALSE. do not use PARAMETER_VALUE_TABLE.INPUT
C   LOG_MV_SOIL       = data for soil moisture
C   LOG_ROUGH_SOIL(1) = data for RMS height
C   LOG_ROUGH_SOIL(2) = data for correlation length
C   LOG_SOIL_TEXT  = data for soil percent sand, silt, clay
C   LOG_SALT       = data for standing water salt
C
C-------------- look-up table vectors ----------------------------------
C
C   MV_SOIL_VECT  = vector containing values for soil moisture
C   RMS_SOIL_VECT = vector containing values for soil RMS heoght
C   LS_SOIL_VECT  = vector containing values for soil correlation length
C   SAND_VECT     = vector containing values for soil % sand
C   CLAY_VECT     = vector containing values for soil % clay
C   SALT_VECT     = vector containing values for staning water salt
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER MV_SOIL_NUMBER, RMS_SOIL_NUMBER, LS_SOIL_NUMBER
        INTEGER SAND_NUMBER, CLAY_NUMBER, SALT_NUMBER
        INTEGER T_SNOW_NUMBER
        INTEGER T_RMS_SOIL, T_MV_SOIL, T_LS_SOIL,T_SAND, T_CLAY, T_SALT
        INTEGER LOOP, SURFACE_TYPE, SNOW_FLAG
        INTEGER SNOW_T_TABLE, i_surf
C
        REAL MV_SOIL_START, MV_SOIL_STOP, MV_SOIL_DELTA
        REAL RMS_SOIL_START, RMS_SOIL_STOP, RMS_SOIL_DELTA
        REAL LS_SOIL_START, LS_SOIL_STOP, LS_SOIL_DELTA
        REAL SAND_START, SAND_STOP, SAND_DELTA
        REAL CLAY_START, CLAY_STOP, CLAY_DELTA
        REAL SALT_START, SALT_STOP, SALT_DELTA
        REAL MV_SOIL_VECT(I_VECT_LEN)
        REAL RMS_SOIL_VECT(I_VECT_LEN)
        REAL LS_SOIL_VECT(I_VECT_LEN)
        REAL SAND_VECT(I_VECT_LEN), CLAY_VECT(I_VECT_LEN)
        REAL SALT_VECT(I_VECT_LEN)
        REAL SNOW_T_VECT(I_VECT_LEN)
        REAL MV_SOIL, SAND, CLAY, RMS_SOIL, LS_SOIL
        REAL SALT, T_SNOW
        REAL SNOW_T_START, SNOW_T_STOP, SNOW_T_DELTA
C
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
        LOGICAL LOG_MV_SOIL, LOG_ROUGH_SOIL(N_ROUGH_SOIL)
        LOGICAL LOG_SOIL_TEXT, LOG_SALT,  LOG_SNOW_T
C
C
        COMMON /R_SURFACE/ RMS_SOIL, LS_SOIL
        COMMON /R_SOIL/ MV_SOIL, SAND, CLAY
        COMMON /R_WATER/ SALT
        COMMON /R_SNOW/ T_SNOW
C
        COMMON /R_SURFACE_DELTA/ RMS_SOIL_DELTA, LS_SOIL_DELTA
        COMMON /R_SOIL_DELTA/ MV_SOIL_DELTA, SAND_DELTA, CLAY_DELTA
        COMMON /R_WATER_DELTA/ SALT_DELTA
        COMMON /R_SNOW_DELTA/ SNOW_T_DELTA
C
        COMMON /I_SURFACE/ RMS_SOIL_NUMBER, LS_SOIL_NUMBER
        COMMON /I_SOIL/ MV_SOIL_NUMBER, SAND_NUMBER, CLAY_NUMBER
        COMMON /I_WATER/ SALT_NUMBER
        COMMON /I_SNOW/ T_SNOW_NUMBER
C
        COMMON /R_SURFACE_VECT/ RMS_SOIL_VECT, LS_SOIL_VECT
        COMMON /R_SOIL_VECT/ MV_SOIL_VECT,SAND_VECT,CLAY_VECT
        COMMON /R_WATER_VECT/ SALT_VECT
        COMMON /R_SNOW_VECT/ SNOW_T_VECT
C
        COMMON /L_SURFACE/ LOG_ROUGH_SOIL
        COMMON /L_SOIL/ LOG_MV_SOIL, LOG_SOIL_TEXT
        COMMON /L_WATER/ LOG_SALT
        COMMON /L_SNOW/ LOG_SNOW_T
C
        COMMON /L_SURFACE_TYPE/ SOIL_SURFACE,WATER_SURFACE,
     &                          SNOW_SURFACE,ICE_SURFACE
        common /surface_mod/ i_surf
C
        COMMON /SURFACE_START/ RMS_SOIL_START, LS_SOIL_START
        COMMON /SOIL_START/ MV_SOIL_START ,SAND_START , CLAY_START
        COMMON /WATER_START/ SALT_START
        COMMON /SNOW_START/ SNOW_T_START
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='ground.input',STATUS='UNKNOWN')
C
       READ(2,10)
       READ(2,20) SURFACE_TYPE, SNOW_FLAG
       READ(2,30)
       READ(2,40) MV_SOIL_START, MV_SOIL_STOP, MV_SOIL_DELTA, T_MV_SOIL
       READ(2,50)
       READ(2,60) RMS_SOIL_START,RMS_SOIL_STOP,RMS_SOIL_DELTA,T_RMS_SOIL
       READ(2,50)
       READ(2,70) LS_SOIL_START, LS_SOIL_STOP, LS_SOIL_DELTA, T_LS_SOIL
       READ(2,80) 
       READ(2,90) SAND_START, SAND_STOP, SAND_DELTA ,T_SAND
       READ(2,120) 
       READ(2,90) CLAY_START, CLAY_STOP, CLAY_DELTA, T_CLAY
       READ(2,100) 
       READ(2,110) SALT_START, SALT_STOP, SALT_DELTA, T_SALT
       READ(2,30)
       READ(2,40) SNOW_T_START, SNOW_T_STOP, SNOW_T_DELTA, SNOW_T_TABLE
       read(2,80) 
       read(2,130) i_surf
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        IF(SURFACE_TYPE.EQ.1)THEN
            SOIL_SURFACE = .TRUE. 
            WATER_SURFACE = .FALSE.
            ICE_SURFACE = .FALSE.
        ELSE IF(SURFACE_TYPE.EQ.2)THEN
            SOIL_SURFACE = .FALSE. 
            WATER_SURFACE = .TRUE.
            ICE_SURFACE = .FALSE.
        ELSE IF(SURFACE_TYPE.EQ.3)THEN
            SOIL_SURFACE = .FALSE. 
            WATER_SURFACE = .FALSE.
            ICE_SURFACE = .TRUE.
        ELSE
            CALL WRITE_ERROR(25)
            STOP
        ENDIF
C
        CALL SET_FLAG(SNOW_FLAG, SNOW_SURFACE)
C
        IF((T_SAND.EQ.1).OR.(T_CLAY.EQ.1))THEN
            LOG_SOIL_TEXT = .TRUE.
        ENDIF
C
        CALL SET_FLAG(T_MV_SOIL, LOG_MV_SOIL)
        CALL SET_FLAG(T_RMS_SOIL, LOG_ROUGH_SOIL(1))
        CALL SET_FLAG(T_LS_SOIL,  LOG_ROUGH_SOIL(2))
        CALL SET_FLAG(T_SALT, LOG_SALT)
        CALL SET_FLAG(SNOW_T_TABLE, LOG_SNOW_T)
C
        IF(LOG_MV_SOIL)THEN
            CALL READ_TABLE(MV_SOIL_VECT(1), MV_SOIL_NUMBER,6)
C            print *,'6: MV_SOIL'  
            MV_SOIL = MV_SOIL_VECT(1)
            MV_SOIL_START = MV_SOIL_VECT(1)
        ELSE
            MV_SOIL_NUMBER = LOOP(MV_SOIL_START, MV_SOIL_STOP,
     &                            MV_SOIL_DELTA)
            MV_SOIL = MV_SOIL_START
        ENDIF
C
        IF(LOG_ROUGH_SOIL(1))THEN
            CALL READ_TABLE(RMS_SOIL_VECT(1), RMS_SOIL_NUMBER,7)
C            print *,'7: RMS_SOIL'
            RMS_SOIL = RMS_SOIL_VECT(1)
            RMS_SOIL_START = RMS_SOIL_VECT(1)
        ELSE
            RMS_SOIL_NUMBER = LOOP(RMS_SOIL_START, RMS_SOIL_STOP,
     &                             RMS_SOIL_DELTA)
            RMS_SOIL = RMS_SOIL_START
        ENDIF
C
        IF(LOG_ROUGH_SOIL(2))THEN
            CALL READ_TABLE(LS_SOIL_VECT (1), LS_SOIL_NUMBER,8)
C            print *,'8: LS_SOIL'
            LS_SOIL = LS_SOIL_VECT(1)
            LS_SOIL_START = LS_SOIL_VECT(1)
        ELSE
            LS_SOIL_NUMBER = LOOP(LS_SOIL_START, LS_SOIL_STOP,
     &                            LS_SOIL_DELTA)
            LS_SOIL = LS_SOIL_START
        ENDIF
C
        IF(LOG_SOIL_TEXT)THEN
            PRINT*,' NOT YET SET UP TO READ SOIL TEXTURE VECTOR'
            STOP
C            CALL READ_TABLE(SAND_VECT(1),SAND_NUMBER)
C            SAND = SAND_VECT(1)
        ELSE
            SAND_NUMBER = LOOP(SAND_START, SAND_STOP, SAND_DELTA)
            SAND = SAND_START
            CLAY_NUMBER = LOOP(CLAY_START, CLAY_STOP, CLAY_DELTA)
            CLAY = CLAY_START
        ENDIF
C
        IF(LOG_SALT)THEN
            CALL READ_TABLE(SALT_VECT(1),SALT_NUMBER,9)
C            print *,'9: SALT'
            SALT = SALT_VECT(1)
            SALT_START = SALT_VECT(1)
        ELSE
            SALT_NUMBER = LOOP(SALT_START, SALT_STOP, SALT_DELTA)
            SALT = SALT_START
        ENDIF
C
        IF(SNOW_SURFACE)THEN
          IF(LOG_SNOW_T)THEN
            CALL READ_TABLE(SNOW_T_VECT(1),T_SNOW_NUMBER,10)
C            print *,'10: T_SNOW'
            T_SNOW = SNOW_T_VECT(1)
            SNOW_T_START = SNOW_T_VECT(1)
          ELSE
            T_SNOW_NUMBER = LOOP(SNOW_T_START,SNOW_T_STOP,SNOW_T_DELTA)
            T_SNOW = SNOW_T_START
          ENDIF
        ELSE
            T_SNOW_NUMBER = 1
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(24X,I1,22X,I1)
30      FORMAT(7(/))
40      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
50      FORMAT(5(/))
60      FORMAT(17X,F7.4,18X,F7.4,18X,F7.4,6X,I2)
70      FORMAT(17X,F8.3,17X,F8.3,17X,F8.3,4X,I2)
80      FORMAT(7(/))
90      FORMAT(17X,F6.2,19X,F6.2,19X,F6.2,6X,I1)
100     FORMAT(7(/))
110     FORMAT(17X,F5.2,20X,F5.2,20X,F5.2,7X,I1)
120     FORMAT(3(/))
130     format(17x,i2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read TRUNK_AND_GROSS_CANOPY.INPUT *******
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_TRUNK
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   MG_TRUNK_START =   starting value of trunk gravimetric moisture
C   MG_TRUNK_STOP  =   stopping value of trunk gravimetric moisture
C   MG_TRUNK_DELTA =   increment for stepping trunk grav. moisture
C
C   RHO_TRUNK_START =   starting value of trunk material dry density
C   RHO_TRUNK_STOP  =   stopping value of trunk material dry density
C   RHO_TRUNK_DELTA =   increment for stepping material dry density
C
C   DENSITY_START =   starting value of canopy density (trees/meter**2)
C   DENSITY_STOP  =   stopping value of canopy density (trees/meter**2)
C   DENSITY_DELTA =   increment for stepping canopy density (trees/m**2)
C
C   CROWN_HGHT_START =   starting value of crown height (meters)
C   CROWN_HGHT_STOP  =   stopping value of crown height (meters)
C   CROWN_HGHT_DELTA =   increment for stepping crown height (meters)
C
C   TRUNK_DIAM_START =   starting value of trunk diameter (centimeters)
C   TRUNK_DIAM_STOP  =   stopping value of trunk diameter (centimeters)
C   TRUNK_DIAM_DELTA =   increment for stepping diameter (centimeters)
C
C   TRUNK_HGHT_START =   starting value of trunk height (meters)
C   TRUNK_HGHT_STOP  =   stopping value of trunk height (meters)
C   TRUNK_HGHT_DELTA =   increment for stepping height (meters)
C
C   TRUNK_Dsig_START =   starting value of trunk diameter sigma (cm)
C   TRUNK_Dsig_STOP  =   stopping value of trunk diameter sigma (cm)
C   TRUNK_Dsig_DELTA =   increment for stepping diameter sigma (cm)
C
C   TRUNK_Hsig_START =   starting value of trunk height sigma (meters)
C   TRUNK_Hsig_STOP  =   stopping value of trunk height sigma (meters)
C   TRUNK_Hsig_DELTA =   increment for stepping height sigma (meters)
C
C   MG_TRUNK_NUM    = number of trunk moisture steps
C   RHO_TRUNK_NUM   = number of trunk dry density steps
C   DENSITY_NUM     = number of tree density steps
C   CROWN_HGHT_NUM  = number of crown height steps
C   TRUNK_DIAM_NUM  = number of trunk diameter steps
C   TRUNK_HGHT_NUM  = number of trunk height steps
C   TRUNK_Dsig_NUM  = number of trunk diameter sigma (stand. dev.) steps
C   TRUNK_Hsig_NUM  = number of trunk height sigma (standard dev.) steps
C
C   LOG_MG_TRUNK    = logical variable indicating look-up table data
C   LOG_RHO_TRUNK   = logical variable indicating look-up table data
C   LOG_DENSITY     = logical variable indicating look-up table data
C   LOG_CROWN_HGHT  = logical variable indicating look-up table data
C   LOG_TRUNK_DIAM  = logical variable indicating look-up table data
C   LOG_TRUNK_HGHT  = logical variable indicating look-up table data
C   LOG_TRUNK_Dsig  = logical variable indicating look-up table data
C   LOG_TRUNK_Hsig  = logical variable indicating look-up table data
C
C   MG_TRUNK_VECT   = vector containing look-up table values
C   RHO_TRUNK_VECT  = vector containing look-up table values
C   DENSITY_VECT    = vector containing look-up table values
C   CROWN_HGHT_VECT = vector containing look-up table values
C   TRUNK_DIAM_VECT = vector containing look-up table values
C   TRUNK_HGHT_VECT = vector containing look-up table values
C   TRUNK_Dsig_VECT = vector containing look-up table values
C   TRUNK_Hsig_VECT = vector containing look-up table values
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_MG_TRUNK, T_RHO_TRUNK, T_DENSITY   
        INTEGER T_CROWN_HGHT, T_TRUNK_DIAM, T_TRUNK_HGHT
        INTEGER T_TRUNK_Dsig, T_TRUNK_Hsig
        INTEGER MG_TRUNK_NUM, RHO_TRUNK_NUM, DENSITY_NUM   
        INTEGER CROWN_HGHT_NUM, TRUNK_DIAM_NUM, TRUNK_HGHT_NUM
        INTEGER TRUNK_Dsig_NUM, TRUNK_Hsig_NUM
        INTEGER LOOP
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,I_PDF_NDL,I_PDF_TRUNK
        INTEGER I_PDF_LEAF_SIZE,I_PDF_BR_1_SIZE
        INTEGER I_PDF_BR_2_SIZE,I_PDF_NDL_SIZE,I_PDF_TRUNK_SIZE

C
        REAL MG_TRUNK_START, MG_TRUNK_STOP, MG_TRUNK_DELTA
        REAL RHO_TRUNK_START, RHO_TRUNK_STOP, RHO_TRUNK_DELTA
        REAL DENSITY_START, DENSITY_STOP, DENSITY_DELTA
        REAL CROWN_HGHT_START, CROWN_HGHT_STOP, CROWN_HGHT_DELTA
        REAL TRUNK_DIAM_START, TRUNK_DIAM_STOP, TRUNK_DIAM_DELTA
        REAL TRUNK_HGHT_START, TRUNK_HGHT_STOP, TRUNK_HGHT_DELTA
        REAL TRUNK_Dsig_START, TRUNK_Dsig_STOP, TRUNK_Dsig_DELTA
        REAL TRUNK_Hsig_START, TRUNK_Hsig_STOP, TRUNK_Hsig_DELTA
        REAL MG_TRUNK_VECT(I_VECT_LEN)
        REAL RHO_TRUNK_VECT(I_VECT_LEN)
        REAL DENSITY_VECT(I_VECT_LEN)
        REAL CROWN_HGHT_VECT(I_VECT_LEN)
        REAL TRUNK_DIAM_VECT(I_VECT_LEN)
        REAL TRUNK_HGHT_VECT(I_VECT_LEN)
        REAL TRUNK_Dsig_VECT(I_VECT_LEN)
        REAL TRUNK_Hsig_VECT(I_VECT_LEN)
        REAL MG_TRUNK, RHO_TRUNK, DENSITY, CROWN_HGHT
        REAL TRUNK_DIAM,TRUNK_HGHT
        REAL TRUNK_Dsig,TRUNK_Hsig
C
        LOGICAL LOG_MG_TRUNK, LOG_RHO_TRUNK, LOG_DENSITY
        LOGICAL LOG_CROWN_HGHT, LOG_TRUNK_DIAM, LOG_TRUNK_HGHT
        LOGICAL LOG_TRUNK_Dsig, LOG_TRUNK_Hsig
C
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_TRUNK_DELTA/ MG_TRUNK_DELTA, RHO_TRUNK_DELTA, 
     &                         TRUNK_DIAM_DELTA
        COMMON /I_TRUNK/ MG_TRUNK_NUM, RHO_TRUNK_NUM, TRUNK_DIAM_NUM
        COMMON /R_TRUNK_VECT/MG_TRUNK_VECT, RHO_TRUNK_VECT,
     &                         TRUNK_DIAM_VECT
        COMMON /L_TRUNK/ LOG_MG_TRUNK,LOG_RHO_TRUNK,LOG_TRUNK_DIAM
        COMMON /TRUNK_START/ MG_TRUNK_START, RHO_TRUNK_START,
     &                       TRUNK_DIAM_START 
C
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_CANOPY_DELTA/ DENSITY_DELTA, CROWN_HGHT_DELTA,
     &                          TRUNK_HGHT_DELTA
        COMMON /I_CANOPY/ DENSITY_NUM, CROWN_HGHT_NUM, TRUNK_HGHT_NUM
        COMMON /R_CANOPY_VECT/DENSITY_VECT, CROWN_HGHT_VECT,
     &                        TRUNK_HGHT_VECT
        COMMON /L_CANOPY/ LOG_DENSITY, LOG_CROWN_HGHT, LOG_TRUNK_HGHT
        COMMON /CANOPY_START/ DENSITY_START,CROWN_HGHT_START,
     &                        TRUNK_HGHT_START
C
        COMMON /I_TRUNK_SIG/ TRUNK_Dsig_NUM, TRUNK_Hsig_NUM
        COMMON /R_TRUNK_SIG/ TRUNK_Dsig, TRUNK_Hsig
        COMMON /TRUNK_SIG_DELTA/ TRUNK_Dsig_DELTA, TRUNK_Hsig_DELTA
        COMMON /TRUNK_SIG_VECT/ TRUNK_Dsig_VECT,  TRUNK_Hsig_VECT
        COMMON /TRUNK_SIG_START/ TRUNK_Dsig_START, TRUNK_Hsig_START
        COMMON /L_TRUNK_SIG/ LOG_TRUNK_Dsig, LOG_TRUNK_Hsig
C
        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_NDL,I_PDF_TRUNK
        COMMON /I_PDF_SIZE/ I_PDF_LEAF_SIZE,I_PDF_BR_1_SIZE,
     &               I_PDF_BR_2_SIZE,I_PDF_NDL_SIZE,I_PDF_TRUNK_SIZE
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='trunk_and_gross_canopy.input',
     &        STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) MG_TRUNK_START, MG_TRUNK_STOP, MG_TRUNK_DELTA,
     &             T_MG_TRUNK
        READ(2,30)
        READ(2,40) RHO_TRUNK_START, RHO_TRUNK_STOP, RHO_TRUNK_DELTA,
     &             T_RHO_TRUNK
        READ(2,30)
        READ(2,50) DENSITY_START, DENSITY_STOP, DENSITY_DELTA, T_DENSITY   
        READ(2,30)
        READ(2,60) CROWN_HGHT_START, CROWN_HGHT_STOP, CROWN_HGHT_DELTA,
     &             T_CROWN_HGHT
        READ(2,30)
        READ(2,70) TRUNK_DIAM_START, TRUNK_DIAM_STOP,TRUNK_DIAM_DELTA,
     &             T_TRUNK_DIAM
        READ(2,30)
        READ(2,80) TRUNK_HGHT_START, TRUNK_HGHT_STOP,TRUNK_HGHT_DELTA,
     &             T_TRUNK_HGHT
        READ(2,30)
        READ(2,70) TRUNK_Dsig_START, TRUNK_Dsig_STOP,TRUNK_Dsig_DELTA,
     &             T_TRUNK_Dsig
        READ(2,30)
        READ(2,80) TRUNK_Hsig_START, TRUNK_Hsig_STOP,TRUNK_Hsig_DELTA,
     &             T_TRUNK_Hsig
        READ(2,90)
        READ(2,100) I_PDF_TRUNK

        READ(2,30)
        READ(2,100) I_PDF_TRUNK_SIZE
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_MG_TRUNK,LOG_MG_TRUNK)
        CALL SET_FLAG(T_RHO_TRUNK,LOG_RHO_TRUNK)
        CALL SET_FLAG(T_DENSITY,LOG_DENSITY)
        CALL SET_FLAG(T_CROWN_HGHT,LOG_CROWN_HGHT)
        CALL SET_FLAG(T_TRUNK_DIAM,LOG_TRUNK_DIAM)
        CALL SET_FLAG(T_TRUNK_HGHT,LOG_TRUNK_HGHT)
        CALL SET_FLAG(T_TRUNK_Dsig,LOG_TRUNK_Dsig)
        CALL SET_FLAG(T_TRUNK_Hsig,LOG_TRUNK_Hsig)
C
        IF( LOG_MG_TRUNK  )THEN
            CALL READ_TABLE(MG_TRUNK_VECT,MG_TRUNK_NUM,11)
C            print *,'11: MG_TRUNK'
            MG_TRUNK = MG_TRUNK_VECT(1)
            MG_TRUNK_START = MG_TRUNK_VECT(1)
        ELSE
            MG_TRUNK_NUM = LOOP(MG_TRUNK_START, MG_TRUNK_STOP,
     &                          MG_TRUNK_DELTA)
            MG_TRUNK = MG_TRUNK_START
        ENDIF
C
        IF( LOG_RHO_TRUNK )THEN
            CALL READ_TABLE(RHO_TRUNK_VECT,RHO_TRUNK_NUM,12)
C            print *,'12: RHO_TRUNK'
            RHO_TRUNK = RHO_TRUNK_VECT(1)
            RHO_TRUNK_START = RHO_TRUNK_VECT(1)
        ELSE
            RHO_TRUNK_NUM = LOOP(RHO_TRUNK_START, RHO_TRUNK_STOP,
     &                           RHO_TRUNK_DELTA)
            RHO_TRUNK = RHO_TRUNK_START
        ENDIF
C
        IF( LOG_DENSITY   )THEN
            CALL READ_TABLE(DENSITY_VECT,DENSITY_NUM,13) 
C            print *,'13: DENSITY'  
            DENSITY = DENSITY_VECT(1)
            DENSITY_START = DENSITY_VECT(1)
        ELSE
            DENSITY_NUM = LOOP(DENSITY_START, DENSITY_STOP,
     &                         DENSITY_DELTA)
            DENSITY = DENSITY_START
        ENDIF
C
        IF( LOG_CROWN_HGHT)THEN
            CALL READ_TABLE(CROWN_HGHT_VECT,CROWN_HGHT_NUM,14)
C            print *,'14: CROWN_HGHT'
            CROWN_HGHT = CROWN_HGHT_VECT(1)
            CROWN_HGHT_START = CROWN_HGHT_VECT(1)
        ELSE
            CROWN_HGHT_NUM = LOOP(CROWN_HGHT_START, CROWN_HGHT_STOP,
     &                            CROWN_HGHT_DELTA)
            CROWN_HGHT = CROWN_HGHT_START
        ENDIF
C
        IF( LOG_TRUNK_DIAM)THEN
            CALL READ_TABLE(TRUNK_DIAM_VECT,TRUNK_DIAM_NUM,15)
C            print *,'15: TRUNK_DIAM'
            TRUNK_DIAM = TRUNK_DIAM_VECT(1)
            TRUNK_DIAM_START = TRUNK_DIAM_VECT(1)
        ELSE
            TRUNK_DIAM_NUM = LOOP(TRUNK_DIAM_START, TRUNK_DIAM_STOP,
     &                            TRUNK_DIAM_DELTA)
            TRUNK_DIAM = TRUNK_DIAM_START
        ENDIF
C
        IF( LOG_TRUNK_HGHT)THEN
            CALL READ_TABLE(TRUNK_HGHT_VECT,TRUNK_HGHT_NUM,16)
C            print *,'16: TRUNK_HGHT'
            TRUNK_HGHT = TRUNK_HGHT_VECT(1)
            TRUNK_HGHT_START = TRUNK_HGHT_VECT(1)
        ELSE
            TRUNK_HGHT_NUM = LOOP(TRUNK_HGHT_START, TRUNK_HGHT_STOP,
     &                            TRUNK_HGHT_DELTA)
            TRUNK_HGHT = TRUNK_HGHT_START
        ENDIF
C
        IF( LOG_TRUNK_Dsig)THEN
            CALL READ_TABLE(TRUNK_Dsig_VECT,TRUNK_Dsig_NUM,17)
C            print *,'17: TRUNK_Dsig'
            TRUNK_Dsig = TRUNK_Dsig_VECT(1)
            TRUNK_Dsig_START = TRUNK_Dsig_VECT(1)
        ELSE
            TRUNK_Dsig_NUM = LOOP(TRUNK_Dsig_START, TRUNK_Dsig_STOP,
     &                            TRUNK_Dsig_DELTA)
            TRUNK_Dsig = TRUNK_Dsig_START
        ENDIF
C
        IF( LOG_TRUNK_Hsig)THEN
            CALL READ_TABLE(TRUNK_Hsig_VECT,TRUNK_Hsig_NUM,18)
C            print *,'18: TRUNK_Hsig'
            TRUNK_Hsig = TRUNK_Hsig_VECT(1)
            TRUNK_Hsig_START = TRUNK_Hsig_VECT(1)
        ELSE
            TRUNK_Hsig_NUM = LOOP(TRUNK_Hsig_START, TRUNK_Hsig_STOP,
     &                            TRUNK_Hsig_DELTA)
            TRUNK_Hsig = TRUNK_Hsig_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(13(/))
20      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
30      FORMAT(5(/))
40      FORMAT(16X,F6.4,19X,F6.4,19X,F6.4,8X,I1)
50      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
60      FORMAT(17X,F6.3,19X,F6.3,19X,F6.3,7X,I1)
70      FORMAT(17X,F6.2,19X,F6.2,19X,F6.2,7X,I1)
80      FORMAT(17X,F6.2,19X,F6.2,19X,F6.2,7X,I1)
90      FORMAT(10(/))
100     FORMAT(17X,I2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read LEAF.INPUT         *****************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_LEAF
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   MG_LEAF_START =   starting value of leaf gravimetric moisture
C   MG_LEAF_STOP  =   stopping value of leaf gravimetric moisture
C   MG_LEAF_DELTA =   increment for stepping leaf grav. moisture
C
C   RHO_LEAF_START =   starting value of leaf material dry density
C   RHO_LEAF_STOP  =   stopping value of leaf material dry density
C   RHO_LEAF_DELTA =   increment for stepping material dry density
C
C   LEAF_DENS_START = starting value of leaf number density(leaves/m**3)
C   LEAF_DENS_STOP  = stopping value of leaf number density(leaves/m**3)
C   LEAF_DENS_DELTA = increment for stepping leaf number density
C
C   LEAF_DIAM_START =   starting value of leaf diameter (centimeters)
C   LEAF_DIAM_STOP  =   stopping value of leaf diameter (centimeters)
C   LEAF_DIAM_DELTA =   increment for stepping diameter (centimeters)
C
C   LEAF_TAU_START =   starting value of leaf thickness (centimeters)
C   LEAF_TAU_STOP  =   stopping value of leaf thickness (centimeters)
C   LEAF_TAU_DELTA =   increment for stepping thickness (centimeters)
C
C   MG_LEAF_NUM    = number of leaf moisture steps
C   RHO_LEAF_NUM   = number of leaf dry density steps
C   LEAF_DENS_NUM  = number of leaf number density steps
C   LEAF_DIAM_NUM  = number of leaf diameter steps
C   LEAF_TAU_NUM   = number of leaf thickness steps
C
C   LOG_MG_LEAF    = logical variable indicating look-up table data
C   LOG_RHO_LEAF   = logical variable indicating look-up table data
C   LOG_DENS_LEAF  = logical variable indicating look-up table data
C   LOG_LEAF_DIAM  = logical variable indicating look-up table data
C   LOG_LEAF_TAU   = logical variable indicating look-up table data
C
C   MG_LEAF_VECT   = vector containing look-up table values
C   RHO_LEAF_VECT  = vector containing look-up table values
C   LEAF_DENS_VECT  = vector containing look-up table values
C   LEAF_DIAM_VECT = vector containing look-up table values
C   LEAF_TAU_VECT  = vector containing look-up table values
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_MG_LEAF, T_RHO_LEAF, T_LEAF_DENS   
        INTEGER T_LEAF_DIAM, T_LEAF_TAU
        INTEGER MG_LEAF_NUM, RHO_LEAF_NUM, LEAF_DENS_NUM   
        INTEGER LEAF_DIAM_NUM, LEAF_TAU_NUM
        INTEGER LOOP
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,I_PDF_NDL,I_PDF_TRUNK
C
        REAL MG_LEAF_START, MG_LEAF_STOP, MG_LEAF_DELTA
        REAL RHO_LEAF_START, RHO_LEAF_STOP, RHO_LEAF_DELTA
        REAL LEAF_DENS_START, LEAF_DENS_STOP, LEAF_DENS_DELTA
        REAL LEAF_DIAM_START, LEAF_DIAM_STOP, LEAF_DIAM_DELTA
        REAL LEAF_TAU_START, LEAF_TAU_STOP, LEAF_TAU_DELTA
        REAL MG_LEAF_VECT(I_VECT_LEN)
        REAL RHO_LEAF_VECT(I_VECT_LEN)
        REAL LEAF_DENS_VECT(I_VECT_LEN)
        REAL LEAF_DIAM_VECT(I_VECT_LEN)
        REAL LEAF_TAU_VECT(I_VECT_LEN)
        REAL MG_LEAF, RHO_LEAF, LEAF_DENS
        REAL LEAF_DIAM,LEAF_TAU
C
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
        COMMON /LEAF_START/ MG_LEAF_START, RHO_LEAF_START,
     &             LEAF_DENS_START,LEAF_DIAM_START,LEAF_TAU_START
        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='leaf.input',STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) MG_LEAF_START, MG_LEAF_STOP, MG_LEAF_DELTA,
     &             T_MG_LEAF
        READ(2,30)
        READ(2,40) RHO_LEAF_START, RHO_LEAF_STOP, RHO_LEAF_DELTA,
     &             T_RHO_LEAF
        READ(2,30)
        READ(2,60) LEAF_TAU_START, LEAF_TAU_STOP,LEAF_TAU_DELTA,
     &             T_LEAF_TAU
        READ(2,30)
        READ(2,50) LEAF_DIAM_START, LEAF_DIAM_STOP,LEAF_DIAM_DELTA,
     &             T_LEAF_DIAM

        READ(2,30)
        READ(2,70) LEAF_DENS_START, LEAF_DENS_STOP, LEAF_DENS_DELTA,
     &             T_LEAF_DENS   

        READ(2,80)
        READ(2,90) I_PDF_LEAF
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_MG_LEAF,LOG_MG_LEAF)
        CALL SET_FLAG(T_RHO_LEAF,LOG_RHO_LEAF)
        CALL SET_FLAG(T_LEAF_DENS,LOG_LEAF_DENS)
        CALL SET_FLAG(T_LEAF_DIAM,LOG_LEAF_DIAM)
        CALL SET_FLAG(T_LEAF_TAU,LOG_LEAF_TAU)
C
        IF( LOG_MG_LEAF  )THEN
            CALL READ_TABLE(MG_LEAF_VECT,MG_LEAF_NUM,19)
C            print *,'19: MG_LEAF'
            MG_LEAF = MG_LEAF_VECT(1)
            MG_LEAF_START = MG_LEAF_VECT(1)
        ELSE
            MG_LEAF_NUM = LOOP(MG_LEAF_START, MG_LEAF_STOP,
     &                          MG_LEAF_DELTA)
            MG_LEAF = MG_LEAF_START
        ENDIF
C
        IF( LOG_RHO_LEAF)THEN
            CALL READ_TABLE(RHO_LEAF_VECT,RHO_LEAF_NUM,20)
C            print *,'20: RHO_LEAF'
            RHO_LEAF = RHO_LEAF_VECT(1)
            RHO_LEAF_START = RHO_LEAF_VECT(1)
        ELSE
            RHO_LEAF_NUM = LOOP(RHO_LEAF_START, RHO_LEAF_STOP,
     &                           RHO_LEAF_DELTA)
            RHO_LEAF = RHO_LEAF_START
        ENDIF
C
        IF( LOG_LEAF_DENS)THEN
            CALL READ_TABLE(LEAF_DENS_VECT,LEAF_DENS_NUM,21)  
C            print *,'21: LEAF_DENS'  
            LEAF_DENS = LEAF_DENS_VECT(1)
            LEAF_DENS_START = LEAF_DENS_VECT(1)
        ELSE
            LEAF_DENS_NUM = LOOP(LEAF_DENS_START, LEAF_DENS_STOP,
     &                         LEAF_DENS_DELTA)
            LEAF_DENS = LEAF_DENS_START
        ENDIF
C
        IF( LOG_LEAF_DIAM)THEN
            CALL READ_TABLE(LEAF_DIAM_VECT,LEAF_DIAM_NUM,22)
C            print *,'22: LEAF_DIAM'
            LEAF_DIAM = LEAF_DIAM_VECT(1)
            LEAF_DIAM_START = LEAF_DIAM_VECT(1)
        ELSE
            LEAF_DIAM_NUM = LOOP(LEAF_DIAM_START, LEAF_DIAM_STOP,
     &                            LEAF_DIAM_DELTA)
            LEAF_DIAM = LEAF_DIAM_START
        ENDIF
C
        IF( LOG_LEAF_TAU)THEN
            CALL READ_TABLE(LEAF_TAU_VECT,LEAF_TAU_NUM,23)
C            print *,'23: LEAF_TAU'
            LEAF_TAU = LEAF_TAU_VECT(1)
            LEAF_TAU_START = LEAF_TAU_VECT(1)
        ELSE
            LEAF_TAU_NUM = LOOP(LEAF_TAU_START, LEAF_TAU_STOP,
     &                            LEAF_TAU_DELTA)
            LEAF_TAU = LEAF_TAU_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
30      FORMAT(5(/))
40      FORMAT(16X,F6.4,19X,F6.4,19X,F6.4,8X,I1)
50      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
60      FORMAT(17X,F5.2,20X,F5.2,20X,F4.2,9X,I1)
70      FORMAT(17X,F7.2,18X,F7.2,18X,F7.2,5X,I1)
80      FORMAT(7(/))
90      FORMAT(17X,I2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read NEEDLE.INPUT       *****************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_NEEDLE
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   MG_NDL_START =   starting value of needle gravimetric moisture
C   MG_NDL_STOP  =   stopping value of needle gravimetric moisture
C   MG_NDL_DELTA =   increment for stepping needle grav. moisture
C
C   RHO_NDL_START =   starting value of needle material dry density
C   RHO_NDL_STOP  =   stopping value of needle material dry density
C   RHO_NDL_DELTA =   increment for stepping material dry density
C
C   NDL_DENS_START = start value of needle number density (leaves/m**3)
C   NDL_DENS_STOP  = stop value of needle number density (leaves/m**3)
C   NDL_DENS_DELTA = increment for stepping needle number density
C
C   NDL_DIAM_START =   starting value of needle diameter (centimeters)
C   NDL_DIAM_STOP  =   stopping value of needle diameter (centimeters)
C   NDL_DIAM_DELTA =   increment for stepping diameter (centimeters)
C
C   NDL_LNG_START =   starting value of needle length (centimeters)
C   NDL_LNG_STOP  =   stopping value of needle length (centimeters)
C   NDL_LNG_DELTA =   increment for stepping length (centimeters)
C
C   MG_NDL_NUM    = number of needle moisture steps
C   RHO_NDL_NUM   = number of needle dry density steps
C   NDL_DENS_NUM  = number of needle number density steps
C   NDL_DIAM_NUM  = number of needle diameter steps
C   NDL_TAU_NUM   = number of needle length steps
C
C   LOG_MG_NDL    = logical variable indicating look-up table data
C   LOG_RHO_NDL   = logical variable indicating look-up table data
C   LOG_DENS_NDL  = logical variable indicating look-up table data
C   LOG_NDL_DIAM  = logical variable indicating look-up table data
C   LOG_NDL_TAU   = logical variable indicating look-up table data
C
C   MG_NDL_VECT   = vector containing look-up table values
C   RHO_NDL_VECT  = vector containing look-up table values
C   NDL_DENS_VECT  = vector containing look-up table values
C   NDL_DIAM_VECT = vector containing look-up table values
C   NDL_TAU_VECT  = vector containing look-up table values
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_MG_NDL, T_RHO_NDL, T_NDL_DENS   
        INTEGER T_NDL_DIAM, T_NDL_LNG
        INTEGER MG_NDL_NUM, RHO_NDL_NUM, NDL_DENS_NUM   
        INTEGER NDL_DIAM_NUM, NDL_LNG_NUM
        INTEGER LOOP
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,I_PDF_NDL,I_PDF_TRUNK
C
        REAL MG_NDL_START, MG_NDL_STOP, MG_NDL_DELTA
        REAL RHO_NDL_START, RHO_NDL_STOP, RHO_NDL_DELTA
        REAL NDL_DENS_START, NDL_DENS_STOP, NDL_DENS_DELTA
        REAL NDL_DIAM_START, NDL_DIAM_STOP, NDL_DIAM_DELTA
        REAL NDL_LNG_START, NDL_LNG_STOP, NDL_LNG_DELTA
        REAL MG_NDL_VECT(I_VECT_LEN)
        REAL RHO_NDL_VECT(I_VECT_LEN)
        REAL NDL_DENS_VECT(I_VECT_LEN)
        REAL NDL_DIAM_VECT(I_VECT_LEN)
        REAL NDL_LNG_VECT(I_VECT_LEN)
        REAL MG_NDL, RHO_NDL, NDL_DENS
        REAL NDL_DIAM,NDL_LNG
C
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
        COMMON /NDL_START/ MG_NDL_START, RHO_NDL_START,
     &             NDL_DENS_START, NDL_DIAM_START, NDL_LNG_START
        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='needle.input',STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) MG_NDL_START, MG_NDL_STOP, MG_NDL_DELTA,
     &             T_MG_NDL
        READ(2,30)
        READ(2,40) RHO_NDL_START, RHO_NDL_STOP, RHO_NDL_DELTA,
     &             T_RHO_NDL
        READ(2,30)
        READ(2,50) NDL_LNG_START, NDL_LNG_STOP, NDL_LNG_DELTA,
     &             T_NDL_LNG
        READ(2,30)
        READ(2,60) NDL_DIAM_START, NDL_DIAM_STOP, NDL_DIAM_DELTA,
     &             T_NDL_DIAM
        READ(2,30)
        READ(2,70) NDL_DENS_START, NDL_DENS_STOP, NDL_DENS_DELTA,
     &             T_NDL_DENS   
        READ(2,80)
        READ(2,90) I_PDF_NDL
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_MG_NDL,LOG_MG_NDL)
        CALL SET_FLAG(T_RHO_NDL,LOG_RHO_NDL)
        CALL SET_FLAG(T_NDL_DENS,LOG_NDL_DENS)
        CALL SET_FLAG(T_NDL_DIAM,LOG_NDL_DIAM)
        CALL SET_FLAG(T_NDL_LNG,LOG_NDL_LNG)
C
        IF( LOG_MG_NDL)THEN
            CALL READ_TABLE(MG_NDL_VECT,MG_NDL_NUM,24)
C            print *,'24: MG_NDL'
            MG_NDL = MG_NDL_VECT(1)
            MG_NDL_START = MG_NDL_VECT(1)
        ELSE
            MG_NDL_NUM = LOOP(MG_NDL_START, MG_NDL_STOP,
     &                          MG_NDL_DELTA)
            MG_NDL = MG_NDL_START
        ENDIF
C
        IF( LOG_RHO_NDL)THEN
            CALL READ_TABLE(RHO_NDL_VECT,RHO_NDL_NUM,25)
C            print *,'25: RHO_NDL'
            RHO_NDL = RHO_NDL_VECT(1)
            RHO_NDL_START = RHO_NDL_VECT(1)
        ELSE
            RHO_NDL_NUM = LOOP(RHO_NDL_START, RHO_NDL_STOP,
     &                           RHO_NDL_DELTA)
            RHO_NDL = RHO_NDL_START
        ENDIF
C
        IF( LOG_NDL_DENS)THEN
            CALL READ_TABLE(NDL_DENS_VECT,NDL_DENS_NUM,26) 
C            print *,'26: NDL_DENS'  
            NDL_DENS = NDL_DENS_VECT(1)
            NDL_DENS_START = NDL_DENS_VECT(1)
        ELSE
            NDL_DENS_NUM = LOOP(NDL_DENS_START, NDL_DENS_STOP,
     &                         NDL_DENS_DELTA)
            NDL_DENS = NDL_DENS_START
        ENDIF
C
        IF( LOG_NDL_DIAM)THEN
            CALL READ_TABLE(NDL_DIAM_VECT,NDL_DIAM_NUM,27)
C            print *,'27: NDL_DIAM'
            NDL_DIAM = NDL_DIAM_VECT(1)
            NDL_DIAM_START = NDL_DIAM_VECT(1)
        ELSE
            NDL_DIAM_NUM = LOOP(NDL_DIAM_START, NDL_DIAM_STOP,
     &                            NDL_DIAM_DELTA)
            NDL_DIAM = NDL_DIAM_START
        ENDIF
C
        IF( LOG_NDL_LNG)THEN
            CALL READ_TABLE(NDL_LNG_VECT,NDL_LNG_NUM,28)
C            print *,'28: NDL_LNG'
            NDL_LNG = NDL_LNG_VECT(1)
            NDL_LNG_START = NDL_LNG_VECT(1)
        ELSE
            NDL_LNG_NUM = LOOP(NDL_LNG_START, NDL_LNG_STOP,
     &                            NDL_LNG_DELTA)
            NDL_LNG = NDL_LNG_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
30      FORMAT(5(/))
40      FORMAT(16X,F6.4,19X,F6.4,19X,F6.4,8X,I1)
50      FORMAT(16X,F7.4,18X,F7.4,18X,F7.4,7X,I1)
60      FORMAT(17X,F5.2,20X,F5.2,20X,F4.2,9X,I1)
70      FORMAT(17X,F8.1,17X,F8.1,17X,F8.1,5X,I1)
80      FORMAT(9(/))
90      FORMAT(17X,I2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read BRANCH_PRIMARY.INPUT ***************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_PRIMARY_BRANCH
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   MG_BR1_START =   starting value of branch gravimetric moisture
C   MG_BR1_STOP  =   stopping value of branch gravimetric moisture
C   MG_BR1_DELTA =   increment for stepping branch grav. moisture
C
C   RHO_BR1_START =   starting value of branch material dry density
C   RHO_BR1_STOP  =   stopping value of branch material dry density
C   RHO_BR1_DELTA =   increment for stepping material dry density
C
C   BR1_DENS_START = start value of branch number density(branches/m**3)
C   BR1_DENS_STOP  = stop value of branch number density (branches/m**3)
C   BR1_DENS_DELTA = increment for stepping branch number density
C
C   BR1_DIAM_START =   starting value of branch diameter (centimeters)
C   BR1_DIAM_STOP  =   stopping value of branch diameter (centimeters)
C   BR1_DIAM_DELTA =   increment for stepping diameter (centimeters)
C
C   BR1_LNG_START =   starting value of branch length (meters)
C   BR1_LNG_STOP  =   stopping value of branch length (meters)
C   BR1_LNG_DELTA =   increment for stepping length (meters)
C
C   MG_BR1_NUM    = number of branch moisture steps
C   RHO_BR1_NUM   = number of branch dry density steps
C   BR1_DENS_NUM  = number of branch number density steps
C   BR1_DIAM_NUM  = number of branch diameter steps
C   BR1_TAU_NUM   = number of branch length steps
C
C   LOG_MG_BR1    = logical variable indicating look-up table data
C   LOG_RHO_BR1   = logical variable indicating look-up table data
C   LOG_DENS_BR1  = logical variable indicating look-up table data
C   LOG_BR1_DIAM  = logical variable indicating look-up table data
C   LOG_BR1_TAU   = logical variable indicating look-up table data
C
C   MG_BR1_VECT   = vector containing look-up table values
C   RHO_BR1_VECT  = vector containing look-up table values
C   BR1_DENS_VECT = vector containing look-up table values
C   BR1_DIAM_VECT = vector containing look-up table values
C   BR1_LNG_VECT  = vector containing look-up table values
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_MG_BR1, T_RHO_BR1, T_BR1_DENS   
        INTEGER T_BR1_DIAM, T_BR1_LNG
        INTEGER MG_BR1_NUM, RHO_BR1_NUM, BR1_DENS_NUM   
        INTEGER BR1_DIAM_NUM, BR1_LNG_NUM
        INTEGER LOOP
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_NDL,I_PDF_TRUNK
C
        REAL MG_BR1_START, MG_BR1_STOP, MG_BR1_DELTA
        REAL RHO_BR1_START, RHO_BR1_STOP, RHO_BR1_DELTA
        REAL BR1_DENS_START, BR1_DENS_STOP, BR1_DENS_DELTA
        REAL BR1_DIAM_START, BR1_DIAM_STOP, BR1_DIAM_DELTA
        REAL BR1_LNG_START, BR1_LNG_STOP, BR1_LNG_DELTA
        REAL MG_BR1_VECT(I_VECT_LEN)
        REAL RHO_BR1_VECT(I_VECT_LEN)
        REAL BR1_DENS_VECT(I_VECT_LEN)
        REAL BR1_DIAM_VECT(I_VECT_LEN)
        REAL BR1_LNG_VECT(I_VECT_LEN)
        REAL MG_BR1, RHO_BR1, BR1_DENS
        REAL BR1_DIAM, BR1_LNG
C
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
        COMMON /BR1_START/ MG_BR1_START, RHO_BR1_START,
     &            BR1_DENS_START, BR1_DIAM_START, BR1_LNG_START
        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='branch_primary.input',STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) MG_BR1_START, MG_BR1_STOP, MG_BR1_DELTA,
     &             T_MG_BR1
        READ(2,30)
        READ(2,40) RHO_BR1_START, RHO_BR1_STOP, RHO_BR1_DELTA,
     &             T_RHO_BR1
        READ(2,30)
        READ(2,50) BR1_LNG_START, BR1_LNG_STOP, BR1_LNG_DELTA,
     &             T_BR1_LNG
        READ(2,30)
        READ(2,60) BR1_DIAM_START, BR1_DIAM_STOP, BR1_DIAM_DELTA,
     &             T_BR1_DIAM
        READ(2,30)
        READ(2,70) BR1_DENS_START, BR1_DENS_STOP, BR1_DENS_DELTA,
     &             T_BR1_DENS   
        READ(2,80)
        READ(2,90) I_PDF_BR_1
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_MG_BR1,LOG_MG_BR1)
        CALL SET_FLAG(T_RHO_BR1,LOG_RHO_BR1)
        CALL SET_FLAG(T_BR1_DENS,LOG_BR1_DENS)
        CALL SET_FLAG(T_BR1_DIAM,LOG_BR1_DIAM)
        CALL SET_FLAG(T_BR1_LNG,LOG_BR1_LNG)
C
        IF( LOG_MG_BR1)THEN
            CALL READ_TABLE(MG_BR1_VECT,MG_BR1_NUM,29)
C            print *,'29: MG_BR1'
            MG_BR1 = MG_BR1_VECT(1)
            MG_BR1_START = MG_BR1_VECT(1)
        ELSE
            MG_BR1_NUM = LOOP(MG_BR1_START, MG_BR1_STOP,
     &                          MG_BR1_DELTA)
            MG_BR1 = MG_BR1_START
        ENDIF
C
        IF( LOG_RHO_BR1)THEN
            CALL READ_TABLE(RHO_BR1_VECT,RHO_BR1_NUM,30)
C            print *,'30: RHO_BR1'
            RHO_BR1 = RHO_BR1_VECT(1)
            RHO_BR1_START = RHO_BR1_VECT(1)
        ELSE
            RHO_BR1_NUM = LOOP(RHO_BR1_START, RHO_BR1_STOP,
     &                           RHO_BR1_DELTA)
            RHO_BR1 = RHO_BR1_START
        ENDIF
C
        IF( LOG_BR1_DENS)THEN
            CALL READ_TABLE(BR1_DENS_VECT,BR1_DENS_NUM,31) 
C            print *,'31: BR1_DENS'   
            BR1_DENS = BR1_DENS_VECT(1)
            BR1_DENS_START = BR1_DENS_VECT(1)
        ELSE
            BR1_DENS_NUM = LOOP(BR1_DENS_START, BR1_DENS_STOP,
     &                         BR1_DENS_DELTA)
            BR1_DENS = BR1_DENS_START
        ENDIF
C
        IF( LOG_BR1_DIAM)THEN
            CALL READ_TABLE(BR1_DIAM_VECT,BR1_DIAM_NUM,32)
C            print *,'32: BR1_DIAM'
            BR1_DIAM = BR1_DIAM_VECT(1)
            BR1_DIAM_START = BR1_DIAM_VECT(1)
        ELSE
            BR1_DIAM_NUM = LOOP(BR1_DIAM_START, BR1_DIAM_STOP,
     &                            BR1_DIAM_DELTA)
            BR1_DIAM = BR1_DIAM_START
        ENDIF
C
        IF( LOG_BR1_LNG)THEN
            CALL READ_TABLE(BR1_LNG_VECT,BR1_LNG_NUM,33)
C            print *,'33: BR1_LNG'
            BR1_LNG = BR1_LNG_VECT(1)
            BR1_LNG_START = BR1_LNG_VECT(1)
        ELSE
            BR1_LNG_NUM = LOOP(BR1_LNG_START, BR1_LNG_STOP,
     &                         BR1_LNG_DELTA)
            BR1_LNG = BR1_LNG_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(13(/))
20      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
30      FORMAT(5(/))
40      FORMAT(16X,F6.4,19X,F6.4,19X,F6.4,8X,I1)
50      FORMAT(17X,F5.2,20X,F5.2,20X,F5.2,8X,I1)
60      FORMAT(17X,F6.3,19X,F6.3,19X,F6.3,7X,I1)
70      FORMAT(17X,F6.2,19X,F6.2,19X,F6.2,7X,I1)
80      FORMAT(9(/))
90      FORMAT(17X,I2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read BRANCH_SECONDARY.INPUT *************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   SET_FLAG                               ***
C***    Called functions:     LOOP                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_SECONDARY_BRANCH
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   MG_BR2_START =   starting value of branch gravimetric moisture
C   MG_BR2_STOP  =   stopping value of branch gravimetric moisture
C   MG_BR2_DELTA =   increment for stepping branch grav. moisture
C
C   RHO_BR2_START =   starting value of branch material dry density
C   RHO_BR2_STOP  =   stopping value of branch material dry density
C   RHO_BR2_DELTA =   increment for stepping material dry density
C
C   BR2_DENS_START = start value of branch number density(branches/m**3)
C   BR2_DENS_STOP  = stop value of branch number density (branches/m**3)
C   BR2_DENS_DELTA = increment for stepping branch number density
C
C   BR2_DIAM_START =   starting value of branch diameter (centimeters)
C   BR2_DIAM_STOP  =   stopping value of branch diameter (centimeters)
C   BR2_DIAM_DELTA =   increment for stepping diameter (centimeters)
C
C   BR2_LNG_START =   starting value of branch length (meters)
C   BR2_LNG_STOP  =   stopping value of branch length (meters)
C   BR2_LNG_DELTA =   increment for stepping length (meters)
C
C   MG_BR2_NUM    = number of branch moisture steps
C   RHO_BR2_NUM   = number of branch dry density steps
C   BR2_DENS_NUM  = number of branch number density steps
C   BR2_DIAM_NUM  = number of branch diameter steps
C   BR2_TAU_NUM   = number of branch length steps
C
C   LOG_MG_BR2    = logical variable indicating look-up table data
C   LOG_RHO_BR2   = logical variable indicating look-up table data
C   LOG_DENS_BR2  = logical variable indicating look-up table data
C   LOG_BR2_DIAM  = logical variable indicating look-up table data
C   LOG_BR2_TAU   = logical variable indicating look-up table data
C
C   MG_BR2_VECT   = vector containing look-up table values
C   RHO_BR2_VECT  = vector containing look-up table values
C   BR2_DENS_VECT = vector containing look-up table values
C   BR2_DIAM_VECT = vector containing look-up table values
C   BR2_LNG_VECT  = vector containing look-up table values
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER T_MG_BR2, T_RHO_BR2, T_BR2_DENS   
        INTEGER T_BR2_DIAM, T_BR2_LNG
        INTEGER MG_BR2_NUM, RHO_BR2_NUM, BR2_DENS_NUM   
        INTEGER BR2_DIAM_NUM, BR2_LNG_NUM
        INTEGER LOOP
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,I_PDF_NDL,I_PDF_TRUNK
C
        REAL MG_BR2_START, MG_BR2_STOP, MG_BR2_DELTA
        REAL RHO_BR2_START, RHO_BR2_STOP, RHO_BR2_DELTA
        REAL BR2_DENS_START, BR2_DENS_STOP, BR2_DENS_DELTA
        REAL BR2_DIAM_START, BR2_DIAM_STOP, BR2_DIAM_DELTA
        REAL BR2_LNG_START, BR2_LNG_STOP, BR2_LNG_DELTA
        REAL MG_BR2_VECT(I_VECT_LEN)
        REAL RHO_BR2_VECT(I_VECT_LEN)
        REAL BR2_DENS_VECT(I_VECT_LEN)
        REAL BR2_DIAM_VECT(I_VECT_LEN)
        REAL BR2_LNG_VECT(I_VECT_LEN)
        REAL MG_BR2, RHO_BR2, BR2_DENS
        REAL BR2_DIAM, BR2_LNG
C
        LOGICAL LOG_MG_BR2, LOG_RHO_BR2, LOG_BR2_DENS
        LOGICAL LOG_BR2_DIAM, LOG_BR2_LNG
C
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        COMMON /R_BR2_DELTA/ MG_BR2_DELTA,RHO_BR2_DELTA,BR2_DENS_DELTA,
     &                        BR2_DIAM_DELTA, BR2_LNG_DELTA
        COMMON /I_BR2/ MG_BR2_NUM, RHO_BR2_NUM, BR2_DENS_NUM,
     &                  BR2_DIAM_NUM, BR2_LNG_NUM
        COMMON /R_BR2_VECT/MG_BR2_VECT, RHO_BR2_VECT, BR2_DENS_VECT,
     &                      BR2_DIAM_VECT, BR2_LNG_VECT
        COMMON /L_BR2/ LOG_MG_BR2, LOG_RHO_BR2, LOG_BR2_DENS,
     &                      LOG_BR2_DIAM, LOG_BR2_LNG
        COMMON /BR2_START/ MG_BR2_START, RHO_BR2_START,
     &            BR2_DENS_START, BR2_DIAM_START, BR2_LNG_START
        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
C***********************************************************************
C   READ DATA
C***********************************************************************
C
        OPEN(UNIT=2,FILE='branch_secondary.input',STATUS='UNKNOWN')
C
        READ(2,10)
        READ(2,20) MG_BR2_START, MG_BR2_STOP, MG_BR2_DELTA,
     &             T_MG_BR2
        READ(2,30)
        READ(2,40) RHO_BR2_START, RHO_BR2_STOP, RHO_BR2_DELTA,
     &             T_RHO_BR2
        READ(2,30)
        READ(2,50) BR2_LNG_START, BR2_LNG_STOP, BR2_LNG_DELTA,
     &             T_BR2_LNG
        READ(2,30)
        READ(2,60) BR2_DIAM_START, BR2_DIAM_STOP, BR2_DIAM_DELTA,
     &             T_BR2_DIAM
        READ(2,30)
        READ(2,70) BR2_DENS_START, BR2_DENS_STOP, BR2_DENS_DELTA,
     &             T_BR2_DENS   
        READ(2,80)
        READ(2,90) I_PDF_BR_2
C
        CLOSE(2)
C
C***********************************************************************
C   SET FLAGS, COUNTERS, AND INITIALIZE VARIABLES
C***********************************************************************
C
        CALL SET_FLAG(T_MG_BR2,LOG_MG_BR2)
        CALL SET_FLAG(T_RHO_BR2,LOG_RHO_BR2)
        CALL SET_FLAG(T_BR2_DENS,LOG_BR2_DENS)
        CALL SET_FLAG(T_BR2_DIAM,LOG_BR2_DIAM)
        CALL SET_FLAG(T_BR2_LNG,LOG_BR2_LNG)
C
        IF( LOG_MG_BR2)THEN
            CALL READ_TABLE(MG_BR2_VECT,MG_BR2_NUM,34)
C            print *,'34: MG_BR2'
            MG_BR2 = MG_BR2_VECT(1)
            MG_BR2_START = MG_BR2_VECT(1)
        ELSE
            MG_BR2_NUM = LOOP(MG_BR2_START, MG_BR2_STOP,
     &                          MG_BR2_DELTA)
            MG_BR2 = MG_BR2_START
        ENDIF
C
        IF( LOG_RHO_BR2)THEN
            CALL READ_TABLE(RHO_BR2_VECT,RHO_BR2_NUM,35)
C            print *,'35: RHO_BR2'
            RHO_BR2 = RHO_BR2_VECT(1)
            RHO_BR2_START = RHO_BR2_VECT(1)
        ELSE
            RHO_BR2_NUM = LOOP(RHO_BR2_START, RHO_BR2_STOP,
     &                           RHO_BR2_DELTA)
            RHO_BR2 = RHO_BR2_START
        ENDIF
C
        IF( LOG_BR2_DENS)THEN
            CALL READ_TABLE(BR2_DENS_VECT,BR2_DENS_NUM,36) 
C            print *,'36: BR2_DENS'   
            BR2_DENS = BR2_DENS_VECT(1)
            BR2_DENS_START = BR2_DENS_VECT(1)
        ELSE
            BR2_DENS_NUM = LOOP(BR2_DENS_START, BR2_DENS_STOP,
     &                         BR2_DENS_DELTA)
            BR2_DENS = BR2_DENS_START
        ENDIF
C
        IF( LOG_BR2_DIAM)THEN
            CALL READ_TABLE(BR2_DIAM_VECT,BR2_DIAM_NUM,37)
C            print *,'37: BR2_DIAM'
            BR2_DIAM = BR2_DIAM_VECT(1)
            BR2_DIAM_START = BR2_DIAM_VECT(1)
        ELSE
            BR2_DIAM_NUM = LOOP(BR2_DIAM_START, BR2_DIAM_STOP,
     &                            BR2_DIAM_DELTA)
            BR2_DIAM = BR2_DIAM_START
        ENDIF
C
        IF( LOG_BR2_LNG)THEN
            CALL READ_TABLE(BR2_LNG_VECT,BR2_LNG_NUM,38)
C            print *,'38: BR2_LNG'
            BR2_LNG = BR2_LNG_VECT(1)
            BR2_LNG_START = BR2_LNG_VECT(1)
        ELSE
            BR2_LNG_NUM = LOOP(BR2_LNG_START, BR2_LNG_STOP,
     &                         BR2_LNG_DELTA)
            BR2_LNG = BR2_LNG_START
        ENDIF
C
C***********************************************************************
C   FORMAT STATEMENTS
C***********************************************************************
C
10      FORMAT(13(/))
20      FORMAT(17X,F6.4,19X,F6.4,19X,F6.4,7X,I1)
30      FORMAT(5(/))
40      FORMAT(16X,F6.4,19X,F6.4,19X,F6.4,8X,I1)
50      FORMAT(17X,F5.2,20X,F5.2,20X,F5.2,8X,I1)
60      FORMAT(17X,F6.3,19X,F6.3,19X,F6.3,7X,I1)
70      FORMAT(17X,F6.2,19X,F6.2,19X,F6.2,7X,I1)
80      FORMAT(9(/))
90      FORMAT(17X,I2)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read deiectric tables        ************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   none                                   ***
C***    Called functions:     none                                   ***
C***********************************************************************
C
        SUBROUTINE READ_EPS_TABLE(J)
        save
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER J, I, NUMBER
        REAL A, B
        COMPLEX EPS_TR_VECT(I_VECT_LEN) ,EPS_TR_VECT_C(I_VECT_LEN)
        COMPLEX EPS_LF_VECT(I_VECT_LEN) ,EPS_LF_VECT_C(I_VECT_LEN)
        COMPLEX EPS_BR_1_VECT(I_VECT_LEN) ,EPS_BR_1_VECT_C(I_VECT_LEN)
        COMPLEX EPS_BR_2_VECT(I_VECT_LEN) ,EPS_BR_2_VECT_C(I_VECT_LEN)
        COMPLEX EPS_NDL_VECT(I_VECT_LEN) ,EPS_NDL_VECT_C(I_VECT_LEN)
        COMPLEX EPS_GRND_VECT(I_VECT_LEN) ,EPS_GRND_VECT_C(I_VECT_LEN)
        COMPLEX EPS_SNOW_VECT(I_VECT_LEN) ,EPS_SNOW_VECT_C(I_VECT_LEN)
C
        COMMON /C_EPS_VECT/ EPS_TR_VECT, EPS_LF_VECT, EPS_BR_1_VECT,
     &      EPS_BR_2_VECT, EPS_NDL_VECT, EPS_GRND_VECT,EPS_SNOW_VECT
        COMMON /C_EPS_VECT_C/ EPS_TR_VECT_C, EPS_LF_VECT_C,
     &              EPS_BR_1_VECT_C, EPS_BR_2_VECT_C, EPS_NDL_VECT_C,
     &              EPS_GRND_VECT_C, EPS_SNOW_VECT_C
C
        INTEGER EPS_NUM(N_EPS)
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMMON /I_EPS_NUM/ EPS_NUM
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
C***********************************************************************
C
        IF(J.EQ.1) THEN
           OPEN(UNIT=2,FILE='dielectric_trunk_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 110 I=1,NUMBER
                READ(2,30) A,B
                EPS_TR_VECT(I) = CMPLX(A,-B)
                EPS_TR_VECT_C(I) = CMPLX(A,B)
110        CONTINUE
           EPSILONR(5) = EPS_TR_VECT(1)
           EPSILONRC(5) = EPS_TR_VECT_C(1)
           EPS_NUM(5) = NUMBER
C
        ELSE IF(J.EQ.2) THEN
C
         OPEN(UNIT=2,FILE='dielectric_branch1_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 120 I=1,NUMBER
                READ(2,30) A,B
                EPS_BR_1_VECT(I) = CMPLX(A,-B)
                EPS_BR_1_VECT_C(I) = CMPLX(A,B)
120        CONTINUE
           EPSILONR(6) = EPS_BR_1_VECT(1)
           EPSILONRC(6) = EPS_BR_1_VECT_C(1)
           EPS_NUM(6) = NUMBER
C
        ELSE IF(J.EQ.3) THEN
C
         OPEN(UNIT=2,FILE='dielectric_branch2_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 130 I=1,NUMBER
                READ(2,30) A,B
                EPS_BR_2_VECT(I) = CMPLX(A,-B)
                EPS_BR_2_VECT_C(I) = CMPLX(A,B)
130        CONTINUE
           EPSILONR(7) = EPS_BR_2_VECT(1)
           EPSILONRC(7) = EPS_BR_2_VECT_C(1)
           EPS_NUM(7) = NUMBER
C
        ELSE IF(J.EQ.4) THEN
C
         OPEN(UNIT=2,FILE='dielectric_leaf_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 140 I=1,NUMBER
                READ(2,30) A,B
                EPS_LF_VECT(I) = CMPLX(A,-B)
                EPS_LF_VECT_C(I) = CMPLX(A,B)
140        CONTINUE
           EPSILONR(4) = EPS_LF_VECT(1)
           EPSILONRC(4) = EPS_LF_VECT_C(1)
           EPS_NUM(4) = NUMBER
C
        ELSE IF(J.EQ.5) THEN
C
         OPEN(UNIT=2,FILE='dielectric_needle_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 150 I=1,NUMBER
                READ(2,30) A,B
                EPS_NDL_VECT(I) = CMPLX(A,-B)
                EPS_NDL_VECT_C(I) = CMPLX(A,B)
150        CONTINUE
           EPSILONR(3) = EPS_NDL_VECT(1)
           EPSILONRC(3) = EPS_NDL_VECT_C(1)
           EPS_NUM(3) = NUMBER
C
        ELSE IF(J.EQ.6) THEN
C
         OPEN(UNIT=2,FILE='dielectric_ground_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 160 I=1,NUMBER
                READ(2,30) A,B
                EPS_GRND_VECT(I) = CMPLX(A,-B)
                EPS_GRND_VECT_C(I) = CMPLX(A,B)
160        CONTINUE
           EPSILONR(1) = EPS_GRND_VECT(1)
           EPSILONR(2) = EPS_GRND_VECT(1)
           EPSILONRC(1) = EPS_GRND_VECT_C(1)
           EPSILONRC(2) = EPS_GRND_VECT_C(1)
           EPS_NUM(1) = 0
           EPS_NUM(2) = NUMBER
C
        ELSE IF(J.EQ.(N_CONSTITUENTS+1)) THEN
C
         OPEN(UNIT=2,FILE='dielectric_snow_table.input',
     &          STATUS='UNKNOWN')
           READ(2,5)
           READ(2,20) NUMBER
           READ(2,25)
           DO 170 I=1,NUMBER
                READ(2,30) A,B
                EPS_SNOW_VECT(I) = CMPLX(A,-B)
                EPS_SNOW_VECT_C(I) = CMPLX(A,B)
170        CONTINUE
           EPSILONR(8) = EPS_SNOW_VECT(1)
           EPSILONRC(8) = EPS_SNOW_VECT_C(1)
           EPS_NUM(8) = NUMBER
C
        ELSE 
          CALL WRITE_ERROR(120)
          STOP
        ENDIF
C
        CLOSE(2)
C
C***********************************************************************
C
5       FORMAT(/,/,/,/,/,/,/)
20      FORMAT(33X,I2)
25      FORMAT(1X)
30      FORMAT(19X,F7.3,5X,F7.3)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read PARAMETER_NESTING.INPUT ************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   none                                   ***
C***    Called functions:     none                                   ***
C*********************************************************************** 
C
        SUBROUTINE READ_NESTING
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   PARAM_NUMBER() = VECTOR INDICATING THE NUMBER ASSIGNED TO EACH
C                    INPUT PARAMETER
C
C   INEST() =   VECTOR CONTAINING NESTING ORDER OF THE MODEL PARAMETERS.
C               EACH VECTOR ELEMENT REPRESENTS ONE MODEL PARAMENTER. 
C               THE PARAMETER CORRESPONDING TO EACH ELEMENT MAY BE FOUND
C               IN THE INPUT FILE.
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER PARAM_NUM(N_VARIABLES), INEST(N_VARIABLES), I
C
        COMMON /I_NEST/ PARAM_NUM, INEST
C
C***********************************************************************
C
        OPEN(UNIT=2,FILE='parameter_nesting.input',
     &          STATUS='UNKNOWN')
C
        READ(2,10)
        DO 5 I=1,N_VARIABLES
            READ(2,20) PARAM_NUM(I),INEST(I)
5       CONTINUE
C
        CLOSE(2)
C
C***********************************************************************
C
10      FORMAT(12(/))
20      FORMAT(I3,19X,I3)
30      FORMAT(/,/)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*** This subroutine reads data from a list provided in the table    ***
C***              'parameter_value_table.input'                      ***
C***    Calling routines:
C***    called subroutines : none
C***********************************************************************
C
        SUBROUTINE READ_TABLE(VECTOR,NUMBER,I)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   VECTOR() = THE INPUT VECTOR TO BE READ
C   NUMBER   = NUMBER OF ELEMENTS IN VECTOR()
C   TABLE_OPEN = .TRUE. -- INPUT FILE IS OPEN
C                .FALSE. -- INPUT TABLE IS NOT OPEN
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER NUMBER,I
        REAL VECTOR(I_VECT_LEN)
        LOGICAL TABLE_OPEN
C
        COMMON /L_TABLE/ TABLE_OPEN

         CHARACTER*82 stars, nums, fhed, pnum, pp, dash
        CHARACTER*82 phed(100), pnum1, pp1
C
C***********************************************************************
        stars='****************************************'//
     &        '******************************************'

        nums ='1234567890123456789012345678901234567890'//
     &        '123456789012345678901234567890123456789012'

        fhed ='**** This file provides a parameter list'//
     &        ' for given inputs of the MIMICS model ****'

        pnum ='--- number of parameter values = XX ----'//
     &        '------------------------------------------'

        pnum1='--- number of parameter values =  1 ----'//
     &        '--- Change as needed ---------------------'

        pp   ='---   Parameter = XXXX.xxxx ------------'//
     &        '------------------------------------------'

        pp1  ='---   Parameter =           ------------'//
     &        '--- Replicate this line as needed --------'

        dash ='----------------------------------------'//
     &        '------------------------------------------'
        
        phed(1)='----------  frequency of operation data '//
     &          '(GHz) ------------------------------------'
        phed(2)='---------Angle of Measurement data (degr'//
     &          'ees)--------------------------------------'


        phed(3)='----------  temperature of soil (C)    -'//
     &          '------------------------------------------'
        phed(4)='----------  temperature of standing wate'//
     &          'r (C)    ---------------------------------'
        phed(5)='----------  temperature of vegetation (C'//
     &          ')   --------------------------------------'


        phed(6)='------ soil volumetric moisture content '//
     &          '- (waterwght in sample)/(volume of smpl)--'
        phed(7)='-------- rms height (centimeters) ------'//
     &          '------------------------------------------'
        phed(8)='------ Gaussian Autocorrelation Paramete'//
     &          'rs - Correlation length (centimeters) ----'
        phed(9)='-------  salinity of standing water (par'//
     &          'ts per thousand (ppt)) -------------------'
        phed(10)='------ snow layer depth (thickness) -- ('//
     &           'meters) ----------------------------------'


        phed(11)='------ trunk gravimetric moisture conten'//
     &           't (0.0 to 1.0) ---------------------------'
        phed(12)='------ dry density of trunk material  (0'//
     &           '.0 to 1.0)     ---------------------------'
        phed(13)='--------  Canopy density (trees per squa'//
     &           're meter)  -------------------------------'
        phed(14)='--------  crown thickness (meters)   ---'//
     &           '------------------------------------------'
        phed(15)='--  trunk diameter data (centimeters) (a'//
     &           'lso is the mean value MU for Normal Dist) '
        phed(16)='----   trunk length data (meters) (also '//
     &           'is the mean value MU for Normal Dist.)  --'
        phed(17)='--  trunk diameter standard deviation (s'//
     &           'igma) data (for Normal Dist) (centimeters)'
        phed(18)='--  trunk length standard deviation (sig'//
     &           'ma) data (for Normal Dist.) data (meters) '


        phed(19)='------ leaf gravimetric moisture content'//
     &           ' (0.0 to 1.0) (wetwght-drywght)/(wetwght)-'
        phed(20)='------ dry density of leaf material  (0.'//
     &           '0 to 1.0)      ---------------------------'
        phed(23)='------ leaf thickness (centimeters) ----'//
     &           '------------------------------------------'
        phed(22)='------ leaf diameter (centimeters) -----'//
     &           '------------------------------------------'
        phed(21)='------ leaf density in crown (leaves per'//
     &           ' cubic meter)      -----------------------'


        phed(24)='------ needle gravimetric moisture conte'//
     &           'nt (0.0 to 1.0) (wetwght-drywght)/(wetwg'//
     &           'ht)  -------------------------------------'
        phed(25)='------ dry density of needle material  ('//
     &           '0.0 to 1.0)     --------------------------'
        phed(28)='------ needle length  (centimeters) ----'//
     &           '------------------------------------------'
        phed(27)='------ needle diameter (centimeters) ---'//
     &           '------------------------------------------'
        phed(26)='------ needle density in crown (needles '//
     &          'per cubic meter)      --------------------'


        phed(29)='------ primary-branch gravimetric moistu'//
     &          're content (0.0 to 1.0) -------------------'
        phed(30)='------ dry density of primary-branch mat'//
     &          'erial  (0.0 to 1.0)     -------------------'
        phed(31)='------ primary-branch density in crown ('//
     &          'branches per cubic meter)      ------------'
        phed(32)='------ primary-branch diameter (centimet'//
     &          'ers) --------------------------------------'
        phed(33)='------ primary-branch length  (meters) -'//
     &          '-------------------------------------------'


        phed(34)='------ secondary-branch gravimetric mois'//
     &          'ture content (0.0 to 1.0) -----------------'
        phed(35)='------ dry density of secondary-branch m'//
     &          'aterial  (0.0 to 1.0)     -----------------'
        phed(36)='------ secondary-branch density in crown'//
     &          ' (branches per cubic meter)      ----------'
        phed(37)='------ secondary-branch diameter (centim'//
     &          'eters) ------------------------------------'
        phed(38)='------ secondary-branch length  (meters)'//
     &          ' ------------------------------------------'

C
C***********************************************************************
C
C       WRITE THE TABLE FILE, WITH EMPTY VALUES, NUMBER=1.
C
        IF(.NOT.TABLE_OPEN)THEN
            OPEN(UNIT=3,FILE='parameter_value_table.input',
     &          STATUS='UNKNOWN')
            WRITE(3,82)stars
 82         FORMAT(A82)
            WRITE(3,82)nums
            WRITE(3,82)stars
            WRITE(3,82)fhed
            WRITE(3,82)stars
            TABLE_OPEN = .TRUE.
        ENDIF
C
        WRITE(3,82)stars
        WRITE(3,82)phed(i)
        WRITE(3,82)pnum
        WRITE(3,82)pnum1
C
        WRITE(3,82)pp
        WRITE(3,82)pp1
        WRITE(3,82)dash

C
C***********************************************************************
C
20      FORMAT(33X,I2)
25      FORMAT(1X)
30      FORMAT(18X,F9.4)
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*************   Subroutine to read histogram tables        ************
C***    Calling routine:      READ_INPUT                             ***
C***    Called subroutines:   none                                   ***
C***    Called functions:     none                                   ***
C***********************************************************************
C
        SUBROUTINE READ_HISTOGRAM(J)
        save
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER J, I, I_TR_DIAM_HIST
        REAL SUM, TRUNK_DIAM_HIST_X(MAXTRDHIST)
        REAL TR_DIAM_HIST_DX(MAXTRDHIST)
        REAL TRUNK_DIAM_HIST_Y(MAXTRDHIST), Y(MAXTRDHIST)
C
        COMMON /TRUNK_HIST/ TRUNK_DIAM_HIST_X, TRUNK_DIAM_HIST_Y,
     &                      TR_DIAM_HIST_DX, I_TR_DIAM_HIST
C
C***********************************************************************
C
        IF(J.EQ.1) THEN
            SUM = 0.0
            OPEN(UNIT=2,FILE='histogram_trunk_diam.input',
     &          STATUS='UNKNOWN')
            READ(2,20) I_TR_DIAM_HIST
            DO 100 I=1,I_TR_DIAM_HIST
                READ(2,30) TRUNK_DIAM_HIST_X(I),Y(I),TR_DIAM_HIST_DX(I)
                SUM = SUM + Y(I)
100         CONTINUE      
            DO 105 I=1,I_TR_DIAM_HIST
                TRUNK_DIAM_HIST_Y(I) = Y(I)/(TR_DIAM_HIST_DX(I)*SUM)
105         CONTINUE
C
        ELSE
          CALL WRITE_ERROR(125)
          STOP
        ENDIF
C
        CLOSE(2)
C
C***********************************************************************
C
20      format(7(/),35x,i2,5(/))
30      format(8x,f5.1,12x,f5.0,6x,f5.1)
C
        RETURN
        END
C
C***********************************************************************

