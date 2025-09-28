C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*** This subroutine initializes fixed constants, flags and counters ***
C*** at the start of the program.                                    ***
c***    Calling routine:      main MIMICS program                    ***
c***    Called subroutines:   none                                   ***
c*********************************************************************** 
C
        SUBROUTINE INITIALIZE(COMPUTE)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C--- PARAMETERS ---
C   N_EPS       = number of elements in the dielectric constant vector
C   N_VARIABLES = number of input variables for the model
C   N_CALLS     = number of subroutine calls that are flagged with
C                   logical variables in the vector CALL_SUB()
C   N_SUB_CALLS = maximum number of subroutine calls within the 
C                   subroutines
C
C--- CONSTANTS ---
C   PI      =   3.14 etc 
C   MU_0    =   permeability of free space  (H/m)
C   EPS_0   =   permittivity of free space  (F/m)
C   LIGHT   =   velocity of light in free space (m/sec)
C
C--- VARIABLES ---
C   OPEN    =   logical variable 
C           =   .TRUE.  if output files are open
C           =   .FALSE. if output files are not open
C
C   COMPUTE =   logical variable 
C           =   .TRUE.  continue computations
C           =   .FALSE. stop looping through the routine
C
C   CALL_SUB()      =   logical vector indicating whether or not to call
C                       a specific subroutine
C                   = .TRUE.  - call the subroutine
C                   = .FALSE. - do not call the subroutine
C
C   STEP_VARIABLE() =   logical vector indicating which variable
C                       is being looped during this sequence
C                       through the program
C                   = .TRUE.  - step this variable this sequence
C                   = .FALSE. - do not step this variable this sequence
C
C   LOOP_COUNT()=   integer vector indicating how many times a variable
C                   has been stepped
C
C   LOOP_NUM()  =   integer vector indicating how many times a variable
C                   will be stepped (based on input data) 
C
C***********************************************************************
C------------------  variable declarations  ----------------------------
C***********************************************************************
C
C----------------------------
c%include 'parameters.include'
        INCLUDE 'parameters.include'
C----------------------------
C
        INTEGER I, LOOP_NUM(N_VARIABLES), LOOP_COUNT(N_VARIABLES)
C
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL CALL_SUB(N_CALLS,N_SUB_CALLS)
        LOGICAL COMPUTE
C
        COMMON /I_COUNT/ LOOP_NUM, LOOP_COUNT
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
C***********************************************************************
C
        OPEN = .FALSE.
        COMPUTE = .TRUE.
C
        DO 10 I=1,N_VARIABLES
            LOOP_COUNT(I) = 1
            STEP_VARIABLE(I) = .FALSE.
            STEP_THIS_TIME(I) = .FALSE.
10      CONTINUE
C
        RETURN
        END
C
C***********************************************************************
