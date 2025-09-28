C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  This subroutine computes the scattering matrix of a single     ***
C***  leaf using a physical optics approximation.                    ***
C*********************************************************************** 
C***    Calling subroutine:   CROWN_PHASE                            ***
C***    Called subroutines:   none                                   ***
C*********************************************************************** 
C
        SUBROUTINE LEAF_PHYS_OPTICS_SCAT_MAT(THETAi,PHIi,THETAs,PHIs,
     &             THETAd,PHId)
        save
C
C***********************************************************************
C------------------  variable definitions   ----------------------------
C***********************************************************************
C
C   J   = SQRT(-1)
C   PI  = 3.1415 etc.
C   k0  = k0a
C       = WAVE NUMBER (1/M)
C
C   THETAI,PHII = ANGLES OF INCIDENCE ON THE LEAF (RADIANS)
C   THETAS,PHIS = ANGLES OF SCATTERING FROM THE LEAF (RADIANS)
C   THETAd,PHId = ANGLES OF ORIENTATION OF THE LEAF (RADIANS) - input
C   THETAJ,PHIJ = ANGLES OF ORIENTATION OF THE LEAF (RADIANS) - working
C
C   TOL =   tolerance for checking sin(x)/(x) = 1 for small x
C
C   LEAF_DIAM = diameter of leaf (centimeters)
C   LEAF_TAU  = thickness of leaf (centimeters)
C   TAU     = thickness of leaf (meters)
C   A = B   = square dimensions of the leaf 
C             (The square with the same area as the disk.)
C
C   EPSILONRC(4) = EPSR  (complex)
C                = Dielectric constant of leaves (e' + j*e'')
C
C   SMAT(2,2)    = 2x2 complex scattering matrix of branches
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        INCLUDE 'parameters.include'
C----------------------------
C
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL FREQ_HERTZ, WAVELENGTH, k0
C
        REAL TAU,Q, P, A, B, TOL
        REAL COS_PHI_1
        REAL THETAi, PHIi, THETAs, PHIs, THETAj, PHIj, THETAd, PHId
        REAL SIN_THETAi, COS_THETAi, SIN_THETAj, COS_THETAj
        REAL SIN_THETAs, COS_THETAs, COS_PHIj, SIN_PHIij, SIN_PHIji 
        REAL COS_PHIij, COS_PHIji, SIN_PHIsj, SIN_PHIjs, COS_PHIsj 
        REAL COS_BETA, SIN_BETA, SIN_PHI, COS_PHI, COS_BETA_COS_PHI
        REAL SIN_PHIPR, COS_PHIPR, SIN_BETAPR
        REAL U, V, SIN_U, SIN_V, SIN_U_OVER_U, SIN_V_OVER_V
        REAL COS_B_COS_P, SIN_Tj_SIN_Ti, SIN_Ts_SIN_Tj 
        REAL COS_Tj_SIN_Pjs, COS_Ts_SIN_Psj
        REAL SIN_Pij_COS_Psj, C_Pij_C_Ti_C_Tj, C_Psj_C_Ts_C_Tj
        REAL SAVE1, SAVE2, SAVE3, SAVE4, SAVE5, WORK1, WORK2
C
        COMPLEX J, EPSILONR(N_EPS), EPSILONRC(N_EPS), EPSR, RES, CONST
        COMPLEX GAMMA_H, GAMMA_E, GAMMA_HE_C_1, GAMMA_HE_C_2  
        COMPLEX SMAT(2,2)
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /C_DIELECTRIC/ EPSILONR, EPSILONRC
        COMMON /LEAF_SCAT_MAT/ Smat
C
C----------------------------
c%include 'constants.include'
        INCLUDE 'constants.include'
C----------------------------
C***********************************************************************
C
        TOL = 0.0001
        J = CMPLX(0.0,1.0)
        TAU = LEAF_TAU/100.
        epsr = EPSILONRC(4)
        A = 0.5*sqrt(pi)*LEAF_DIAM/100.
        B = A 
C
C***********************************************************************
C
C--- CHECK ORIENTATION ---
C
        COS_PHI_1=-(SIN(THETAd)*SIN(THETAi)*COS(PHId-PHIi)+
     &           COS(THETAi)*COS(THETAd))
C
        IF(COS_PHI_1 .LT. 0.0)THEN
          THETAj = PI - THETAd
          PHIj = PI + PHId
          COS_PHI_1 = -COS_PHI_1
        ELSE
            THETAj = THETAd
            PHIj = PHId
        ENDIF
C
C--- EVALUATE CONSTANTS ---
C
        SIN_THETAi = SIN(THETAi)
        COS_THETAi = COS(THETAi)
        SIN_THETAj = SIN(THETAj)
        COS_THETAj = COS(THETAj)
        SIN_THETAs = SIN(THETAs)
        COS_THETAs = COS(THETAs)
C
        COS_PHIj = COS(PHIj)
C
        SIN_PHIij  = SIN(PHIi-PHIj)
        SIN_PHIji  = -SIN_PHIij
        COS_PHIij  = COS(PHIi-PHIj)
        COS_PHIji  = COS_PHIij
C
        SIN_PHIsj  = SIN(PHIs-PHIj)
        SIN_PHIjs  = -SIN_PHIsj
        COS_PHIsj  = COS(PHIs-PHIj)
C
C---
C
        WORK1 = SIN_THETAi*SIN_PHIji
C
        Q = 1./SQRT(1.- WORK1*WORK1)
        COS_BETA = Q*COS_PHI_1
        SIN_BETA = Q*(-COS_THETAj*SIN_THETAi*COS_PHIji+
     &           COS_THETAi*SIN_THETAj)
        SIN_PHI = SIN_THETAi*SIN_PHIij
        COS_PHI = SQRT(1.-SIN_PHI*SIN_PHI)
C
        COS_BETA_COS_PHI = COS_BETA*COS_PHI
C
        SIN_PHIPR  = SIN_THETAs*SIN_PHIsj
        COS_PHIPR  = SQRT(1.0-SIN_PHIPR*SIN_PHIPR)
        SIN_BETAPR = (COS_THETAs*SIN_THETAj -
     &                COS_THETAj*SIN_THETAs*COS_PHIsj)/COS_PHIPR
C
c        P = 1./SQRT(1.-COS_PHI_1*COS_PHI_1)     {different relationship than in the report
        P= 1.0/SQRT(1-COS_BETA_COS_PHI*COS_BETA_COS_PHI)
C
        RES = J/(k0*TAU*(epsr - 1.0))
C
        GAMMA_H = 1./(1.+2.*RES/COS_PHI_1)
        GAMMA_E = 1./(1.+2.*RES*COS_PHI_1)
C
        U = 0.5*K0*A*(SIN_PHI - SIN_PHIPR)
        V = 0.5*K0*B*(SIN_BETA*COS_PHI - SIN_BETAPR*COS_PHIPR)
C
        SIN_U = SIN(U)
        SIN_V = SIN(V)
C
        IF(ABS(SIN_U).LE.TOL)THEN
            SIN_U_OVER_U = 1.0
        ELSE
            SIN_U_OVER_U = SIN_U/U
        ENDIF
C
        IF(ABS(SIN_V).LE.TOL)THEN
            SIN_V_OVER_V = 1.0
        ELSE
            SIN_V_OVER_V = SIN_V/V
        ENDIF
C
        CONST = -J*k0*A*B*SIN_V_OVER_V*SIN_U_OVER_U*P*P/(2.0*PI)
C
C***********************************************************************
C   EVALUATE SCATTERING MATRIX
C***********************************************************************
C
C--- EVALUATE CONSTANTS ---
C
        COS_B_COS_P = COS_BETA*COS_PHI
C
        GAMMA_HE_C_1  = (GAMMA_H-GAMMA_E)*COS_B_COS_P
        GAMMA_HE_C_2  = GAMMA_H - COS_B_COS_P*COS_B_COS_P*GAMMA_E
C
        SIN_Tj_SIN_Ti = SIN_THETAj*SIN_THETAi
        SIN_Ts_SIN_Tj = SIN_THETAs*SIN_THETAj
C
        COS_Tj_SIN_Pjs = COS_THETAj*SIN_PHIjs
        COS_Ts_SIN_Psj = COS_THETAs*SIN_PHIsj
C
        SIN_Pij_COS_Psj = SIN_PHIij*COS_PHIsj
C
        C_Pij_C_Ti_C_Tj = COS_PHIij*COS_THETAi*COS_THETAj
        C_Psj_C_Ts_C_Tj = COS_PHIsj*COS_THETAs*COS_THETAj
C
        SAVE1 = SIN_Tj_SIN_Ti + C_Pij_C_Ti_C_Tj
        SAVE2 = SIN_Ts_SIN_Tj + C_Psj_C_Ts_C_Tj 
        SAVE3 = SIN_PHIij*COS_Ts_SIN_Psj
        SAVE4 = SIN_PHIij*SAVE2 
        SAVE5 = SIN_PHIij*COS_Tj_SIN_Pjs
C
C--- Svv TERM ---
C
        WORK1 = SAVE1*SAVE2 + COS_THETAi*SAVE3
        WORK2 = COS_PHIij*SAVE2 + COS_THETAj*SAVE3
C
        Smat(1,1) = CONST*(WORK1*GAMMA_HE_C_1 + WORK2*GAMMA_HE_C_2)
C
C--- Svh TERM ---
C
        WORK1 = -COS_THETAj*SAVE4+COS_PHIij*COS_Ts_SIN_Psj
        WORK2 = -COS_THETAi*SAVE4 + SAVE1*COS_Ts_SIN_Psj
C
        Smat(1,2) = CONST*(WORK1*GAMMA_HE_C_1 + WORK2*GAMMA_HE_C_2)
C
C--- Shv TERM ---
C 
        WORK1 = SAVE1*COS_Tj_SIN_Pjs + COS_THETAi*SIN_Pij_COS_Psj
        WORK2 = COS_PHIij*COS_Tj_SIN_Pjs + COS_THETAj*SIN_Pij_COS_Psj
C
        Smat(2,1) = CONST*(WORK1*GAMMA_HE_C_1 + WORK2*GAMMA_HE_C_2)
C
C--- Shh TERM ---
C
        WORK1 = -COS_PHIj*SAVE5 + COS_PHIij*COS_PHIsj
        WORK2 = -COS_THETAi*SAVE5 + SAVE1*COS_PHIsj
C
        Smat(2,2) = CONST*(WORK1*GAMMA_HE_C_1 + WORK2*GAMMA_HE_C_2)
C
        RETURN
        END
C
C***********************************************************************
