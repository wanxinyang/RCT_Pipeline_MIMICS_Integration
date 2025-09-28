C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  This subroutine computes the phase and average M- matrices of  ***
C***  ONE branch.        (3rd branch)                                ***
C                                                                    ***
C     modified to force reciprocity  for resonant branches (3-7-91)  ***
C     (really done by Leland Pierce on 3-10-92,                      ***
C       as a fix to version 1.2)                                     ***
C                                                                    ***
C*********************************************************************** 
C***    Calling subroutine:   CROWN_LAYER                            ***
C***    Called subroutines:   RESONANT_CYL_SCAT_MAT                  ***
C***                          long_thin_cyl_scatmat                  ***
C***                          STOKES_SUB                             ***
C***                          MAT_EIGEN_SUB                          ***
C***                          EXP_MAT                                ***
C***                          PDF_BR3_SETUP                          ***
C***    Called functions:     PDF_BR3                                ***
C*********************************************************************** 
C
        SUBROUTINE BRANCH_PHASE3(P_BRANCH_3)
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
C   P_BRANCH_3(4,4,4) = Four 4x4 Phase Matrices of the branches
C   P_BRANCH_3(L,K,1) = Phase Matrix of ground-direct term
C   P_BRANCH_3(L,K,2) = Phase Matrix of crown-ground term
C   P_BRANCH_3(L,K,3) = Phase Matrix of ground-crown term
C   P_BRANCH_3(L,K,4) = Phase Matrix of direct-crown term
C
C   SMAT_AVG(2,2,2) = The two 2x2 average scattering Matrices of the
C                      branches
C   SMAT_AVG(L,K,1) = average scattering matrix in positive-going
C                      direction
C   SMAT_AVG(L,K,2) = average scattering matrix in negative-going
C                      direction
C
C   M_BRANCH_P3(2,2) = 2x2 complex M matrix in positive-going direction
C                      -- AVERAGE SCATTEING MATRIX TIMES 2*PI*N/k0
C   M_BRANCH_M3(2,2) = 2x2 complex M matrix in negative-going direction
C
C   LAMBDA_BR_P(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the positive-going M-matrix
C   LAMBDA_BR_M(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the negative-going M-matrix
C
C   Q_BR_P(4,4)    = 4x4 COMPLEX EIGENMATRIX of the positive-going 
C                      M-matrix
C   Q_BR_M(4,4)    = 4x4 COMPLEX EIGENMATRIX of the negative-going 
C                      M-matrix
C   Q_BR_P_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_BR_P(4,4)
C   Q_BR_M_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_BR_M(4,4)
C
C   EXP_KAPPA_BR_P3(4,4) = positive-going attenuation of the branches
C                         alone
C   EXP_KAPPA_BR_M3(4,4) = negative-going attenuation of the branches
C                         alone
C
C   Z     = EXTENT OVER WHICH TO EVALUATE THE EXPONENTAL OF THE 
C           EXTINCTION MATRIX; Z = CROWN_HGHT/COS(THETA)
C
C   DELTA = DIFFERENTIAL ANGULAR ELEMENT FOR INTEGRATION
C   NORM  = NORMALIZATION FACTOR FOR NORMALIZING THE INTEGRATION
C           OVER THE P.D.F
C
C   THETAI,PHII = ANGLES OF INCIDENCE ON THE BRANCH (RADIANS)
C   THETAS,PHIS = ANGLES OF SCATTERING FROM THE BRANCH (RADIANS)
C   THETAC,PHIC = ANGLES OF ORIENTATION OF THE BRANCH (RADIANS)
C
C   LFLAG = LOGICAL VARIABLE RETURNED FROM THE SUBROUTINE
C            LONG_CYL_SCAT_MAT
C         = .TRUE. - SCATTERING MARTIX COMPUTED NORMALLY
C         = .FALSE.- SCATTERING MARTIX NOT COMPUTED.
C                       -- elements set to zero for end-on incidence
C                          or for values of sin(x)/(x) < thresh
C
C   PDF_BR3(THETAC) = real-valued P.D.F. function describing the
C                    variation of orientation of the branches with
C                    respect to THETAC
C   PDF  = variable set to PDF 
C        = PDF_BR3(THETAC)
C
C   CROWN_HGHT  = height of crown layer (meters)
C   BR1_DIAM  = diamcm 
C             = diameter of branches (centimeters)
C   BR1_LNG   = lengthm 
C             = length of branches (meters)
C   THETA  = radar look angle (radians)
C   CTHETA = cosine of radar look angle
C
C   EPSILONRC(7) = EPSR  (complex)
C                = Dielectric constant of branches (e' + j*e'')
C
C   STOKES_MAT(4,4) = 4x4 real Stokes matrix for branches
C   SMAT(2,2)       = 2x2 complex scattering matrix of branches
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)

        COMMON /R_CONSTANTS/ PI,MU_0,EPS_0,LIGHT
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA

        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC

        INTEGER II,I,L,K,JJ
C
        REAL diamcm, lengthm, k0a
        REAL THETAI,PHII,THETAS,PHIS,THETAC,PHIC
        REAL P_BRANCH_3(4,4,4),STOKES_MAT(4,4)
        REAL P_BRANCH_tmp_1(4,4,4), P_BRANCH_tmp_2(4,4,4)
        REAL P_BRANCH_tmp_3(4,4,4)
        REAL EXP_KAPPA_BR_P3(4,4),EXP_KAPPA_BR_M3(4,4)
        REAL PDF, PDF_BR3, NORM, WORK
C
        COMPLEX J, CWORK
        COMPLEX SMAT_AVG(2,2,2),M_BRANCH_P3(2,2),M_BRANCH_M3(2,2)
        COMPLEX SMAT_AVG_tmp_1(2,2,2), SMAT_AVG_tmp_2(2,2,2)
        COMPLEX SMAT_AVG_tmp_3(2,2,2), smat_tmp1(2,2), smat_tmp2(2,2)
        COMPLEX epsr
        COMPLEX Smat(2,2)

        real lengthcm, lengthresh, diamthresh, lambdacm
        logical resonant
C   
        LOGICAL LFLAG
C
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
c
        real t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        real t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        real t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        real p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        real p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        real p_rad_start_3, p_rad_stop_3, delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
C
C--- COMMON BLOCKS PASSED TO RESONANT_CYL_SCAT_MAT SUBROUTINE ---
C
        common /rcyldata/ diamcm,lengthm,k0a
        common /ccyldata/ epsr
        common /cylscatmat/ Smat
C
C--- COMMON BLOCKS FOR EXTINCTION AND PHASE MATRICES ----
C
        COMMON /BR3_EXT/ EXP_KAPPA_BR_P3, EXP_KAPPA_BR_M3
        COMMON /BR3_M_MAT/ M_BRANCH_P3, M_BRANCH_M3
C
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C------------------  initialize variables   ----------------------------
C***********************************************************************
C
        J = CMPLX(0.0,1.0)
        EPSR = EPSILONRC(9)
        diamcm = BR3_DIAM
        lengthm = BR3_LNG
        k0a = k0

        lengthcm = 100.*lengthm
        lambdacm = 100.*(2.*pi/k0)
        lengthresh = 2.0*lambdacm
        diamthresh = 0.1*lambdacm/(2.*pi*cabs(csqrt(epsr)))

        if((diamcm.le.diamthresh).and.(lengthcm.ge.lengthresh))then
            resonant = .false.
        else
            resonant = .true.
        endif

        print*,' 3rd branch, resonant = ',resonant
C
C***********************************************************************
C   COMPUTE PHASE AND EXTINCTION MATRICES FOR THE BRANCHES
C***********************************************************************
C
C--- INITIALIZE TO ZERO ---
C
        DO 50 i=1,4
         DO 50 k=1,4
          DO 50 l=1,4
             P_branch_3(l,k,i) = 0.0
             P_branch_tmp_1(l,k,i) = 0.0
             P_branch_tmp_2(l,k,i) = 0.0
             P_branch_tmp_3(l,k,i) = 0.0
50      CONTINUE
c
        DO 60 i=1,2
         DO 60 k=1,2
          DO 60 l=1,2
             Smat_avg(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_1(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_2(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_3(l,k,i) = CMPLX(0.0,0.0)
60      CONTINUE
C
C***********************************************************************
C--- INTEGRATE OVER ORIENTATIONS --- ASSUME AZIMUTHAL SYMMETRY ---------
C***********************************************************************
C
        NORM = 2.0*PI
        PHII = 0.0
C
        call pdf_br3_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
     &       delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,
     &       delta_t_rad_2,n_theta_3,t_rad_start_3,t_rad_stop_3,
     &       delta_t_rad_3,n_phi_1,p_rad_start_1,p_rad_stop_1,
     &  delta_p_rad_1,n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &       n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        delta_p1_t1 = delta_p_rad_1*delta_t_rad_1
        delta_p2_t2 = delta_p_rad_2*delta_t_rad_2
        delta_p3_t3 = delta_p_rad_3*delta_t_rad_3
c
c-----------------------------------------------------------------------
c------------ integrate over first interval of THETAc ------------------
c-----------------------------------------------------------------------
c
      IF((n_theta_1.ne.0).and.(n_phi_1.ne.0))THEN
c
      do 150 jj=1,n_theta_1
       thetac = t_rad_start_1 + delta_t_rad_1*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3a thetac,pdf = ',thetac,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over first interval of PHIc --------------------
c-----------------------------------------------------------------------
c
        DO 145 II = 1,n_phi_1
         PHIC = p_rad_start_1 + delta_p_rad_1*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,1) =
     &          P_BRANCH_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,1) =
     &          P_BRANCH_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) = smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) = smat(2,2)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,2) =
     &          P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,2) =
     &          P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
c-----------------------------------------------------------------------
c--- process specular crown scatter for resonant cylinder solution ---
c        forcing reciprocity to hold
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,2) =
     &         P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,2) =
     &         P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c-----------------------------------------------------------------------
C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,4) =
     &         P_BRANCH_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,4) =
     &         P_BRANCH_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,1) =
     &           Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,1) =
     &           Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,2) =
     &          Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,2) =
     &          Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO

C
145   CONTINUE
C
      ENDIF
150   CONTINUE
c
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over second interval of THETAc -----------------
c-----------------------------------------------------------------------
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
       thetac = t_rad_start_2 + delta_t_rad_2*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3b thetac,pdf = ',thetac,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over second interval of PHIc -------------------
c-----------------------------------------------------------------------
c
        DO 165 II = 1,n_phi_2
         PHIC = p_rad_start_2 + delta_p_rad_2*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,1) =
     &          P_BRANCH_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,1) =
     &          P_BRANCH_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) =  smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) =  smat(2,2)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,2) =
     &          P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,2) =
     &          P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
c-----------------------------------------------------------------------
c--- process specular crown scatter for resonant cylinder solution ---
c        forcing reciprocity to hold
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,2) =
     &         P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,2) =
     &         P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c-----------------------------------------------------------------------
C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,4) =
     &         P_BRANCH_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,4) =
     &         P_BRANCH_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,1) =
     &           Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,1) =
     &           Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,2) =
     &          Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,2) =
     &          Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO
C
165   CONTINUE
C
      ENDIF
170   CONTINUE
c
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over third interval of THETAc -----------------
c-----------------------------------------------------------------------
c
      IF((n_theta_3.ne.0).and.(n_phi_3.ne.0))THEN
c
      do 200 jj=1,n_theta_3
       thetac = t_rad_start_3 + delta_t_rad_3*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3c thetac,pdf = ',thetac,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over third interval of PHIc -------------------
c-----------------------------------------------------------------------
c
        DO 195 II = 1,n_phi_3
         PHIC = p_rad_start_3 + delta_p_rad_3*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,1) =
     &          P_BRANCH_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,1) =
     &          P_BRANCH_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) =  smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,2) =
     &          P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,2) =
     &          P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
c-----------------------------------------------------------------------
c--- process specular crown scatter for resonant cylinder solution ---
c        forcing reciprocity to hold
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,2) =
     &         P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,2) =
     &         P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c-----------------------------------------------------------------------
C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,4) =
     &         P_BRANCH_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,4) =
     &         P_BRANCH_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,1) =
     &           Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,1) =
     &           Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,2) =
     &          Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,2) =
     &          Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
          ENDDO
         ENDDO
C
195   CONTINUE
C
      ENDIF
200   CONTINUE
c
      endif
C
C***********************************************************************
C--- CORRECT TO ACTUAL PHASE AND M - MATRICES --------------------------
C***********************************************************************
C
        WORK = 1.0/NORM
C
        DO 250 i=1,4
         DO 250 k=1,4
          DO 250 l=1,4
           P_branch_3(l,k,i) = WORK*(delta_p1_t1*P_branch_tmp_1(l,k,i) +
     &                               delta_p2_t2*P_branch_tmp_2(l,k,i) +
     &                               delta_p3_t3*P_branch_tmp_3(l,k,i))
250     CONTINUE
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO 260 k=1,2
         DO 260 l=1,2
          M_BRANCH_P3(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
260     CONTINUE
        DO 270 k=1,2
         DO 270 l=1,2
          M_BRANCH_M3(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))

270     CONTINUE
c-----------------------------------------------------------------------
c   correct for reciprocity correction
c-----------------------------------------------------------------------
        if(resonant)then
         DO k=1,4
          DO l=1,4
           P_branch_3(l,k,2) = 0.5*P_branch_3(l,k,2)
           P_branch_3(l,k,3) = 0.5*P_branch_3(l,k,3)
          enddo
         enddo
        endif
c-----------------------------------------------------------------------
C
        RETURN
        END
C
C***********************************************************************
