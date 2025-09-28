C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  This subroutine computes the phase and average M- matrices of  ***
C***  ONE needle.                                                    ***
C*********************************************************************** 
C***    Calling subroutine:   CROWN_LAYER                            ***
C***    Called subroutines:   NEEDLE_RLGH_ELPSD_SCAT_MAT             ***
C***                          STOKES_SUB                             ***
C***                          MAT_EIGEN_SUB                          ***
C***                          EXP_MAT                                ***
C***    Called functions:     PDF_NDL                                ***
C*********************************************************************** 
C
        SUBROUTINE NEEDLE_PHASE(P_NEEDLE)
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
C   P_NEEDLE(4,4,4) = Four 4x4 Phase Matrices of the needles
C   P_NEEDLE(L,K,1) = Phase Matrix of ground-direct term
C   P_NEEDLE(L,K,2) = Phase Matrix of crown-ground term
C   P_NEEDLE(L,K,3) = Phase Matrix of ground-crown term
C   P_NEEDLE(L,K,4) = Phase Matrix of direct-crown term
C
C   SMAT_AVG(2,2,2) = The two 2x2 average scattering Matrices of the
C                      needles
C   SMAT_AVG(L,K,1) = average scattering matrix in positive-going
C                      needles
C   SMAT_AVG(L,K,2) = average scattering matrix in negative-going
C                      needles
C
C   M_NEEDLE_P(2,2) = 2x2 complex M matrix in positive-going direction
C                      -- AVERAGE SCATTEING MATRIX TIMES 2*PI*N/k0
C   M_NEEDLE_M(2,2) = 2x2 complex M matrix in negative-going direction
C
C   LAMBDA_NDL_P(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the positive-going M-matrix
C   LAMBDA_NDL_M(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the negative-going M-matrix
C
C   Q_NDL_P(4,4)    = 4x4 COMPLEX EIGENMATRIX of the positive-going 
C                      M-matrix
C   Q_NDL_M(4,4)    = 4x4 COMPLEX EIGENMATRIX of the negative-going 
C                      M-matrix
C   Q_NDL_P_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_NDL_P(4,4)
C   Q_NDL_M_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_NDL_M(4,4)
C
C   EXP_KAPPA_NDL_P(4,4) = positive-going attenuation of the needles
C                         alone
C   EXP_KAPPA_NDL_M(4,4) = negative-going attenuation of the needles
C                         alone
C
C   Z     = EXTENT OVER WHICH TO EVALUATE THE EXPONENTAL OF THE 
C           EXTINCTION MATRIX; Z = dm/COS(THETA)
C
C   DELTA = DIFFERENTIAL ANGULAR ELEMENT FOR INTEGRATION
C   NORM  = NORMALIZATION FACTOR FOR NORMALIZING THE INTEGRATION
C           OVER THE P.D.F
C
C   THETAI,PHII = ANGLES OF INCIDENCE ON THE NEEDLE (RADIANS)
C   THETAS,PHIS = ANGLES OF SCATTERING FROM THE NEEDLE (RADIANS)
C   THETAN,PHIN = ANGLES OF ORIENTATION OF THE NEEDLE (RADIANS)
C
C   PDF_NDL(THETAN) = real-valued P.D.F. function describing the
C                     variation of orientation of the needles with
C                     respect to THETAN
C   PDF  = variable set to PDF 
C        = PDF_NDL(THETAN)
C
C   CROWN_HGHT = height of crown layer (meters)
C   NDL_DIAM   = diameter of needles (centimeters)
C   NDL_LNG    = length of needles (centimeters)
C   THETA  = radar look angle (radians)
C   CTHETA = cosine of radar look angle
C
C   EPSILONRC(3) = EPSR  (complex)
C                = Dielectric constant of needles (e' + j*e'')
C
C   STOKES_MAT(4,4) = 4x4 real Stokes matrix for needles
C   SMAT(2,2)       = 2x2 complex scattering matrix of needles
C
C***********************************************************************
C------------------  variable declations    ----------------------------
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER II,I,L,K,JJ
                               
        REAL LENGTHCM, diamcm, K0A
        REAL WORK, NORM ,PDF, PDF_NDL
        REAL THETAI,PHII,THETAS,PHIS,THETAN,PHIN
        REAL P_NEEDLE(4,4,4),STOKES_MAT(4,4)
        REAL P_NEEDLE_tmp_1(4,4,4), P_NEEDLE_tmp_2(4,4,4)
        REAL P_NEEDLE_tmp_3(4,4,4)
C
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL EXP_KAPPA_NDL_P(4,4),EXP_KAPPA_NDL_M(4,4)
        real lambdacm, lengthresh
        logical rayleigh
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMPLEX J, CWORK, EPSR
        COMPLEX SMAT_AVG(2,2,2),M_NEEDLE_P(2,2),M_NEEDLE_M(2,2)
        COMPLEX SMAT_AVG_tmp_1(2,2,2), SMAT_AVG_tmp_2(2,2,2)
        COMPLEX SMAT_AVG_tmp_3(2,2,2)
        COMPLEX Smat(2,2)
c
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
        real t_rad_start_1,t_rad_stop_1,delta_t_rad_1
        real t_rad_start_2,t_rad_stop_2,delta_t_rad_2
        real t_rad_start_3,t_rad_stop_3,delta_t_rad_3
        real p_rad_start_1,p_rad_stop_1,delta_p_rad_1
        real p_rad_start_2,p_rad_stop_2,delta_p_rad_2
        real p_rad_start_3,p_rad_stop_3,delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
C
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
        common /ndldatar/ diamcm,lengthcm,k0A
        common /ndldatac/ epsr
        common /ndlscatmat/ Smat
C
C--- COMMON BLOCKS FOR EXTINCTION AND PHASE MATRICES ----
C
        COMMON /NDL_EXT/ EXP_KAPPA_NDL_p, EXP_KAPPA_NDL_m
        COMMON /NDL_M_MAT/ M_NEEDLE_P, M_NEEDLE_M
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
        EPSR = EPSILONRC(3)
        diamcm = NDL_DIAM
        lengthcm = NDL_LNG
        k0a = k0
        lambdacm = 100.*(2.*pi/k0)
        lengthresh = 1.5*lambdacm
        if(lengthcm.lt.lengthresh)then
           rayleigh = .true.
        else
           rayleigh = .false.
        endif
        print*,' needle rayleigh = ',rayleigh
C
C***********************************************************************
C   COMPUTE PHASE AND EXTINCTION MATRICES FOR THE NEEDLES
C***********************************************************************
C
C--- INITIALIZE TO ZERO ---
C
        DO 50 i=1,4
         DO 50 k=1,4
          DO 50 l=1,4
             P_NEEDLE(l,k,i) = 0.0
             P_NEEDLE_tmp_1(l,k,i) = 0.0
             P_NEEDLE_tmp_2(l,k,i) = 0.0
             P_NEEDLE_tmp_3(l,k,i) = 0.0
50      CONTINUE
C
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
c
        call pdf_ndl_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
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
c------------ integrate over first interval of THETAn ------------------
c-----------------------------------------------------------------------
c
      IF((n_theta_1.ne.0).and.(n_phi_1.ne.0))THEN
c
      do 150 jj=1,n_theta_1
       thetan = t_rad_start_1 + delta_t_rad_1*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over first interval of PHIn --------------------
c-----------------------------------------------------------------------
c
        DO 145 II = 1,n_phi_1
         PHIN = p_rad_start_1 + delta_p_rad_1*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,1) =
     &               P_NEEDLE_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)

        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,1) =
     &               P_NEEDLE_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,2) =
     &       P_NEEDLE_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,2) =
     &       P_NEEDLE_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,3) =
     &             P_NEEDLE_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,3) =
     &             P_NEEDLE_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,4) =
     &            P_NEEDLE_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,4) =
     &            P_NEEDLE_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,1) =
     &        Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,1) =
     &        Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
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
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over second interval of THETAn -----------------
c-----------------------------------------------------------------------
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
       thetan = t_rad_start_2 + delta_t_rad_2*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over second interval of PHIn -------------------
c-----------------------------------------------------------------------
c
        DO 165 II = 1,n_phi_2
         phin = p_rad_start_2 + delta_p_rad_2*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,1) =
     &               P_NEEDLE_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,1) =
     &               P_NEEDLE_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,2) =
     &       P_NEEDLE_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,2) =
     &       P_NEEDLE_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,3) =
     &             P_NEEDLE_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,3) =
     &             P_NEEDLE_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,4) =
     &            P_NEEDLE_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,4) =
     &            P_NEEDLE_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,1) =
     &        Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,1) =
     &        Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
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
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over third interval of THETAn ------------------
c-----------------------------------------------------------------------
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 200 jj=1,n_theta_3
       thetan = t_rad_start_3 + delta_t_rad_3*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over third interval of PHIn --------------------
c-----------------------------------------------------------------------
c
        DO 195 II = 1,n_phi_3
         phin = p_rad_start_3 + delta_p_rad_3*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,1) =
     &               P_NEEDLE_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,1) =
     &               P_NEEDLE_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,2) =
     &       P_NEEDLE_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,2) =
     &       P_NEEDLE_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
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
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,3) =
     &             P_NEEDLE_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,3) =
     &             P_NEEDLE_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
C
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,4) =
     &            P_NEEDLE_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,4) =
     &            P_NEEDLE_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,1) =
     &        Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,1) =
     &        Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO
C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
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
           P_NEEDLE(l,k,i) = WORK*(delta_p1_t1*P_NEEDLE_tmp_1(l,k,i) +
     &                             delta_p2_t2*P_NEEDLE_tmp_2(l,k,i) +
     &                             delta_p3_t3*P_NEEDLE_tmp_3(l,k,i))
250     CONTINUE
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO 260 k=1,2
         DO 260 l=1,2
           M_NEEDLE_P(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
260     CONTINUE
        DO 270 k=1,2
         DO 270 l=1,2
           M_NEEDLE_M(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))
270     CONTINUE
C
        RETURN
        END
C
C***********************************************************************
