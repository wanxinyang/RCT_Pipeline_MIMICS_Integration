C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  This subroutine computes the phase and average M- matrices of  ***
C***  ONE leaf.                                                      ***
C*********************************************************************** 
C***    Calling subroutine:   CROWN_LAYER                            ***
C***    Called subroutines:   LEAF_PHYS_OPTICS_SCAT_MAT              ***
C***                          STOKES_SUB                             ***
C***                          MAT_EIGEN_SUB                          ***
C***                          EXP_MAT                                ***
C***                          PDF_LF_SETUP                           ***
C***    Called functions:     PDF_LF                                 ***
C*********************************************************************** 
C
        SUBROUTINE LEAF_PHASE(P_LEAF)
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
C   P_LEAF(4,4,4) = Four 4x4 Phase Matrices of the leaf
C   P_LEAF(L,K,1) = Phase Matrix of ground-direct term
C   P_LEAF(L,K,2) = Phase Matrix of crown-ground term
C   P_LEAF(L,K,3) = Phase Matrix of ground-crown term
C   P_LEAF(L,K,4) = Phase Matrix of direct-crown term
C
C   SMAT_AVG(2,2,2) = The two 2x2 average scattering Matrices of the
C                      leaves
C   SMAT_AVG(L,K,1) = average scattering matrix in positive-going
C                      leaves
C   SMAT_AVG(L,K,2) = average scattering matrix in negative-going
C                      leaves
C
C   M_LEAF_P(2,2) = 2x2 complex M matrix in positive-going direction
C                      -- AVERAGE SCATTEING MATRIX TIMES 2*PI/k0
C   M_LEAF_M(2,2) = 2x2 complex M matrix in negative-going direction
C
C   LAMBDA_LF_P(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the positive-going M-matrix
C   LAMBDA_LF_M(4) = 4 ELEMENT COMPLEX VECTOR CONTAINING THE EIGENVALUES 
C                      of the negative-going M-matrix
C
C   Q_LF_P(4,4)    = 4x4 COMPLEX EIGENMATRIX of the positive-going 
C                      M-matrix
C   Q_LF_M(4,4)    = 4x4 COMPLEX EIGENMATRIX of the negative-going 
C                      M-matrix
C   Q_LF_P_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_LF_P(4,4)
C   Q_LF_M_INV(4,4) = 4x4 COMPLEX INVERVSE OF Q_LF_M(4,4)
C
C   EXP_KAPPA_LF_P(4,4) = positive-going attenuation of the leaves
C                         alone
C   EXP_KAPPA_LF_M(4,4) = negative-going attenuation of the leaves
C                         alone
C
C   Z     = EXTENT OVER WHICH TO EVALUATE THE EXPONENTAL OF THE 
C           EXTINCTION MATRIX; Z = CROWN_HGHT/COS(THETA)
C
C   DELTA = DIFFERENTIAL ANGULAR ELEMENT FOR INTEGRATION
C   NORM  = NORMALIZATION FACTOR FOR NORMALIZING THE INTEGRATION
C           OVER THE P.D.F
C
C   THETAI,PHII = ANGLES OF INCIDENCE ON THE LEAF (RADIANS)
C   THETAS,PHIS = ANGLES OF SCATTERING FROM THE LEAF (RADIANS)
C   THETAD,PHID = ANGLES OF ORIENTATION OF THE LEAF (RADIANS)
C
C   PDF_LF(THETAD) = real-valued P.D.F. function describing the
C                    variation of orientation of the leaves with
C                    respect to THETAD
C   PDF  = variable set to PDF 
C        = PDF_LF(THETAD)
C
C   CROWN_HGHT = height of crown layer (meters)
C   LEAF_DIAM  = diameter of leaves (centimeters)
C   LEAF_TAU   = thickness of leaves (centimeters)
C   THETA  = radar look angle (radians)
C   CTHETA = cosine of radar look angle
C
C   EPSILONRC(4) = EPSR  (complex)
C                = Dielectric constant of leaves (e' + j*e'')
C
C   STOKES_MAT(4,4) = 4x4 real Stokes matrix for leaves
C   SMAT(2,2)       = 2x2 complex scattering matrix of leaves
c
c   L_PO = logical variable -- .true.  = call physical optics routine
c                              .false. = call rayleigh routine 
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

        REAL  WORK, NORM ,PDF, PDF_LF
        REAL THETAI,PHII,THETAS,PHIS,THETAD,PHID
        REAL P_LEAF(4,4,4), STOKES_MAT(4,4)
        REAL P_LEAF_TMP_1(4,4,4),P_LEAF_TMP_2(4,4,4),P_LEAF_TMP_3(4,4,4)
        REAL EXP_KAPPA_LF_P(4,4),EXP_KAPPA_LF_M(4,4)
C
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMPLEX J, CWORK
        COMPLEX SMAT_AVG(2,2,2),M_LEAF_P(2,2),M_LEAF_M(2,2)
        COMPLEX SMAT_AVG_TMP_1(2,2,2), SMAT_AVG_TMP_2(2,2,2)
        COMPLEX SMAT_AVG_TMP_3(2,2,2)
        COMPLEX Smat(2,2)
c
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
        real t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        real t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        real t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        real p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        real p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        real p_rad_start_3, p_rad_stop_3, delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
c
        logical L_PO
C
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
        COMMON /LEAF_SCAT_MAT/ Smat
C
C--- COMMON BLOCKS FOR EXTINCTION AND PHASE MATRICES ----
C
        COMMON /LEAF_EXT/ EXP_KAPPA_LF_p, EXP_KAPPA_LF_m
        COMMON /LEAF_M_MAT/ M_LEAF_P, M_LEAF_M
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
c
        if((WAVELENGTH/(LEAF_DIAM/100.)).lt.1.5)then
            L_PO = .true.
        else
            L_PO = .false.
        endif
        print*,' Leaf Phys Optics = ',L_PO
C
C***********************************************************************
C   COMPUTE PHASE AND EXTINCTION MATRICES FOR THE LEAVES
C***********************************************************************
C
C--- INITIALIZE TO ZERO ---
C
        DO 50 i=1,4
         DO 50 k=1,4
          DO 50 l=1,4
             P_LEAF(l,k,i) = 0.0
             P_LEAF_tmp_1(l,k,i) = 0.0
             P_LEAF_tmp_2(l,k,i) = 0.0
             P_LEAF_tmp_3(l,k,i) = 0.0
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
c
        call pdf_lf_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
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
c------------ integrate over first interval of THETAd ------------------
c-----------------------------------------------------------------------
C
      IF((n_theta_1.ne.0).and.(n_phi_1.ne.0))THEN
c
      do 150 jj=1,n_theta_1
        thetad = t_rad_start_1 + delta_t_rad_1*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad, pdf
C
c-----------------------------------------------------------------------
c------------ integrate over first interval of PHId --------------------
c-----------------------------------------------------------------------
c
          do 145 ii=1,n_phi_1
            PHId = p_rad_start_1 + delta_p_rad_1*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c                {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,1) =
     &               P_LEAF_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,1) =
     &               P_LEAF_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c                 {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,2) =
     &             P_LEAF_tmp_1(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,2) =
     &             P_LEAF_tmp_1(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c                     {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,3) =
     &             P_LEAF_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,3) =
     &             P_LEAF_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

C         
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c               { PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF

            CALL STOKES_SUB(SMAT,STOKES_MAT)

            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,4) =
     &             P_LEAF_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)

            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,4) =
     &             P_LEAF_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO


C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,1) =
     &             Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,1) =
     &             Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,2) =
     &            Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,2) =
     &            Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
145       CONTINUE
C
         ENDIF
150   CONTINUE
c
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over second interval of THETAd -----------------
c-----------------------------------------------------------------------
C
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
        thetad = t_rad_start_2 + delta_t_rad_2*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over second interval of PHId -------------------
c-----------------------------------------------------------------------
c
          do 165 ii=1,n_phi_2
            PHId = p_rad_start_2 + delta_p_rad_2*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c                {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,1) =
     &               P_LEAF_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,1) =
     &               P_LEAF_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c                 {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,2) =
     &             P_LEAF_tmp_2(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,2) =
     &             P_LEAF_tmp_2(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c                     {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,3) =
     &             P_LEAF_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,3) =
     &             P_LEAF_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C         
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c               { PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,4) =
     &             P_LEAF_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,4) =
     &             P_LEAF_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,1) =
     &             Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,1) =
     &             Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,2) =
     &            Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,2) =
     &            Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO
C
165       CONTINUE
C
         ENDIF
170   CONTINUE
c
      endif
c
c-----------------------------------------------------------------------
c------------ integrate over third interval of THETAd ------------------
c-----------------------------------------------------------------------
C
      IF((n_theta_3.ne.0).and.(n_phi_3.ne.0))THEN
c
      do 200 jj=1,n_theta_3
        thetad = t_rad_start_3 + delta_t_rad_3*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad,pdf
C
c-----------------------------------------------------------------------
c------------ integrate over third interval of PHId --------------------
c-----------------------------------------------------------------------
c
          do 195 ii=1,n_phi_3
            PHId = p_rad_start_3 + delta_p_rad_3*float(ii-1)
C
C--- PHASE MATRIX CASE I --- BACKSCATTER LOOKING UP FROM GROUND
C                                 (GROUND-DIRECT)
C
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c                {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,1) =
     &               P_LEAF_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,1) =
     &               P_LEAF_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C
C--- PHASE MATRIX CASE II ---  CANOPY SPECULAR DOWNWARD TOWARD GROUND 
C                               (CROWN-GROUND)
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c                 {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,2) =
     &             P_LEAF_tmp_3(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,2) =
     &             P_LEAF_tmp_3(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C
C--- PHASE MATRIX CASE III --- CANOPY SPECULAR UPWARD FROM GROUND 
C                               (GROUND-CROWN)
C
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c                     {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,3) =
     &             P_LEAF_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,3) =
     &             P_LEAF_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
C         
C--- PHASE MATRIX CASE IV --- DIRECT BACKSCATTER FROM CROWN
C                             (DIRECT-CROWN)
C
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c               { PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,4) =
     &             P_LEAF_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,4) =
     &             P_LEAF_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO


C
C--- EXTINCTION MATRIX CASE I --- POSITIVE (UPWARD) PROPAGATING CASE
C
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,1) =
     &             Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,1) =
     &             Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
C--- EXTINCTION MATRIX CASE II --- NEGATIVE (DOWNWARD) PROPAGATING CASE
C
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,2) =
     &            Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,2) =
     &            Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
195       CONTINUE
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
             P_LEAF(l,k,i) = WORK*(delta_p1_t1*P_LEAF_tmp_1(l,k,i) +
     &                             delta_p2_t2*P_LEAF_tmp_2(l,k,i) +
     &                             delta_p3_t3*P_LEAF_tmp_3(l,k,i))
250     CONTINUE
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO 260 k=1,2
         DO 260 l=1,2
             M_LEAF_P(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
260     CONTINUE
        DO 270 k=1,2
         DO 270 l=1,2
             M_LEAF_M(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))
270     CONTINUE
C
        RETURN
        END
C
C***********************************************************************
