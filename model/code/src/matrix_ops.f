c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the exponential of a 4x4 matrix.       ***
C*** i.e. ANSWER = exp(kappa*z) where kappa is th 4x4 matrix.        ***
c*** The matrix KAPPA is supplied in terms of its eigenvalue matrices***
c*** such that ANSWER = Q*D*QINV  is a 4x4 matrix.                   ***
c***    Calling routine:      TRUNK                                  ***
c***    Called subroutines:   CMATMULT3                              ***
c***********************************************************************
C
        SUBROUTINE EXP_MAT_test(LAMBDA,Q,QINV,z,ANSWER)
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
C
C   LAMBDA(4)   =   4 ELEMENT VECTOR CONTAINING THE EIGENVALUES(COMPLEX)
C   Q(4,4)      =   4x4 COMPLEX EIGENMATRIX
C   QINV(4,4)   =   4x4 COMPLEX INVERVSE OF Q(4,4)
C   Z           =   EXTENT OF PROPAGATION THE VOLUME(METERS)
C               =   (H/COS(THETA))
C   D(4,4)      =   COMPLEX DIAGONAL MATRIX
C   ANSWER(4,4) =   COMPLEX OUTPUT MATRIX
C
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
C
        INTEGER I,J
        REAL Z
        COMPLEX LAMBDA(4),Q(4,4),QINV(4,4),D(4,4),ANSWER(4,4)
C
C***********************************************************************
C
        DO 50 I=1,4
         DO 50 J=1,4
            D(J,I) = CMPLX(0.0,0.0)
50      CONTINUE
C
        DO 60 I=1,4
            D(I,I) = CEXP(LAMBDA(I)*Z)
60      CONTINUE
C
        CALL CMATMULT3(Q,D,QINV,ANSWER)
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the exponential of a 4x4 matrix.       ***
C*** i.e. ANSWER = exp(kappa*z) where kappa is th 4x4 matrix.        ***
c*** The matrix KAPPA is supplied in terms of its eigenvalue matrices***
c*** such that ANSWER = Q*D*QINV  is a 4x4 matrix.                   ***
C*** D, Q, and QINV are 4x4 complex matrices.                        ***
C*** ANSWER is a 4x4 real materix.                                   ***
c***    Calling routine:      TRUNK                                  ***
c***    Called subroutines:   none                                   ***
c*********************************************************************** 
C
        SUBROUTINE EXP_MAT(LAMBDA,Q,QINV,z,ANSWER)
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
C
C   LAMBDA(4)   =   4 ELEMENT VECTOR CONTAINING THE EIGENVALUES(COMPLEX)
C   Q(4,4)      =   4x4 COMPLEX EIGENMATRIX
C   QINV(4,4)   =   4x4 COMPLEX INVERVSE OF Q(4,4)
C   Z           =   EXTENT OF PROPAGATION THE VOLUME(METERS)
C               =   (H/COS(THETA))
C   D(4,4)      =   COMPLEX DIAGONAL MATRIX
C   ANSWER(4,4) =   REAL OUTPUT MATRIX
C
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
C
        INTEGER I,J,K
        REAL Z ,ANSWER(4,4)
        REAL Qr(4,4),QINVr(4,4),Lr(4,4),Qi(4,4),QINVi(4,4),Li(4,4)
        REAL FACTRR(4,4), FACTRI(4,4), FACTII(4,4), FACTIR(4,4)
        REAL DIFF(4,4), SUM(4,4)
        COMPLEX LAMBDA(4),Q(4,4),QINV(4,4),D(4,4)
C
C***********************************************************************
C
        DO 50 I=1,4
         DO 45 J=1,4
            D(J,I) = CMPLX(0.0,0.0)
            Lr(J,I) = 0.0
            Li(J,I) = 0.0
            Qr(J,I) = REAL(Q(J,I))
            Qi(J,I) = AIMAG(Q(J,I))
            QINVr(J,I) = REAL(QINV(J,I))
            QINVi(J,I) = AIMAG(QINV(J,I))
            ANSWER(J,I) = 0.0
45        CONTINUE
50      CONTINUE
C
        DO 60 I=1,4
            D(I,I) = CEXP(LAMBDA(I)*Z)
            Lr(I,I) = REAL(D(I,I))
            Li(I,I) = AIMAG(D(I,I))
60      CONTINUE
C
        DO 100 I=1,4
         DO 100 J=1,4
            FACTRR(J,I) = QINVR(J,I)*LR(J,J) 
            FACTRI(J,I) = QINVR(J,I)*LI(J,J) 
            FACTII(J,I) = QINVI(J,I)*LI(J,J) 
            FACTIR(J,I) = QINVI(J,I)*LR(J,J) 
100     CONTINUE
C
        DO 120 I=1,4
         DO 120 J=1,4
            DIFF(J,I) = FACTRR(J,I) - FACTII(J,I)
            SUM(J,I)  = FACTRI(J,I) + FACTIR(J,I)
120     CONTINUE
C
        DO 200 K=1,4
        DO 200 J=1,4
        DO 200 I=1,4
         ANSWER(K,J) = ANSWER(K,J) + QR(K,I)*DIFF(I,J)-QI(K,I)*SUM(I,J)
200     CONTINUE
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C*** This subroutine computes the inverse of a 4x4 complex matrix    ***
C***********************************************************************
C
        SUBROUTINE INVERSE_MAT(INPUT,OUTPUT)
        save
C
C***********************************************************************
C
        INTEGER IV(4),JOB,I,J
        REAL RC
        COMPLEX SV(4),DET(2)
        COMPLEX INPUT(4,4),OUTPUT(4,4)
C
C***********************************************************************
C

        DO 10 I=1,4
         DO 10 J=1,4
            OUTPUT(J,I)=INPUT(J,I)
10      CONTINUE 
C
C   USE NAAS ROUTINES TO INVERT THE MATRIX (LINPACK)
C     OUTPUT = MATRIX TO BE INVERTED
C     ADIM   = ROW DIMENSION OF THE SYSTEM (ADIM = 4)
C     N      = ORDER OF SYSTEM (N = 4) 
C     IV     = INTEGER VECTOR OF LENGTH 4
C     RC     = ESTIMATE OF 1/COND(OUTPUT)
C     SV     = SCRATH VECTOR OF LENGTH N=4
C     DET(2) = ONE DIMENSIONAL ARRAY REPRESENTING DETERMINANT OF OUTPUT
C     JOB    = INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
        CALL CGECO(OUTPUT,4,4,IV,RC,SV)   
        JOB = 1
        CALL CGEDI(OUTPUT,4,4,IV,DET,SV,JOB) 
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***     This subroutine multiplies two 4x4 complex matrices.       ****
c*********************************************************************** 
C
        SUBROUTINE CMATMULT2(MAT1,MAT2,OUTPUT)
        save
C
C***********************************************************************
C
        INTEGER I,J,K
        COMPLEX MAT1(4,4),MAT2(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        DO 10 I=1,4
         DO 10 J=1,4
            OUTPUT(J,I) = CMPLX(0.0,0.0)
10      CONTINUE
C
        DO 50 K=1,4
          DO 50 J=1,4
            DO 50 I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
50      CONTINUE
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***     This subroutine multiplies three 4x4 complex matrices.     ****
c***    Calling routine:      EXP_KAPPA                              ***
c***    Called subroutines:   CMATMULT2                              ***
c*********************************************************************** 
C
        SUBROUTINE CMATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
C
C***********************************************************************
C
        COMPLEX MAT1(4,4),MAT2(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        CALL CMATMULT2(MAT1,MAT2,MAT4)
        CALL CMATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***    This subroutine multiplies three 4x4  matrices.             ****
C***    They are multiplied as AxBxC  where A and C are complex and  ***
C***    B is real.                                                   ***
c***    Calling routine:      EXP_KAPPA                              ***
c***    Called subroutines:   CMATMULT2                              ***
c*********************************************************************** 
C
        SUBROUTINE CrMATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
C
C***********************************************************************
C
        REAL MAT2(4,4)
        COMPLEX MAT1(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        CALL CrMATMULT2(MAT1,MAT2,MAT4)
        CALL CMATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***     This subroutine multiplies two 4x4 matrices.               ****
C***     MAT1 is complex and MAT2 is real.                          ****
c*********************************************************************** 
C
        SUBROUTINE CrMATMULT2(MAT1,MAT2,OUTPUT)
        save
C
C***********************************************************************
C
        INTEGER I,J,K
        REAL MAT2(4,4)
        COMPLEX MAT1(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        DO 10 I=1,4
         DO 10 J=1,4
            OUTPUT(J,I) = CMPLX(0.0,0.0)
10      CONTINUE
C
        DO 50 K=1,4
          DO 50 J=1,4
            DO 50 I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
50      CONTINUE
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***    This subroutine multiplies three 4x4 real matrices.         ****
c***    Calling routine:      SOLVE_CANOPY                           ***
c***    Called subroutines:   MATMULT2                               ***
c*********************************************************************** 
C
        SUBROUTINE MATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
C
C***********************************************************************
C
        REAL MAT1(4,4),MAT2(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        CALL MATMULT2(MAT1,MAT2,MAT4)
        CALL MATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C***     This subroutine multiplies two 4x4 real matrices.          ****
c*********************************************************************** 
C
        SUBROUTINE MATMULT2(MAT1,MAT2,OUTPUT)
        save
C
C***********************************************************************
C
        INTEGER I,J,K
        REAL MAT1(4,4),MAT2(4,4),OUTPUT(4,4)
C
C***********************************************************************
C
        DO 10 I=1,4
         DO 10 J=1,4
            OUTPUT(J,I) = 0.0
10      CONTINUE
C
        DO 50 K=1,4
          DO 50 J=1,4
            DO 50 I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
50      CONTINUE
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the product of three 4x4 matrices.     ***
C*** The ANSWER = Q*D*QINV where ANSWER is real-valued.              ***
C*** D, Q, and QINV are 4x4 complex matrices.                        ***
C*** ANSWER is a 4x4 real materix.                                   ***
c***    Calling routine:      TRUNK                                  ***
c***    Called subroutines:   none                                   ***
c*********************************************************************** 
C
        SUBROUTINE MATMULT3Cr(Q,A,QINV,ANSWER)
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
C
C   Q(4,4)      =   4x4 COMPLEX EIGENMATRIX
C   QINV(4,4)   =   4x4 COMPLEX INVERVSE OF Q(4,4)
C   D(4,4)      =   4x4 COMPLEX MATRIX
C   ANSWER(4,4) =   REAL OUTPUT MATRIX
C
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
C
        INTEGER I,J,K
        REAL ANSWER(4,4)
        REAL Qr(4,4),QINVr(4,4),Lr(4,4),Qi(4,4),QINVi(4,4),Li(4,4)
        REAL FACTRR(4,4), FACTRI(4,4), FACTII(4,4), FACTIR(4,4)
        REAL DIFF(4,4), SUM(4,4)
        COMPLEX Q(4,4),QINV(4,4),A(4,4)
C
C***********************************************************************
C
        DO 50 I=1,4
         DO 50 J=1,4
            Lr(J,I) = REAL(A(J,I))
            Li(J,I) = AIMAG(A(J,I))
            Qr(J,I) = REAL(Q(J,I))
            Qi(J,I) = AIMAG(Q(J,I))
            QINVr(J,I) = REAL(QINV(J,I))
            QINVi(J,I) = AIMAG(QINV(J,I))
            ANSWER(J,I) = 0.0
            FACTRR(I,J) = 0.0
            FACTII(I,J) = 0.0
            FACTRI(I,J) = 0.0
            FACTIR(I,J) = 0.0
50      CONTINUE
C
C***********************************************************************
C
        DO 120 J = 1, 4
         DO 110 I = 1, 4
          DO 100 K = 1, 4
            FACTRR(I,J) = FACTRR(I,J) + QINVR(K,J)*LR(I,K)
            FACTII(I,J) = FACTII(I,J) + QINVI(K,J)*LI(I,K)
            FACTRI(I,J) = FACTRI(I,J) + QINVR(K,J)*LI(I,K)
            FACTIR(I,J) = FACTIR(I,J) + QINVI(K,J)*LR(I,K)
100       CONTINUE
110      CONTINUE
120     CONTINUE
C
        DO 150 J = 1, 4
         DO 140 I = 1, 4
          DIFF(I,J) = FACTRR(I,J) - FACTII(I,J)
          SUM(I,J)  = FACTRI(I,J) + FACTIR(I,J)
140      CONTINUE
150     CONTINUE
C
        DO 220 J = 1,4
         DO 210 I = 1, 4
          DO 200 K = 1, 4
           ANSWER(I,J) = ANSWER(I,J) +
     &                   QR(I,K)*DIFF(K,J) - QI(I,K)*SUM(K,J)
200       CONTINUE
210      CONTINUE
220     CONTINUE
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the 4x4 real extinction matrix KAPPA   ***
c*** of the complex 2x2 input M-matrix M                             ***
c***    Calling routines:      TRUNK                                 ***
c***    Called subroutines:    none                                  ***
c*********************************************************************** 
C
        SUBROUTINE EXTINCT_SUB(M,KAPPA)
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
c
c   M(2,2)      =   2x2 complex matrix (input)
c   KAPPA(4,4)  =   4x4 real extinction matrix (output)
c
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
c
        REAL KAPPA(4,4)
c
        COMPLEX M(2,2)
c
c***********************************************************************
c
        KAPPA(1,1) = -2.0*REAL(M(1,1))
        KAPPA(2,1) =  0.0
        KAPPA(3,1) = -2.0*REAL(M(2,1))
        KAPPA(4,1) =  2.0*AIMAG(M(2,1))
        KAPPA(1,2) =  0.0
        KAPPA(2,2) = -2.0*REAL(M(2,2))
        KAPPA(3,2) = -2.0*REAL(M(1,2))
        KAPPA(4,2) = -2.0*AIMAG(M(1,2))
        KAPPA(1,3) = -REAL(M(1,2))
        KAPPA(2,3) = -REAL(M(2,1))
        KAPPA(3,3) = -(REAL(M(1,1)) + REAL(M(2,2))) 
        KAPPA(4,3) = -(AIMAG(M(1,1)) - AIMAG(M(2,2))) 
        KAPPA(1,4) = -2.0*REAL(M(1,2))
        KAPPA(2,4) = -2.0*REAL(M(2,1))
        KAPPA(3,4) =  (AIMAG(M(1,1)) - AIMAG(M(2,2))) 
        KAPPA(4,4) = -(REAL(M(1,1)) + REAL(M(2,2))) 
C
        RETURN
        END
c
c***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the eigenvalue solution for the        ***
c*** propagation of the mean field for the complex 2x2 input         ***
c*** M-matrix M.                                                     ***
c***    Calling routines:      TRUNK                                 ***
c***                           CROWN_LAYER                           ***
c***    Called subroutines:    INVERSE_MAT                           ***
c***    Called functions:      CSQRT                                 ***
c*********************************************************************** 
C
        SUBROUTINE MAT_EIGEN_SUB(M,LAMBDA,Q,QINV)
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
c
c   M(2,2)      =   2x2 complex M matrix (input) 
C                   -- AVERAGE SCATTEING MATRIX TIMES 2*PI*N/k0
C   LAMBDA(4)   =   4 ELEMENT VECTOR CONTAINING THE EIGENVALUES(COMPLEX)
C   Q(4,4)      =   4x4 COMPLEX EIGENMATRIX
C   QINV(4,4)   =   4x4 COMPLEX INVERVSE OF Q(4,4)
C   K1          =   EIGENVALUE OF MEAN PROPAGATING v-FIELD
C   KV1         =   EIGENVECTOR CORRESPONDING TO K1
C   K2          =   EIGENVALUE OF MEAN PROPAGATING h-FIELD
C   KV2         =   EIGENVECTOR CORRESPONDING TO K2
C   r           =   FACTOR USED TO COMPUTE K1 AND K2
c   b1          =   SECOND ELEMENT OF KV1
c   b2          =   FIRST ELEMENT OF KV2
C   TOL         =   TOLLERANCE USED TO DETERMINE THE IMPORTANCE OF 
C                   CROSS-POLARIZATION
C   J           =  SQRT(-1)
C   k0          = WAVE NUMBER (1/METERS)
c
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
c
        REAL TOL, WORK1, WORK2, CHECK1, CHECK2
        REAL FREQ_HERTZ, WAVELENGTH, k0
c
        COMPLEX M(2,2),LAMBDA(4),Q(4,4),QINV(4,4),KV1(2),KV2(2),K1,K2
        COMPLEX r,b1,b2,CDIFF,CSUM,CWORK1,CWORK2
        COMPLEX J, MYCSQRT
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
c
c***********************************************************************
C
        J = CMPLX(0.0,1.0)
        TOL = 1.0E-3
C
C***********************************************************************
C
        IF(CABS(M(1,2)).EQ.0.0)THEN
            CHECK1 = 0.0
        ELSE IF(CABS(M(1,1)).NE.0.0)THEN
            CHECK1 = CABS(M(1,2)/M(1,1))
        ELSE
            CHECK1 = TOL + 1.0
        ENDIF
C
        IF(CABS(M(2,1)).EQ.0.0)THEN
            CHECK2 = 0.0
        ELSE IF(CABS(M(2,2)).NE.0.0)THEN
            CHECK2 = CABS(M(2,1)/M(2,2))
        ELSE
            CHECK2 = TOL + 1.0
        ENDIF
C
        IF((CHECK1.LT.TOL).AND.(CHECK2.LT.TOL)) THEN
C
C--- CASE FOR NO CROSS POLARIZATION ---
C
            K1 = k0 - J*M(1,1)
            K2 = k0 - J*M(2,2)
C
            KV1(1) = CMPLX(1.0,0.0)
            KV1(2) = CMPLX(0.0,0.0)
            KV2(1) = CMPLX(0.0,0.0)
            KV2(2) = CMPLX(1.0,0.0)
C
            WORK1 = -REAL(M(1,1) + M(2,2))
            WORK2 = AIMAG(M(1,1) - M(2,2))
C
            LAMBDA(1) = -2.0*REAL(M(1,1))
            LAMBDA(2) =  CMPLX(WORK1,-WORK2)
            LAMBDA(3) =  CMPLX(WORK1,WORK2)
            LAMBDA(4) = -2.0*REAL(M(2,2))
C
            Q(1,1) = CMPLX(1.0,0.0)
            Q(2,1) = CMPLX(0.0,0.0)
            Q(3,1) = CMPLX(0.0,0.0)
            Q(4,1) = CMPLX(0.0,0.0)
            Q(1,2) = CMPLX(0.0,0.0)
            Q(2,2) = CMPLX(0.0,0.0)
            Q(3,2) = CMPLX(1.0,0.0)
            Q(4,2) = CMPLX(0.0,-1.0)
            Q(1,3) = CMPLX(0.0,0.0)
            Q(2,3) = CMPLX(0.0,0.0)
            Q(3,3) = CMPLX(1.0,0.0)
            Q(4,3) = CMPLX(0.0,1.0)
            Q(1,4) = CMPLX(0.0,0.0)
            Q(2,4) = CMPLX(1.0,0.0)
            Q(3,4) = CMPLX(0.0,0.0)
            Q(4,4) = CMPLX(0.0,0.0)
c
            QINV(1,1) = CMPLX(1.0,0.0)
            QINV(2,1) = CMPLX(0.0,0.0)
            QINV(3,1) = CMPLX(0.0,0.0)
            QINV(4,1) = CMPLX(0.0,0.0)
            QINV(1,2) = CMPLX(0.0,0.0)
            QINV(2,2) = CMPLX(0.0,0.0)
            QINV(3,2) = CMPLX(0.0,0.0)
            QINV(4,2) = CMPLX(1.0,0.0)
            QINV(1,3) = CMPLX(0.0,0.0)
            QINV(2,3) = CMPLX(0.5,0.0)
            QINV(3,3) = CMPLX(0.5,0.0)
            QINV(4,3) = CMPLX(0.0,0.0)
            QINV(1,4) = CMPLX(0.0,0.0)
            QINV(2,4) = CMPLX(0.0,0.5)
            QINV(3,4) = CMPLX(0.0,-0.5)
            QINV(4,4) = CMPLX(0.0,0.0)
c
        else
C
C--- GENERAL CASE ---
C
            CDIFF = M(1,1)-M(2,2)
            CSUM =  M(1,1)+M(2,2)
C
C   note: must CSQRT be in the upper half-space?
C
            r = MYCSQRT(CDIFF*CDIFF +4.0*M(2,1)*M(1,2))
            K1 = k0 - (J/2.)*(CSUM + r)
            K2 = k0 - (J/2.)*(CSUM - r)
C
            b1 = 2.0*M(2,1)/(CDIFF+r)
            b2 = 2.0*M(1,2)/(-CDIFF-r)
C
            K1 = k0 - J*(M(1,1)+M(2,2)+r)/2.0
            K2 = k0 - J*(M(1,1)+M(2,2)-r)/2.0
C
            KV1(1) = CMPLX(1.0,0.0)
c            KV1(2) = CMPLX(0.0,b1)
            KV1(2) = b1
c            KV2(1) = CMPLX(b2,0.0)
            KV2(1) = b2
            KV2(2) = CMPLX(1.0,0.0)
C
            LAMBDA(1) = 2.0*AIMAG(K1)
            LAMBDA(2) = J*(CONJG(K2) - K1)
            LAMBDA(3) = J*(CONJG(K1) - K2)
            LAMBDA(4) = 2.0*AIMAG(K2)
C
            WORK1 = CABS(b1)
            WORK2 = CABS(b2)
            CWORK1 = CONJG(b1)
            CWORK2 = CONJG(b2)
C
            Q(1,1) =  CMPLX(1.0,0.0)
            Q(2,1) =  WORK1*WORK2
            Q(3,1) =  2.0*REAL(b1)
            Q(4,1) = -2.0*AIMAG(b1)
            Q(1,2) =  CWORK2
            Q(2,2) =  b1
            Q(3,2) =  1.0 + b1*CWORK2
            Q(4,2) = -J*(1.0 - b1*CWORK2)
            Q(1,3) =  b2
            Q(2,3) =  CWORK1
            Q(3,3) =  1.0 + b2*CWORK1
            Q(4,3) =  J*(1.0 - b2*CWORK1)
            Q(1,4) =  WORK2*WORK2
            Q(2,4) =  CMPLX(1.0,0.0)
            Q(3,4) =  2.0*REAL(b2)
            Q(4,4) =  2.0*AIMAG(b2)
C
            CALL INVERSE_MAT(Q,QINV)
        ENDIF
C
        RETURN
        END
c
c***********************************************************************
c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c   This subroutine computes the modified Stokes scattering operator
c   or the modified Mueller matrix based on
c        M_mod = [U]'[M][U] 
c   where means [U] transpose.
c***********************************************************************
c
        subroutine modify_mat(mat) 
        save
c
c***********************************************************************
c
        integer i, j
        real mat(4,4), u(4,4), utranspose(4,4), mat_mod(4,4)
c
c***********************************************************************
c
        u(1,1) =  1.0
        u(2,1) =  1.0
        u(3,1) =  0.0
        u(4,1) =  0.0
        u(1,2) =  1.0
        u(2,2) = -1.0
        u(3,2) =  0.0
        u(4,2) =  0.0
        u(1,3) =  0.0
        u(2,3) =  0.0
        u(3,3) =  1.0
        u(4,3) =  0.0
        u(1,4) =  0.0
        u(2,4) =  0.0
        u(3,4) =  0.0
        u(4,4) =  1.0
c
        do i=1,4
          do j=1,4
            utranspose(j,i) = u(i,j)
          enddo
        enddo
c
        call MATMULT3(utranspose,mat,u,mat_mod)
c
        do i=1,4
          do j=1,4
            mat(j,i) = mat_mod(j,i)
          enddo
        enddo
c
        return
        end
c
c***********************************************************************


