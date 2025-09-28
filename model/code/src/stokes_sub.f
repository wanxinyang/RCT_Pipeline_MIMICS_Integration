C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C   Subroutine to compute the Stokes matrix of a target given its      c
C   scattering matrix.                                                 c
C***********************************************************************
C   Initial routine written 31 Aug 1988                                c
C***********************************************************************
C***    Calling routine:      TRUNK_PHASE                            ***
C***                          NEEDLE_PHASE                           ***
C***    Called subroutines:   none                                   ***
C***********************************************************************
C
        SUBROUTINE STOKES_SUB(SCAT,STOKES)
        save
C
C***********************************************************************
C                                                                      c
C   scat(2,2)   scattering matrix of the object (complex)              c
C   stokes(4,4) Stokes matrix (real) as in Eq.2-85b p.35 of Ishimaru   c
C   cw*****     Working variables (complex)                            c
C   work        Working variable (real)                                c
C                                                                      c
C***********************************************************************
C
        INTEGER I, J
        COMPLEX SCAT(2,2)
        COMPLEX CW1121C,CW1122C,CW1112C,CW2122C,CW1221C,CW1222C,CW1,CW2
        REAL STOKES(4,4),WORK
C                                                                      
C***********************************************************************
C
        DO 20 I=1,2
            DO 20 J=1,2
                WORK = CABS(SCAT(J,I))
                STOKES(J,I) = WORK*WORK
20      CONTINUE
C
        CW1121C = SCAT(1,1)*CONJG(SCAT(2,1))
        CW1222C = SCAT(1,2)*CONJG(SCAT(2,2))
        CW1112C = SCAT(1,1)*CONJG(SCAT(1,2))
        CW2122C = SCAT(2,1)*CONJG(SCAT(2,2))
C
        CW1122C = SCAT(1,1)*CONJG(SCAT(2,2))
        CW1221C = SCAT(1,2)*CONJG(SCAT(2,1))
C
        CW1 =  CW1122C + CW1221C
        CW2 =  CW1122C - CW1221C
C
        STOKES(3,1) = 2.0*REAL(CW1121C)
        STOKES(4,1) = 2.0*AIMAG(CW1121C)
C
        STOKES(3,2) = 2.0*REAL(CW1222C)
        STOKES(4,2) = 2.0*AIMAG(CW1222C)
C
        STOKES(1,3) = REAL(CW1112C)
        STOKES(2,3) = REAL(CW2122C)
        STOKES(3,3) = REAL(CW1)
        STOKES(4,3) = AIMAG(CW1)
C
        STOKES(1,4) = -AIMAG(CW1112C)
        STOKES(2,4) = -AIMAG(CW2122C)
        STOKES(3,4) = -AIMAG(CW2)
        STOKES(4,4) = REAL(CW2)
C
        RETURN
        END
C                                                                      
C***********************************************************************
