C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C   Function to compute PDF(thetac) where thetac is the angle describing
C   the orientation of the branches. (assume azimuthal symmetry)
C   (3rd branches)
C***********************************************************************
C
        REAL FUNCTION PDF_BR3(THETAc)
        save
C
C***********************************************************************
C
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C
        IF(I_PDF_BR_3.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR3 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR3 = 0.5
            ELSE
                PDF_BR3 = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR3 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*SIN(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*SIN(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3/4)*SIN(THETA)**3  between 0 < theta < pi
C   or
C          (3/2)*SIN(THETA)**3  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   PDF is (8/(3*pi))*SIN(THETA)**4  between 0 < theta < pi
C   or
C          (16/(3*pi))*SIN(THETA)**4  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   PDF is (15/16)*SIN(THETA)**5  between 0 < theta < pi
C   or
C          (15/8)*SIN(THETA)**5  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3.2/pi)*SIN(THETA)**6  between 0 < theta < pi
C   or
C          (6.4/pi)*SIN(THETA)**6  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.8)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.093750)*SIN(THETA)**7  between 0 < theta < pi
C   or
C          (2.187499)*SIN(THETA)**7  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR3 = (1.23046)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR3 = (1.23)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR3 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR3 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.14)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA+30 degrees)**9  between 0 < theta < pi
C   or  
C          (1.80656)*(SIN(THETA+30)**20 + SIN(180-THETA+30))**20
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR3 = (1.80656)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.17)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.62088)*SIN(THETA+64 degrees)**16  between 0 < theta < pi
C   or  
C          (1.62088)*(SIN(THETA+64)**16 + SIN(180-THETA+64))**16
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR3 = (1.62088)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_3.EQ.18)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.1641)*SIN(THETA+41 degrees)**8  between 0 < theta < pi
C   or  
C          (1.1641)*(SIN(THETA+41)**8 + SIN(180-THETA+41))**8
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR3 = (1.1641)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
C
C***********************************************************************

        ELSE
            CALL WRITE_ERROR(110)
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
C   Function to compute PDF(thetac) where thetac is the angle describing
C   the orientation of the branches. (assume azimuthal symmetry)
C   (4th branches)
C***********************************************************************
C
        REAL FUNCTION PDF_BR4(THETAc)
        save
C
C***********************************************************************
C
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C
        IF(I_PDF_BR_4.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR4 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR4 = 0.5
            ELSE
                PDF_BR4 = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR4 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*SIN(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*SIN(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3/4)*SIN(THETA)**3  between 0 < theta < pi
C   or
C          (3/2)*SIN(THETA)**3  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   PDF is (8/(3*pi))*SIN(THETA)**4  between 0 < theta < pi
C   or
C          (16/(3*pi))*SIN(THETA)**4  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   PDF is (15/16)*SIN(THETA)**5  between 0 < theta < pi
C   or
C          (15/8)*SIN(THETA)**5  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3.2/pi)*SIN(THETA)**6  between 0 < theta < pi
C   or
C          (6.4/pi)*SIN(THETA)**6  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.8)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.093750)*SIN(THETA)**7  between 0 < theta < pi
C   or
C          (2.187499)*SIN(THETA)**7  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR4 = (1.23046)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR4 = (1.23)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR4 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR4 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.14)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA+30 degrees)**9  between 0 < theta < pi
C   or  
C          (1.80656)*(SIN(THETA+30)**20 + SIN(180-THETA+30))**20
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR4 = (1.80656)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.17)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.62088)*SIN(THETA+64 degrees)**16  between 0 < theta < pi
C   or  
C          (1.62088)*(SIN(THETA+64)**16 + SIN(180-THETA+64))**16
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR4 = (1.62088)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_4.EQ.18)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.1641)*SIN(THETA+41 degrees)**8  between 0 < theta < pi
C   or  
C          (1.1641)*(SIN(THETA+41)**8 + SIN(180-THETA+41))**8
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR4 = (1.1641)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
C
C***********************************************************************

        ELSE
            CALL WRITE_ERROR(110)
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
C   Function to compute PDF(thetac) where thetac is the angle describing
C   the orientation of the branches. (assume azimuthal symmetry)
C   (5th branches)
C***********************************************************************
C
        REAL FUNCTION PDF_BR5(THETAc)
        save
C
C***********************************************************************
C
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C
        IF(I_PDF_BR_5.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR5 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR5 = 0.5
            ELSE
                PDF_BR5 = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR5 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*SIN(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*SIN(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3/4)*SIN(THETA)**3  between 0 < theta < pi
C   or
C          (3/2)*SIN(THETA)**3  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   PDF is (8/(3*pi))*SIN(THETA)**4  between 0 < theta < pi
C   or
C          (16/(3*pi))*SIN(THETA)**4  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   PDF is (15/16)*SIN(THETA)**5  between 0 < theta < pi
C   or
C          (15/8)*SIN(THETA)**5  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3.2/pi)*SIN(THETA)**6  between 0 < theta < pi
C   or
C          (6.4/pi)*SIN(THETA)**6  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.8)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.093750)*SIN(THETA)**7  between 0 < theta < pi
C   or
C          (2.187499)*SIN(THETA)**7  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR5 = (1.23046)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR5 = (1.23)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR5 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR5 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.14)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA+30 degrees)**9  between 0 < theta < pi
C   or  
C          (1.80656)*(SIN(THETA+30)**20 + SIN(180-THETA+30))**20
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR5 = (1.80656)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.17)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.62088)*SIN(THETA+64 degrees)**16  between 0 < theta < pi
C   or  
C          (1.62088)*(SIN(THETA+64)**16 + SIN(180-THETA+64))**16
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR5 = (1.62088)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_5.EQ.18)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.1641)*SIN(THETA+41 degrees)**8  between 0 < theta < pi
C   or  
C          (1.1641)*(SIN(THETA+41)**8 + SIN(180-THETA+41))**8
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR5 = (1.1641)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
C
C***********************************************************************
        ELSE
            CALL WRITE_ERROR(110)
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
C   Function to compute PDF(thetac) where thetac is the angle describing
C   the orientation of the branches. (assume azimuthal symmetry)
C   (6th branches)
C***********************************************************************
C
        REAL FUNCTION PDF_BR6(THETAc)
        save
C
C***********************************************************************
C
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C
        IF(I_PDF_BR_6.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR6 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR6 = 0.5
            ELSE
                PDF_BR6 = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR6 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*SIN(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*SIN(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3/4)*SIN(THETA)**3  between 0 < theta < pi
C   or
C          (3/2)*SIN(THETA)**3  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   PDF is (8/(3*pi))*SIN(THETA)**4  between 0 < theta < pi
C   or
C          (16/(3*pi))*SIN(THETA)**4  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   PDF is (15/16)*SIN(THETA)**5  between 0 < theta < pi
C   or
C          (15/8)*SIN(THETA)**5  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3.2/pi)*SIN(THETA)**6  between 0 < theta < pi
C   or
C          (6.4/pi)*SIN(THETA)**6  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.8)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.093750)*SIN(THETA)**7  between 0 < theta < pi
C   or
C          (2.187499)*SIN(THETA)**7  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR6 = (1.23046)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60))
C                                              between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR6 = (1.23)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR6 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR6 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.14)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA+30 degrees)**9  between 0 < theta < pi
C   or  
C          (1.80656)*(SIN(THETA+30)**20 + SIN(180-THETA+30))**20
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR6 = (1.80656)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.17)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.62088)*SIN(THETA+64 degrees)**16  between 0 < theta < pi
C   or  
C          (1.62088)*(SIN(THETA+64)**16 + SIN(180-THETA+64))**16
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR6 = (1.62088)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_6.EQ.18)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.1641)*SIN(THETA+41 degrees)**8  between 0 < theta < pi
C   or  
C          (1.1641)*(SIN(THETA+41)**8 + SIN(180-THETA+41))**8
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR6 = (1.1641)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
C
C***********************************************************************
        ELSE
            CALL WRITE_ERROR(110)
            STOP
        ENDIF
C  
        RETURN
        END
