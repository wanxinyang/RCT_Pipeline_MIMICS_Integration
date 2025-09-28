C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C   Function to compute PDF(thetan) where thetan is the angle describing
C   the orientation of the needles. (assume azimuthal symmetry)
C----UNIFORM ORIENTATION----
C latest revision: 4-23-91 by Kyle.
C (actually implemented by LEP, 3-10-92)
C***********************************************************************
C
        REAL FUNCTION PDF_NDL(THETAN)
        save
C
C***********************************************************************
C
        REAL THETAN

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
        IF(I_PDF_NDL.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(thetan.GT.PI/2.)THEN
                PDF_NDL = 0.0
            ELSE IF(thetan.EQ.PI/2.)THEN
                PDF_NDL = 0.5
            ELSE
                PDF_NDL = SIN(thetan)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0.0).AND.(thetan.LT.(PI/2.0)))THEN
                  DUM = SIN(2*thetan)
                  DUM = DUM*DUM
                  PDF_NDL = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*SIN(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*SIN(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (4.0/PI)*DUM*DUM
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3/4)*SIN(THETA)**3  between 0 < theta < pi
C   or
C          (3/2)*SIN(THETA)**3  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.5)*DUM*DUM*DUM
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   PDF is (8/(3*pi))*SIN(THETA)**4  between 0 < theta < pi
C   or
C          (16/(3*pi))*SIN(THETA)**4  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.697653)*(DUM**4)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (0.848826)*(DUM**4)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   PDF is (15/16)*SIN(THETA)**5  between 0 < theta < pi
C   or
C          (15/8)*SIN(THETA)**5  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (15./8.)*(DUM**5)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (15./16.)*(DUM**5)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   PDF is (3.2/pi)*SIN(THETA)**6  between 0 < theta < pi
C   or
C          (6.4/pi)*SIN(THETA)**6  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (6.4/pi)*(DUM**6)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.8)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.093750)*SIN(THETA)**7  between 0 < theta < pi
C   or
C          (2.187499)*SIN(THETA)**7  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.187499)*(DUM**7)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.093750)*(DUM**7)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.328210)*(DUM**8)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.164105)*(DUM**8)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.460938)*(DUM**9)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.230469)*(DUM**9)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(thetan-OFFSET))**9 
     &               + ABS(SIN(PI-thetan-OFFSET))**9 
                 PDF_NDL = (1.23046)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60)**9)
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(thetan+OFFSET))**9 
     &               + ABS(SIN(PI-thetan+OFFSET))**9 
                 PDF_NDL = (1.23)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 PDF_NDL = (3.61311)*(sin(thetan))**20
            ELSE IF((thetan.EQ.PI/2.))THEN
                 PDF_NDL = (1.80656)*(sin(thetan))**20
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.14)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA+30 degrees)**9  between 0 < theta < pi
C   or  
C          (1.80656)*(SIN(THETA+30)**20 + SIN(180-THETA+30))**20
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(thetan+OFFSET))**20
     &               + ABS(SIN(PI-thetan+OFFSET))**20
                 PDF_NDL = (1.80656)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_NDL.EQ.17)THEN
C
C-----------------------------------------------------------------------
C   PDF = (1.0/0.616948)*(COS(THETA))**16  for 0 <  theta < PI
C       OR
C       = (2.0/0.616948)*(COS(THETA))**16 for 0 <  theta < PI/2
C                               in azimuthally symmetric case
C-----------------------------------------------------------------------
C
         IF(thetan.EQ.0.0)THEN
            dum = abs(cos(thetan))
            dum = dum**16.0
            PDF_NDL = 1.62088*dum
         ELSE IF((thetan.GT.0.0).AND.(thetan.LT.PI/2.))THEN
            dum = abs(cos(thetan))
            dum = dum**16.0
            PDF_NDL = 3.24176*dum
         ELSE
            PDF_NDL = 0.0
         ENDIF
C
C
C***********************************************************************
C

        ELSE
            print*,' *** Bad value for needle pdf function ***'
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
C   Function to compute PDF(thetad) where thetad is the angle describing
C   the orientation of the leaves. (assume azimuthal symmetry)
C----UNIFORM ORIENTATION----
C***********************************************************************
C
        REAL FUNCTION PDF_LF(THETAd)
        save
C
C***********************************************************************
C
        REAL THETAd, MU, NU ,TOL
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
        TOL = 8.727E-7
c
        IF(I_PDF_LEAF.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION
C   PDF is (0.5)*SIN(THETA)  between 0 < theta < pi
C   or
C          SIN(THETA)        between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF(THETAd.GT.PI/2.)THEN
                PDF_LF = 0.0
            ELSE IF(THETAd.EQ.PI/2.)THEN
                PDF_LF = 0.5
            ELSE
                PDF_LF = SIN(THETAd)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PLANOPHILE DISTRIBUTION   (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 2.770
                mu = 1.172
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   ERECTOPHILE DISTRIBUTION  (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.1720
                mu = 2.7700
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
C
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PLAGIOPHILE DISTRIBUTION  (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 3.326
                mu = 3.326
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
C
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.5)THEN
C
C-----------------------------------------------------------------------
C   EXTREMEOPHILE DISTRIBUTION   (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 0.433
                mu = 0.433
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
C
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.6)THEN
C
C-----------------------------------------------------------------------
C   UNIFORM DISTRIBUTION      (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.000
                mu = 1.000
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
C
C***********************************************************************
C
        ELSE IF(I_PDF_LEAF.EQ.7)THEN
C
C-----------------------------------------------------------------------
C   SPHEREICAL DISTRIBUTION   (BETA DISTRIBUTION)
C-----------------------------------------------------------------------
C
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.101
                mu = 1.930
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
C
C***********************************************************************
C
        ELSE
            CALL WRITE_ERROR(111)
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
C***********************************************************************
C
        REAL FUNCTION PDF_BR(THETAc)
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
        IF(I_PDF_BR_1.EQ.1)THEN
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
                PDF_BR = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR = 0.5
            ELSE
                PDF_BR = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.3)THEN
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
                 PDF_BR = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.4)THEN
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
                 PDF_BR = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.5)THEN
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
                 PDF_BR = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.6)THEN
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
                 PDF_BR = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.7)THEN
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
                 PDF_BR = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.8)THEN
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
                 PDF_BR = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.9)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.164105)*SIN(THETA)**8  between 0 < theta < pi
C   or
C          (2.328210)*SIN(THETA)**8  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.10)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.230469)*SIN(THETA)**9  between 0 < theta < pi
C   or
C          (2.460938)*SIN(THETA)**9  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.11)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23046)*SIN(THETA-30 degrees)**9  between 0 < theta < pi
C   or
C          (1.23046)*(SIN(THETA-30)**9 + SIN(180-THETA-30))
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR = (1.23046)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.12)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.23)*SIN(THETA+60 degrees)**9  between 0 < theta < pi
C   or  
C          (1.23)*(SIN(THETA+60)**9 + SIN(180-THETA+60))**9
C                                              between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR = (1.23)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.14)THEN
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
                 PDF_BR = (1.80656)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
C
        ELSE IF(I_PDF_BR_1.EQ.17)THEN
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
                 PDF_BR = (1.62088)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.18)THEN
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
                 PDF_BR = (1.1641)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_1.EQ.23)THEN
C
C-----------------------------------------------------------------------
C   PDF is (2/pi)*COS(THETA)**2  between 0 < theta < pi
C   or
C          (4/pi)*COS(THETA)**2  between 0 < theta < pi/2
C   for an azimthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = COS(THETAc)
                 PDF_BR = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = COS(THETAc)
                 PDF_BR = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
C
C***********************************************************************
C

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
C   (secondary branches)
C***********************************************************************
C
        REAL FUNCTION PDF_BR2(THETAc)
        save
C
C***********************************************************************
C
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF, I_PDF_BR_1, I_PDF_BR_2,
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
        IF(I_PDF_BR_2.EQ.1)THEN
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
                PDF_BR2 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR2 = 0.5
            ELSE
                PDF_BR2 = SIN(THETAc)
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (16/(3*pi))*SIN(2*THETA)**4  between 0 < theta < pi/2
C       =  0                                  OTHERWISE
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR2 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.3)THEN
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
                 PDF_BR2 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.4)THEN
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
                 PDF_BR2 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.5)THEN
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
                 PDF_BR2 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.6)THEN
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
                 PDF_BR2 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.7)THEN
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
                 PDF_BR2 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.8)THEN
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
                 PDF_BR2 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.9)THEN
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
                 PDF_BR2 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.10)THEN
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
                 PDF_BR2 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.11)THEN
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
                 PDF_BR2 = (1.23046)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.12)THEN
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
                 PDF_BR2 = (1.23)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.13)THEN
C
C-----------------------------------------------------------------------
C   PDF is (1.80656)*SIN(THETA)**20  between 0 < theta < pi
C   or  
C          (3.61311)*(SIN(THETA)**20  between 0 < theta < pi/2
C   for an azimuthally symmetric canopy.
C-----------------------------------------------------------------------
C
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR2 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR2 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.14)THEN
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
                 PDF_BR2 = (1.80656)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
C
        ELSE IF(I_PDF_BR_2.EQ.17)THEN
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
                 PDF_BR2 = (1.62088)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C
        ELSE IF(I_PDF_BR_2.EQ.18)THEN
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
                 PDF_BR2 = (1.1641)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
C
C***********************************************************************
C

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
C   the orientation of the trunks. (assume azimuthal symmetry)
C***********************************************************************
C
        REAL FUNCTION PDF_TR(THETAC)
        save
C
C***********************************************************************
C
        REAL THETAC, dum
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
        IF(I_PDF_TRUNK.EQ.1)THEN
C
C-----------------------------------------------------------------------
C   VERTICAL ORIENTATION
C   PDF = (1.0)/delta  for thetac = 0.0
C       = (0.0)        otherwise    
C-----------------------------------------------------------------------
C
            IF(THETAc.eq.0.0)THEN
                PDF_TR = 11.4592
            ELSE 
                PDF_TR = 0.0
            ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.2)THEN
C
C-----------------------------------------------------------------------
C   PDF = (1.0/0.85903)*(COS(THETAC))**8  for 0 <  thetac < PI
C       OR
C       = (2.0/0.85903)*(COS(THETAC))**8 for 0 <  thetac < PI/2
C                               in azimuthally symmetric case
C-----------------------------------------------------------------------
C
         IF(THETAC.EQ.0.0)THEN
            dum = abs(cos(thetac))
            dum = dum**8.0
            PDF_TR = 1.1641*dum
         ELSE IF((THETAC.GT.0.0).AND.(THETAC.LT.PI/2.))THEN
            dum = abs(cos(thetac))
            dum = dum**8.0
            PDF_TR = 2.32821*dum
         ELSE
            PDF_TR = 0.0
         ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.3)THEN
C
C-----------------------------------------------------------------------
C   PDF = (1.0/0.98175)*(COS(THETAC))**6  for 0 <  thetac < PI
C       OR
C       = (2.0/0.98175)*(COS(THETAC))**6 for 0 <  thetac < PI/2
C                               in azimuthally symmetric case
C-----------------------------------------------------------------------
C
         IF(THETAC.EQ.0.0)THEN
            dum = abs(cos(thetac))
            dum = dum**6.0
            PDF_TR = 1.01859*dum
         ELSE IF((THETAC.GE.0.0).AND.(THETAC.LT.PI/2.))THEN
            dum = abs(cos(thetac))
            dum = dum**6.0
            PDF_TR = 2.03718*dum
         ELSE
            PDF_TR = 0.0
         ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.4)THEN
C
C-----------------------------------------------------------------------
C   PDF = (8.0)*(EXP(-8.0*THETAC))  for 0 <  thetac < PI
C       ALSO
C       = (8.0)*(EXP(-8.0*THETAC))  for 0 <  thetac < PI/2
C                               in azimuthally symmetric case
C-----------------------------------------------------------------------
C
         IF((THETAC.GE.0.0).AND.(THETAC.LT.PI/2.))THEN
            PDF_TR = 8.0*(EXP(-8.0*THETAC))
         ELSE
            PDF_TR = 0.0
         ENDIF
C
C-----------------------------------------------------------------------
C
        ELSE
            PRINT*,'IMPROPER TRUNK LAYER PDF'
            STOP
        ENDIF
C
        RETURN
        END
C
C***********************************************************************
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C  SUBROUTINE for pdf of leaf angle distribution using a beta
c  distribution.
c***********************************************************************
c
        subroutine leaf_beta_dist(mu,nu,theta,pdf)
c
c***********************************************************************
c
c   The various distributions are from:
c   "Simple Beta Distribution Representation of Leaf Orientation
c        in Vegetation Canopies"
c     -- Narendra S. Goel and Donald E. Strebel
c     -- Agronomy Journalm Vol. 76 Sept.-Oct. 1984
c
c   In this routine, theta is the leaf angle from zenith (radians)
c
c   The values of i_choice correspond to:
c
c       i_choice = 1  -- planophile
c            nu = 2.770
c            mu = 1.172
c       i_choice = 2  -- erectophile
c            nu = 1.1720
c            mu = 2.7700
c       i_choice = 3  -- plagiophile
c            nu = 3.326
c            mu = 3.326
c       i_choice = 4  -- extremophile
c            nu = 0.433
c            mu = 0.433
c       i_choice = 5  -- uniform
c            nu = 1.000
c            mu = 1.000
c       i_choice = 6  -- spherical
c            nu = 1.101
c            mu = 1.930
c
c***********************************************************************
c
        real mu,nu,pi,pdf,theta,gammln
c
c***********************************************************************
c
        pi = 3.141592654
c
        pdf = (2./pi)*exp(gammln(mu+nu) - gammln(mu) - gammln(nu))
     &        *((1.0-theta/(pi/2.))**(nu-1.))*
     &          ((theta/(pi/2.))**(mu-1.))
c
        return
        end
c
c***********************************************************************
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
C**********************************************************************
      FUNCTION GAMMLN(XX)
      SAVE
C**********************************************************************
C
C  Function which returns the value of the Gamma function
C  (from Numerical Recipes by W. H. Press et al.)
C**********************************************************************
      INTEGER J
      REAL*4 GAMMLN,XX
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

c***********************************************************************
