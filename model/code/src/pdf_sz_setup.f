C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C   Subroutine to set up trunk size pdf integration parameters
C***********************************************************************
c
        subroutine pdf_tr_size_setup(s1,s2,PDFs1,PDFs2,ds1,ds2,Ns1,Ns2) 
c
        save
C***********************************************************************
C
C   PDFs1,PDFs2 = P.D.F. of size parameters s1 and s2
C   Ns1,Ns2  = number of integration steps for size parameters s1 and s2
C   ds1,ds2  = integration step sizes for size parameters s1 and s2 (cm)
C   s1,s2    = values of integration variables s1 and s2
C   s1 = cylinder length (cm)
C   s2 = cylinder diameter (cm)
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        INCLUDE 'parameters.include'
C----------------------------
c
        integer I
        integer Ns1,Ns2

        real s1(maxtrs1), s2(maxtrs2),ds1(maxtrs1), ds2(maxtrs2)
        real PDFs1(maxtrs1),PDFs2(maxtrs2),TRUNK_HGHT_FUNC

        INTEGER I_PDF_LEAF_SIZE, I_PDF_BR_1_SIZE, I_PDF_BR_2_SIZE,
     &          I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, I_PDF_BR_5_SIZE, 
     &          I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE, I_PDF_TRUNK_SIZE

        COMMON /I_PDF_SIZE/ I_PDF_LEAF_SIZE,I_PDF_BR_1_SIZE,
     &               I_PDF_BR_2_SIZE, I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, 
     &               I_PDF_BR_5_SIZE, I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE,
     &               I_PDF_TRUNK_SIZE


c
        REAL MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL TRUNK_Dsig, TRUNK_Hsig 

        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_TRUNK_SIG/ TRUNK_Dsig, TRUNK_Hsig 
C
        REAL TRUNK_DIAM_HIST_X(MAXTRDHIST),TRUNK_DIAM_HIST_Y(MAXTRDHIST)
        REAL TR_DIAM_HIST_DX(MAXTRDHIST)
        INTEGER I_TR_DIAM_HIST
C
        COMMON /TRUNK_HIST/ TRUNK_DIAM_HIST_X, TRUNK_DIAM_HIST_Y,
     &                      TR_DIAM_HIST_DX, I_TR_DIAM_HIST
C
C***********************************************************************
C
        IF(I_PDF_TRUNK_SIZE.EQ.0)THEN
            Ns1 = 1
            Ns2 = 1
            s1(1) = TRUNK_HGHT
            s2(1) = TRUNK_DIAM
            ds1(1) = 1.0
            ds2(1) = 1.0
            PDFs1(1) = 1.0
            PDFs2(1) = 1.0
        ELSE IF(I_PDF_TRUNK_SIZE.EQ.1)THEN
            Ns1 = 1
            Ns2 = I_TR_DIAM_HIST
            DO 10 I=1,Ns2
C
                s2(I) = TRUNK_DIAM_HIST_X(I)
                PDFs2(I) = TRUNK_DIAM_HIST_Y(I)
                ds2(I) = TR_DIAM_HIST_DX(I)
C
                PDFs1(I) = 1.0
                ds1(I) = 1.0
                s1(I) = TRUNK_HGHT_FUNC(s2(I))
c
10          CONTINUE
C
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.2)THEN
CC
C            Ns1 = INT(TRUNK_Hsig/TRUNK_HGHT)
C            Ns2 = INT(TRUNK_Dsig/TRUNK_DIAM)
C            do 30 i=1,Ns1
C                x = 
C                do 20 j=1,Ns2
C                    y = 
C                    f(i,j) = bivnorm(x,y,TRUNK_HGHT,TRUNK_DIAM,
C     &                           TRUNK_Hsig,TRUNK_Dsig,RHOtrHD)
C20              continue
C30          continue
C
C
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.3)THEN
C            GAUSS_DIST(X,TRUNK_HGHT,TRUNK_Hsig)
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.4)THEN
C            GAUSS_DIST(X,TRUNK_DIAM,TRUNK_Dsig)
        ELSE
            CALL WRITE_ERROR(115)
            STOP
        ENDIF
c
        return
        end
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C   function to compute bivariate normal distribution
C***********************************************************************
Cc
C        real function bivnorm(x,y,mux,muy,sigx,sigy,rho)
Cc
C        save
CC***********************************************************************
Cc
C        real x,y,mux,muy,sigx,sigy,rho
C        real dum1,dum2,dum3,dum4,dum5
Cc
CC----------------------------
Cc%include 'constants.include'
C        INCLUDE 'constants.include'
CC----------------------------
CC***********************************************************************
Cc
C        dum1 = 1.0 - rho*rho
C        dum2 = 1.0/(2.*pi*sigx*sigy*sqrt(dum1))
C        dum3 = -1.0/(2.0*dum1)
C        dum4 = (x-mux)/sigx
C        dum5 = (y-muy)/sigy
Cc
C        bivnorm = exp(dum3*(dum4*dum4-2.0*rho*dum4*dum5+dum5*dum5)
Cc
C        return
C        end
Cc
C***********************************************************************



