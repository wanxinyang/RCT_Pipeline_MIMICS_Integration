c***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
C*** This subroutine computes the reflectivity matrix of the         ***
c*** specular ground layer and the backscatter from a                ***
c*** rough ground layer.                                             ***
c***    Calling routine:      main MIMICS program                    ***
c***    Called subroutines:   GROUND_SPEC                            ***
C***                          GROUND_BACK                            ***
c*********************************************************************** 
C
        SUBROUTINE GROUND_LAYER
        save
C
c***********************************************************************
c------------------  variable definitions   ----------------------------
c***********************************************************************
C
C   THETA               = angle of incidence (radians)
C   EPSILONR(numeps)    = complex vector containing the relative 
C                         dielectrics of the canopy constituents
C   GRND_REFLECT(4,4)   = specular surface reflectivity matrix of ground
C   GRND_BACK_MAT(4,4)  = backscattering matrix for the ground surface.
C
c***********************************************************************
c------------------  variable declations    ----------------------------
c***********************************************************************
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        REAL GRND_REFLECT_MAT(4,4),GRND_BACK_MAT(4,4)
        REAL THETA, CTHETA, STHETA
        REAL FREQ_HERTZ, WAVELENGTH, k0
        real phi_hv_vv,phi_vh_vv,phi_hh_vv
        real  gnd_spec_p_diff(3), gnd_back_p_diff(3)
        real t_snow
        real p, q, alpha, beta, theta_loc, x, y
C
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL  CALL_SUB(N_CALLS,N_SUB_CALLS)
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
        COMMON /L_SURFACE_TYPE/ SOIL_SURFACE,WATER_SURFACE,
     &                          SNOW_SURFACE,ICE_SURFACE
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        common /r_snow/ t_snow
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
        complex eps_rel
        real two_way_tau
        integer i, j
C
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
        COMMON /R_GROUND_MATS/ GRND_REFLECT_MAT, GRND_BACK_MAT
        common /grnd_p_diff/ gnd_spec_p_diff, gnd_back_p_diff
C
C***********************************************************************
C   COMPUTE REFLECTIVITY MATRIX OF SPECULAR GROUND SURFACE
C***********************************************************************
C
        if((snow_surface).and.(t_snow.gt.0.0))then
            if(soil_surface)then
                eps_rel = EPSILONRC(2)/EPSILONRC(8)
            else if(water_surface)then
                eps_rel = EPSILONRC(1)/EPSILONRC(8)
            else
                print*,' ERROR in ground_layer --- bad surface type (a)'
            endif
c
            if(aimag(eps_rel).lt.0.) eps_rel = conjg(eps_rel)
c
c            theta_loc = arcsin(stheta/sqrt(real(epsilonrc(8))))
            alpha = abs(k0*aimag(csqrt(EPSILONRC(8))))
            beta = abs(k0*real(csqrt(EPSILONRC(8))))
            p = 2.0*alpha*beta
            q = beta*beta - alpha*alpha - k0*k0*stheta*stheta
            y = k0*stheta 
            x = (sqrt(sqrt(p*p + q*q) + q))/sqrt(2.0)
            theta_loc = atan2(y,x)
c
            call snow_layer(EPSILONRC(8),theta_loc,two_way_tau)
c
        else
c
            if(soil_surface)then
                eps_rel = EPSILONRC(2)
            else if(water_surface)then
                eps_rel = EPSILONRC(1)
            else
                print*,' ERROR in ground_layer --- bad surface type (b)'
            endif
            theta_loc = theta
        endif
c
        CALL GROUND_SPEC(eps_rel,theta_loc,GRND_REFLECT_MAT)
C
C***********************************************************************
C   COMPUTE BACKSCATTERING  MATRIX OF ROUGH GROUND SURFACE
C***********************************************************************
C
        CALL GROUND_BACK(THETA_loc,eps_rel,GRND_BACK_MAT)
c
        if((snow_surface).and.(t_snow.gt.0.0))then
          do i=1,4
            do j=1,4
               GRND_REFLECT_MAT(j,i) = two_way_tau*GRND_REFLECT_MAT(j,i)
               GRND_BACK_MAT(j,i) = two_way_tau*GRND_BACK_MAT(j,i)
            enddo
          enddo
        endif
c
        call phase_diff(GRND_REFLECT_MAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        gnd_spec_p_diff(1) = phi_hh_vv
        gnd_spec_p_diff(2) = phi_vh_vv
        gnd_spec_p_diff(3) = phi_hv_vv
c
        call phase_diff(GRND_BACK_MAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        gnd_back_p_diff(1) = phi_hh_vv
        gnd_back_p_diff(2) = phi_vh_vv
        gnd_back_p_diff(3) = phi_hv_vv
C
        RETURN
        END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  THIS SUBROUTINE COMPUTE THE REFLECTIVITY MATRIX OF A           ***
C***  SPECULAR GROUND SURFACE.                                       ***
C***    Calling routine:      GROUND_LAYER                           ***
C***    Called subroutines:   none                                   ***
C***********************************************************************
C                       
        SUBROUTINE GROUND_SPEC(eps_rel,theta_loc,R_MAT)
        save
C
C***********************************************************************
C-----------------Variable definitions----------------------------------
C***********************************************************************
c       
c   eps_rel     = dielectric of soil relative to snow or to free space
c               = epsrr + j*epsri
C   THETA_loc   = local incident angle - measured from nadir - radians
C   CTHETA_loc  = COS(THETA_loc)
C   STHETA_loc  = SIN(THETA_loc)
C
c   rperp       = fresnel reflection coefficient 
c                        - perpendicular polarization
c   rpar         = fresnel reflection coefficient 
c                        - parallel polarization
c   R_MAT(4,4)   = specular reflection matrix for the ground (real)
c
C***********************************************************************
C-----------------Variable declarations---------------------------------
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C           
        REAL R_MAT(4,4),WORK1,WORK2,WORK3,WORK4
        real theta_loc, stheta_loc, ctheta_loc
        COMPLEX RPERP,RPAR,CWORK,eps_rel
        REAL factor

        REAL FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0

        REAL RMS_SOIL, LS_SOIL
        COMMON /R_SURFACE/ RMS_SOIL, LS_SOIL
C
c
C***********************************************************************
c---------- Compute reflection coefficients ----------------------------
C***********************************************************************
C
        STHETA_loc = sin(theta_loc)
        CTHETA_loc = cos(theta_loc)
c
        CWORK = CSQRT(eps_rel - STHETA_loc*STHETA_loc)
C
        RPERP = (CTHETA_loc - CWORK)/(CTHETA_loc + CWORK)
        RPAR = (eps_rel*CTHETA_loc - CWORK)
     &         /(eps_rel*CTHETA_loc + CWORK)
C
C***********************************************************************
c----------compute reflectivity matrix----------------------------------
C***********************************************************************
C
        WORK1 = CABS(RPAR)
        WORK2 = CABS(RPERP)
        CWORK = RPAR*CONJG(RPERP)
        WORK3 = REAL(CWORK)
        WORK4 = AIMAG(CWORK)
C
        R_MAT(1,1) =  WORK1*WORK1
        R_MAT(2,1) =  0.0
        R_MAT(3,1) =  0.0
        R_MAT(4,1) =  0.0
        R_MAT(1,2) =  0.0
        R_MAT(2,2) =  WORK2*WORK2
        R_MAT(3,2) =  0.0
        R_MAT(4,2) =  0.0
        R_MAT(1,3) =  0.0
        R_MAT(2,3) =  0.0
        R_MAT(3,3) = -WORK3
        R_MAT(4,3) = -WORK4
        R_MAT(1,4) =  0.0
        R_MAT(2,4) =  0.0
        R_MAT(3,4) =  WORK4
        R_MAT(4,4) = -WORK3
c
C
C***********************************************************************
c----------modify reflectivity matrix to account for rough ground-------
C modified by Dr. Leland Pierce, 1-14-93.
C***********************************************************************
C
C      The equation for the modification is:
C
C      Gamma(new-specular) = Gamma(old_specular)*exp(-(2*Ks*cos(theta))**2)
C
C      where Ks = 2*pi*S/lambda,
C            where S is RMS height in meters,
C                  lambda is free-space wavelength in meters;
C      or Ks =  k0*S,
C            where k0 = 2*pi/lambda.
C
C      The variable RMS_SOIL is in cm, so we "*0.01" to convert to meters.
C
       factor = EXP( -(2.*k0*0.01*RMS_SOIL*CTHETA_loc)**2 )
       R_MAT(1,1)=R_MAT(1,1)*factor
       R_MAT(2,2)=R_MAT(2,2)*factor
       R_MAT(3,3)=R_MAT(3,3)*factor
       R_MAT(4,3)=R_MAT(4,3)*factor
       R_MAT(3,4)=R_MAT(3,4)*factor
       R_MAT(4,4)=R_MAT(4,4)*factor

       RETURN
       END      
c
C***********************************************************************
C***********************************************************************
C***  THIS SUBROUTINE COMPUTE THE BACKSCATTER MATRIX OF A            ***
C***  ROUGH GROUND SURFACE.                                          ***
C***    Calling routine:      GROUND_LAYER                           ***
C***    Called subroutines:   RUFGND                                 ***
C***********************************************************************
C
        SUBROUTINE GROUND_BACK(THETA,eps_rel,GRND_BACK_MAT)
        save
C
C***********************************************************************
C-----------------Variable definitions----------------------------------
C***********************************************************************
C
C   I = Model type for calculations
C         I = 0   Specular
C         I = 1   Geometrical Optic
C         I = 2   Physical Optics
C         I = 3   Small Perturbation
C         I = 4   UMich Empirical
C
C    eps_rel - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space.
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    M       - Rms surface slope
C    S       - Standard deviation of surface heights
C    LS      - Large scale surface correlation length
C    BACK    - .TRUE. if backscatter is desired and .FALSE. if not
C    WFTN    - Function representing the Fourier transform of the
C                correlation function of the surface.
C                (function of K1X,K1Y,K1XI,K1YI,LS)
C
C  OUTPUT:
C
C    MAT - 4X4 real matrix relating the incident and scattered
C            intensity vectors. (without range dependence)
C-----------------------------------------------------------------------
C     RMS_SOIL     =   surface RMS height  (centimeters)
C     LS_SOIL      =   surface correlation length (centimeters)
C    
C***********************************************************************
C-----------------Variable declarations---------------------------------
C***********************************************************************
C----------------------------
c%include 'parameters.include'
        include 'parameters.include'
C----------------------------
C
        INTEGER I, i_surf, j, jj
C
        REAL THETA, GRND_BACK_MAT(4,4), FREQ_HERTZ, WAVELENGTH, k0
        REAL RMS_SOIL, LS_SOIL
        REAL S, LS, M, THI, THS, PHS, PHI
        real mueller(4,4), prob, wftn
C
        COMPLEX eps_rel, WW(4,4)
C
        LOGICAL BACK
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
        EXTERNAL WFTN,PROB
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /R_SURFACE/ RMS_SOIL, LS_SOIL
        COMMON /L_SURFACE_TYPE/ SOIL_SURFACE,WATER_SURFACE,
     &                          SNOW_SURFACE,ICE_SURFACE
        common /surface_mod/ i_surf
C----------------------------
c%include 'constants.include'
        include 'constants.include'
C----------------------------
C***********************************************************************
C
C   Set up for the routine RUFGND.
C   This subroutine combines the five polarimetric rough surface
C   scattering models (Specular, Geometrical Optics, Physical Optics, 
C   Small Perturbation, and UMich Empirical).  The model type 
C   (one of the above) is selected with the call parameter I.
C
C         I = 0   Specular
C         I = 1   Geometrical Optics
C         I = 2   Physical Optics
C         I = 3   Small Perturbation
C         I = 4   UMich Empirical
C
C   Not all of the parameters in the call statement are used with
C   a given model type.  For example, only the geometrical optics
C   model uses the parameters BACK and M, and only the small perturbation
C   model uses the function W.  The call statements for each 
C   model are listed below to indicate the parameters used with each
C   model.
C
C***********************************************************************
C
C     For specular sufaces only:
      IF(i_surf.eq.0)THEN
        do jj=1,4
            do j=1,4
               GRND_BACK_MAT(j,jj) = 0.0
            enddo
        enddo 
        RETURN
      ENDIF
C
      IF(SOIL_SURFACE)THEN
        S = RMS_SOIL/100.
        LS = LS_SOIL/100.
C
        I = i_surf
        M = SQRT(2.0)*S/LS                                               {Use this value of M (actually you only need S and LS)
        BACK = .TRUE.                   
        THI = PI-THETA                                                   {Notice that this is PI-THETA for backscatter
        THS = THETA                                                      {Notice that this is TH0 for backscatter
        PHS = PI                                                         {Backscatter
        PHI = 0.0                                                        {Backscatter
C
c------ old call -------------------
c        CALL RUFGND(I,eps_rel,WAVELENGTH,THI,PHI,THS,PHS,M,S,LS,WFTN,
c     &                PROB,BACK,mueller,GRND_BACK_MAT)
c-------------------------------

        call RUFGND(I,1,eps_rel,WAVELENGTH,THI,PHI,THS,PHS,
     &        M,S,LS,1.,WFTN,PROB,BACK,.FALSE.,WW,1,GRND_BACK_MAT)

        call flip_stokes(GRND_BACK_MAT,GRND_BACK_MAT)
C
        RETURN
      else if (WATER_SURFACE) then
        do jj=1,4
            do j=1,4
               GRND_BACK_MAT(j,jj) = 0.0
            enddo
        enddo 
        RETURN
      ELSE
        CALL WRITE_ERROR(50)
        STOP
      ENDIF
      END
C
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C***  This subroutine computes backscattering from a soil surface.  ****
C***  This routine was originally written by Mike Whitt and was     ****
C***  obtained on 18 October 1988.  The updated/corrected version   ****
C***  below was installed 2-12-91.                                  ****
C***********************************************************************
      SUBROUTINE RUFGND(I,J,EPSR,WVL,THI,PHI,THS,PHS,M,S,L,LS,W,PROB,
     &                  BACK,SPEC,WW,IMM,MM)
      SAVE
C*********************************************************************
C
C   This subroutine combines the five polarimetric rough surface
C   scattering models (Specular, Geometrical Optics, Physical Optics, 
C   Small Perturbation, and UMich Empirical).  The model type 
C   (one of the above) is selected with the call parameter I.
C
C         I = 0   Specular
C         I = 1   Geometrical Optics
C         I = 2   Physical Optics
C         I = 3   Small Perturbation
C         I = 4   UMich Empirical
C
C  Not all of the parameters in the call statement are used with
C  a given model type.  For example, only the geometrical optics
C  model uses the parameters BACK and M, and only the small perturbation
C  model uses the function W.  The call statements for each 
C  model are listed below to indicate the parameters used with each
C  model.
C
C  Currently, only a Gaussian distributed surface is considered,
C  but other surface distributions can be examined with slight
C  modifications.  Also, the FSA convention is assumed throughout
C  the entire program.
C
C  The REAL*4 functions W and PROB must be explicitly typed and
C  declared EXTERNAL in the main driver.
C
C  Geometrical optics -- This routine has been checked against
C      Figure 12.2 in Ulaby et al., Vol. II.
C  
C      CALL GND_GO(EPSR,WVL,THI,PHI,THS,PHS,M,PROB,BACK,WW,IMM,MM)
C
C  Physical optics -- This routine has been checked against
C      Figures 12.4 and 12.3(a)-(c) in Ulaby et al., Vol. II.
C      In Fig. 12.4, the plot is the result of setting 
C      RPER1=RPAR1=0.  When CABS(EPSR) is small, RPER1 and
C      RPAR1 should be set to zero in the subroutine COEF_PO.
C
C      CALL GND_PO(J,EPSR,WVL,THI,PHI,THS,PHS,S,L,LS,SPEC,WW,IMM,MM)
C
C  Small perturbation -- This routine has been checked against
C      Figures 12.6(a), 12.6(b), 12.7(a), and 12.7(b) in Ulaby
C      et al., Vol. II.
C
C      CALL GND_SP(EPSR,WVL,THI,PHI,THS,PHS,S,L,W,WW,IMM,MM)
C
C
C  UMich Empirical -- This routine has been checked against
C      figs 13 and 14 in: 
C      Yisok Oh, K. Sarabandi, Fawwaz Ulaby: An Empirical Model
C      and an Inversion Technique for Radar Scattering from Bare
C      Soil Surfaces, IEEE Transactions on Geoscience and
C      Remote Sensing, vol 30., No. 2, 1992. pages 370-381.
C      Works only in backscatter.
C
C      CALL GND_EM(EPSR,WVL,THI,S,WW,MM)
C  
C---------------------------------------------------------------------
C
C  A description of all call parameters is given below.
C
C  INPUTS:
C
C    I       - Model type (1=geom. optics, 2=phys. optics, 3=small pert.)
C    J       - Gives type of PO model.  1=gaussian correlation,
C                2=arbitrary correlation (must create subroutine CORR)
C                3=corr #3 in Polarimetry book
C    EPSR    - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space. (METERS)
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    M       - Rms surface slope
C    S       - Standard deviation of surface heights (RMS SURF HEIGHT, METERS)
C    L       - Large scale surface correlation length
C    LS      - Small scale surface correlation length
C    W       - Function representing the Fourier transform of the
C                correlation function of the surface.
C                (function of K1X,K1Y,K1XI,K1YI,L)
C    PROB    - Represents the probability density function for
C                the surface slopes.
C    BACK    - .TRUE. if backscatter is desired and .FALSE. if not
C    SPEC    - .TRUE. if specular scatter is desired and .FALSE. if not
C    IMM     - IMM=0 means do not compute MM
C              IMM=1 means compute MM
C
C  OUTPUT:
C
C    WW  - 4x4 complex W matrix 
C    MM  - 4X4 real modified Mueller matrix relating the incident and 
C            scattered intensity vectors. (without range dependence)
C    
C---------------------------------------------------------------------
C
C  The following is a list of all the subroutines used by this
C  routine.
C
C  GEOMETRICAL OPTICS SUBROUTINES:
C
C  GND_GO
C  REFCOMP
C  DOTCOMP
C  PROB
C  
C  PHYSICAL OPTICS SUBROUTINES:  
C
C  GND_PO
C  COEF_PO
C  INIS_GAUSS
C  
C  SMALL PERTURBATION SUBROUTINES:
C
C  GND_SP
C  COEF_SP
C  W
C  
C  UMich Empirical SUBROUTINES:
C
C  GND_EM
C  
C  GENERAL SUBROUTINES:
C
C  W_TO_MM     - link to polar_subs.o and math_subs.o
C
C*********************************************************************
C
      INTEGER I,IMM,J
      REAL*4 WVL,THI,PHI,THS,PHS,M,S,L,LS,W,MM(4,4),PROB
      COMPLEX EPSR,WW(4,4)
      LOGICAL BACK,SPEC
      EXTERNAL W,PROB
C
      IF (BACK.AND.SPEC) THEN
        WRITE(*,*) 'Can''t have both backscatter and specular scatter'
        STOP
      ENDIF
C
      IF (I.EQ.1) THEN
        CALL GND_GO(EPSR,WVL,THI,PHI,THS,PHS,M,PROB,BACK,SPEC,
     &                  WW,IMM,MM)

      ELSE IF (I.EQ.2) THEN
        CALL GND_PO(J,EPSR,WVL,THI,PHI,THS,PHS,S,L,LS,SPEC,WW,IMM,MM)
      ELSE IF (I.EQ.3) THEN
        CALL GND_SP(EPSR,WVL,THI,PHI,THS,PHS,S,L,W,WW,IMM,MM)
      ELSE IF (I.EQ.4) THEN
        CALL GND_EM(EPSR,WVL,THI,S,WW,MM)
      ELSE
        WRITE(*,*) 'IMPROPER MODEL TYPE'
      ENDIF
C
      RETURN
      END


C*********************************************************************
C
C  UM Empirical SUBROUTINE
C
C  Written by Dr. Leland Pierce, 10-92
C
C*********************************************************************

      SUBROUTINE GND_EM(EPSR,WVL,THI,S,WW,MM)
      REAL*4 WVL,THI,S,MM(4,4)
      COMPLEX EPSR,WW(4,4)


C--------------------------------------------------------------------
C  UMich Empirical -- This routine has been checked against
C      figs 13 and 14 in: 
C      Yisok Oh, K. Sarabandi, Fawwaz Ulaby: An Empirical Model
C      and an Inversion Technique for Radar Scattering from Bare
C      Soil Surfaces, IEEE Transactions on Geoscience and
C      Remote Sensing, vol 30., No. 2, 1992. pages 370-381.
C      Works only in backscatter.
C
C--------------------------------------------------------------------
C  A description of all call parameters is given below.
C
C  INPUTS:
C
C    EPSR    - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space. -- IN METERS (LEP)
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    S       - Standard deviation of surface heights 
C                           -- RMS SURFACE HEIGHT IN METERS! (LEP)
C
C---------------------------------------------------------------------
C Local Variables:
C
C    lamcm   - wavelength in cm   (R)
C    scm     - RMS Surf Hgt in cm (R)
C    ks      - 2pi/lamcm * scm    (R)
C    n       - sqrt(epsr)         (CMPLX)
C    gam0    - gamma0 from paper  (R)
C    gamv    - gamma_v from paper (R)
C    gamh    - gamma_h from paper (R)
C    g,p,q   - from paper         (R)
C    rootp   - from paper         (R)
C    hhmag   - |S_hh|             (R)
C    vvmag   - |S_vv|             (R)
C    hvmag   - |S_hv|             (R)
C    ct      - cos(thi)           (R)
C    ct2     - cos(theta2)        (CMPLX)
C    sig0vv  - sigma^0 (vv)       (R)
C    incang  - incidence angle, radians (R)
C---------------------------------------------------------------------
C
      REAL lamcm, scm, ks, gam0, gamv, gamh
      REAL g, p, rootp, q
      REAL hhmag, vvmag, hvmag
      COMPLEX n, ctmp
      REAL ct
      COMPLEX ct2
      REAL sig0vv
      REAL incang

      include 'constants.include'

C---------------------------------------------------------------------
C      print *,'in ground_em:'
      incang   = PI - THI
C      print *,'EPSR,WVL,incang,S:',EPSR,WVL,incang/PI*180.,S

      lamcm = 100.*WVL
      scm   = 100.*S
      ks    = (2.*PI/lamcm)*scm
      n     = CSQRT(EPSR)
      ctmp  = (n-1.)/(n+1.)
      gam0  = CABS(ctmp)**2

C      print *,'lamcm, scm, ks, n, gam0:',lamcm, scm, ks, n, gam0

      ct    = cos(incang)
      ct2   = CSQRT( 1. - SIN(incang)**2/EPSR )

      ctmp  = (ct - n*ct2)/(ct + n*ct2)
      gamh  = CABS(ctmp)**2

      ctmp  = (n*ct - ct2)/(n*ct + ct2)
      gamv  = CABS(ctmp)**2

C      print *,'ct, ct2, gamh, gamv:',ct, ct2, gamh, gamv

      g     = 0.7*( 1. - exp( -0.65*(ks**1.8) ) )

      rootp = 1. - ((2.*incang/PI)**(1./(3.*gam0)))*exp(-ks)
      p     = rootp*rootp

      q     = 0.23*SQRT(gam0)*(1.-exp(-ks))

C      print *,'g,p,q:',g,p,q

      sig0vv=g*(ct**3)*(gamv+gamh)/rootp

C      print *,'sig0vv: lin, dB:',sig0vv,10.*ALOG10(sig0vv)
      
      vvmag = SQRT(  sig0vv/(4.*PI))
      hhmag = SQRT(p*sig0vv/(4.*PI))
      hvmag = SQRT(q*sig0vv/(4.*PI))

C      print *,'hh,vv,hv linear magnitudes:', hhmag, vvmag, hvmag
C      print *,'hh,vv,hv powe rin dB', 20.*ALOG10(hhmag), 
C     &             20.*ALOG10(vvmag), 20.*ALOG10(hvmag)

C
C     Compute W matrix (eqn 4.139 of book)
C
      WW(1,1) = CMPLX(vvmag**2,0.)
      WW(1,2) = CMPLX(hvmag**2,0.)
      WW(1,3) = CMPLX(0.,0.)
      WW(1,4) = CMPLX(0.,0.)

      WW(2,1) = CMPLX(hvmag**2,0.)
      WW(2,2) = CMPLX(hhmag**2,0.)
      WW(2,3) = CMPLX(0.,0.)
      WW(2,4) = CMPLX(0.,0.)

      WW(3,1) = CMPLX(0.,0.)
      WW(3,2) = CMPLX(0.,0.)
      WW(3,3) = CMPLX(hhmag*vvmag,0.)
      WW(3,4) = CMPLX(0.,0.)

      WW(4,1) = CMPLX(0.,0.)
      WW(4,2) = CMPLX(0.,0.)
      WW(4,3) = CMPLX(0.,0.)
      WW(4,4) = CMPLX(hhmag*vvmag,0.)


C
C  Compute the average modified Mueller matrix relating the incident and
C  scattered intensities for the rough surface.
C
      CALL W_TO_MM(WW,MM)
C
      RETURN
      END
C*********************************************************************
C
C  GEOMETRICAL OPTICS SUBROUTINES
C
C*********************************************************************
      SUBROUTINE GND_GO(EPSR,WVL,THI,PHI,THS,PHS,M,PROB,BACK,SPEC,
     &                  WW,IMM,MM)

      SAVE
C*********************************************************************
C
C  Program to compute the transformation matrix relating the
C  incident and scattered intensity vectors for a rough surface
C  under the Kirchhoff and stationary phase approximations 
C  (sometimes called the geometrical optics model).  
C  The incident medium is free space.
C
C  The routine is currently set up to use a surface height which 
C  is Gaussian distributed with zero mean.  However, it can be
C  modified (by changing the function PROB(HX,HY)) to use any 
C  other type of distribution.
C
C  INPUTS
C    
C    EPSR    - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space.
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    M       - Rms surface slope
C    PROB    - Represents the probability density function for
C                the surface slopes.
C    BACK    - .TRUE. if backscatter is desired and .FALSE. if not
C    SPEC    - .TRUE. if specular scatter is desired and .FALSE. if not
C    IMM     - IMM=0 means do not compute MM
C              IMM=1 means compute MM
C
C  OUTPUT
C
C    WW  - 4x4 complex W matrix 
C    MM  - 4X4 real modified Mueller matrix relating the incident and 
C            scattered intensity vectors. (without range dependence)
C
C*********************************************************************
C
      INTEGER IMM
      REAL*4 WVL,THI,PHI,THS,PHS,M,MM(4,4),PROB
      COMPLEX EPSR,WW(4,4)
      LOGICAL BACK,SPEC
C
      REAL*4 K1,PI,QX,QY,QZ,QT,Q,U0,D2,HXS,HYS
      REAL*4 STHI,CTHI,STHS,CTHS,SPHI,SPHS,CPHI,CPHS
      REAL*4 HINS,HSNI,VINS,VSNI,M1P
      COMPLEX UPQ(2,2),K2,RPER,RPAR
      EXTERNAL PROB
      DATA PI/3.141592654/
C
C  Initialize variables
C
      K1=2.0*PI/WVL
      K2=K1*CSQRT(EPSR)
      STHI=SIN(THI)
      STHS=SIN(THS)
      CTHI=COS(THI)
      CTHS=COS(THS)
      SPHI=SIN(PHI)
      SPHS=SIN(PHS)
      CPHI=COS(PHI)
      CPHS=COS(PHS)
      QX=K1*(STHI*CPHI-STHS*CPHS)
      QY=K1*(STHI*SPHI-STHS*SPHS)
      QZ=K1*(CTHI-CTHS)
      QT=SQRT(QX*QX+QY*QY)
      Q=SQRT(QX*QX+QY*QY+QZ*QZ)
C
C  Compute dot products and reflection coefficients
C
      CALL DOTCOMP(K1,K2,THI,PHI,THS,PHS,VINS,VSNI,HINS,HSNI)
      CALL REFCOMP(K1,K2,Q,QZ,BACK,RPER,RPAR)
C
C  Compute the terms UVV, UVH, UHV, UHH
C  These correspond to UPQ(1,1), UPQ(1,2), UPQ(2,1), UPQ(2,2)
C
      D2=HSNI**2+VSNI**2
      M1P=Q*ABS(QZ)/(K1*QZ*D2)
      UPQ(1,1)=-M1P*(RPER*HINS*HSNI+RPAR*VINS*VSNI)
      UPQ(1,2)= M1P*(RPER*VINS*HSNI-RPAR*HINS*VSNI)
      UPQ(2,1)= M1P*(RPER*HINS*VSNI-RPAR*VINS*HSNI)
      UPQ(2,2)=-M1P*(RPER*VINS*VSNI+RPAR*HINS*HSNI)
C
C  Compute the W matrix
C
      HXS=-QX/QZ
      HYS=-QY/QZ
      U0=Q*Q*K1*K1*PROB(HXS,HYS,M)/(4.0*QZ**4)
      WW(1,1)=U0*UPQ(1,1)*CONJG(UPQ(1,1))
      WW(1,2)=U0*UPQ(1,2)*CONJG(UPQ(1,2))
      WW(1,3)=U0*UPQ(1,2)*CONJG(UPQ(1,1))
      WW(1,4)=U0*UPQ(1,1)*CONJG(UPQ(1,2))
      WW(2,1)=U0*UPQ(2,1)*CONJG(UPQ(2,1))
      WW(2,2)=U0*UPQ(2,2)*CONJG(UPQ(2,2))
      WW(2,3)=U0*UPQ(2,2)*CONJG(UPQ(2,1))
      WW(2,4)=U0*UPQ(2,1)*CONJG(UPQ(2,2))
      WW(3,1)=U0*UPQ(2,1)*CONJG(UPQ(1,1))
      WW(3,2)=U0*UPQ(2,2)*CONJG(UPQ(1,2))
      WW(3,3)=U0*UPQ(2,2)*CONJG(UPQ(1,1))
      WW(3,4)=U0*UPQ(2,1)*CONJG(UPQ(1,2))
      WW(4,1)=U0*UPQ(1,1)*CONJG(UPQ(2,1))
      WW(4,2)=U0*UPQ(1,2)*CONJG(UPQ(2,2))
      WW(4,3)=U0*UPQ(1,2)*CONJG(UPQ(2,1))
      WW(4,4)=U0*UPQ(1,1)*CONJG(UPQ(2,2))
C
C  Compute the average modified Mueller matrix relating the incident and
C  scattered intensities for the rough surface.
C
      IF (IMM.EQ.1) THEN
        CALL W_TO_MM(WW,MM)
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE REFCOMP(K1,K2,Q,QZ,BACK,RPER,RPAR)
      SAVE
C*********************************************************************
C
C  This subroutine computes the reflection coefficients for
C  perpendicular and parallel polarizations.  The parameters
C  QX, QY, and QZ are the x, y, and z components of the term
C  -k1*(nshat-nihat).
C
C  INPUTS:
C
C    K1    - wave number in free space
C    K2    - wave number in dielectric
C    Q     - SQRT(QX**2 + QY**2 + QZ**2)
C    QZ    - z component of the term -k1*(nshat-nihat)
C    BACK  - .TRUE. if backscatter is desired
C            .FALSE. if not
C
C  OUTPUTS:
C
C    RPER - Reflection coeff. for perpendicular polarization
C    RPAR - Reflection coeff. for parallel polarization
C
C********************************************************************
      REAL*4 K1,Q,QZ
      COMPLEX K2,RPER,RPAR
      LOGICAL BACK
C
      REAL*4 CTHIL,STHIL
      COMPLEX STHTL,CTHTL,DEN1,DEN2
C
      IF (BACK) THEN
        CTHIL=1.0
        STHIL=0.0
        STHTL=0.0
        CTHTL=1.0
      ELSE
        CTHIL=-Q*ABS(QZ)/(2.0*K1*QZ)
        STHIL=SQRT(1.0-CTHIL*CTHIL)
        STHTL=(K1/K2)*STHIL
        CTHTL=CSQRT(1.0-STHTL*STHTL)
      ENDIF
C
      DEN1=K1*CTHIL+K2*CTHTL
      DEN2=K2*CTHIL+K1*CTHTL
      RPER=(K1*CTHIL-K2*CTHTL)/DEN1
      RPAR=(K2*CTHIL-K1*CTHTL)/DEN2
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE DOTCOMP(K1,K2,THI,PHI,THS,PHS,VINS,VSNI,HINS,HSNI)
      SAVE
C*********************************************************************
C
C  This subroutine computes the dot products between polarization
C  unit vectors and the unit vectors in the incident and scattered
C  directions.
C
C  INPUTS:
C
C    K1  - wave number in free space
C    K2  - wave number in dielectric
C    THI - Incident elevation angle in radians (measured from z-axis 
C            to the unit vector in the direction of 
C            propagation of the incident wave)
C    PHI - Incident azimuth angle in radians (measured from the
C            x-axis to the projection on the xy-plane of the unit
C            vector in the direction of propagation of the incident
C            wave)
C    THS - Scattered elevation angle in radians (measured from z-axis 
C            to the unit vector in the direction of 
C            propagation of the incident wave)
C    PHS - Scattered azimuth angle in radians (measured from the
C            x-axis to the projection on the xy-plane of the unit
C            vector in the direction of propagation of the incident
C            wave)
C
C  OUTPUTS:
C
C    VINS - vi dot ns
C    VSNI - vs dot ni
C    HINS - hi dot ns
C    HSNI - hs dot ni
C
C********************************************************************
      REAL*4 K1,THI,PHI,THS,PHS,VINS,VSNI,HINS,HSNI
      COMPLEX K2
C
      REAL*4 CTHI,STHI,CTHS,STHS,CPHSI,SPHSI
C
      CTHI=COS(THI)
      STHI=SIN(THI)
      CTHS=COS(THS)
      STHS=SIN(THS)
      CPHSI=COS(PHS-PHI)
      SPHSI=SIN(PHS-PHI)
C
      VINS=CTHI*STHS*CPHSI-STHI*CTHS
      VSNI=CTHS*STHI*CPHSI-STHS*CTHI
      HINS=STHS*SPHSI
      HSNI=-STHI*SPHSI
C 
      RETURN
      END
CC*********************************************************************
CC*********************************************************************
C      FUNCTION PROB(HX,HY,M)
C      SAVE
CC*********************************************************************
CC
CC  This function returns the value of the surface slope probability
CC  density at the X and Y slopes HX and HY, respectively.
CC
CC  INPUTS:
CC
CC    HX - Slope in the X-direction
CC    HY - Slope in the Y-direction
CC    M  - RMS surface slope
CC
CC********************************************************************
C      REAL*4 HX,HY,PROB,PI,M
C      DATA PI/3.141592654/
CC
C      PROB=EXP(-(HX*HX+HY*HY)/(2.0*M*M))
C      PROB=PROB/(2.0*PI*M*M)
CC
C      RETURN
C      END
C*********************************************************************
C
C  PHYSICAL OPTICS SUBROUTINES
C
C*********************************************************************
      SUBROUTINE GND_PO(POTYPE,EPSR,WVL,THI,PHI,THS,PHS,S,L,LS,
     &                  SPEC,WW,IMM,MM)
      SAVE
C*********************************************************************
C
C  Program to compute the transformation matrix relating the
C  incident and scattered intensity vectors for a rough surface
C  under the Kirchhoff and scalar approximations (sometimes called
C  the physical optics model).  The incident medium is free
C  space.
C
C  INPUTS
C    
C    POTYPE  - 1=gaussian correlation 
C              2=arbitrary correlation (must create subroutine CORR)
C    EPSR    - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space.
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    S       - Standard deviation of surface heights
C    L       - Large scale surface correlation length
C    LS      - Small scale surface correlation length
C    SPEC    - .TRUE. if specular scatter is desired and .FALSE. if not
C    IMM     - IMM=0 means do not compute MM
C              IMM=1 means compute MM
C
C  OUTPUT
C
C    WW  - 4x4 complex W matrix 
C    MM  - 4X4 real modified Mueller matrix relating the incident and 
C            scattered intensity vectors. (without range dependence)
C
C*********************************************************************
C
      INTEGER IMM,POTYPE
      REAL*4 WVL,THI,PHI,THS,PHS,S,L,LS,MM(4,4)
      COMPLEX EPSR,WW(4,4)
C
      INTEGER I,J,M,N
      REAL*4 K1,PI,QX,QY,QZ,QT,C0,TMP
      REAL*4 STHI,CTHI,STHS,CTHS,SPHI,SPHS,CPHI,CPHS
      COMPLEX IN(2,2,2,2),IS(2,2,2,2),IC(2,2,2,2)
      COMPLEX K2,APQ(2,2),BPQ(2,2),CPQ(2,2)
      LOGICAL SPEC
      DATA PI/3.141592654/
C
C  Initialize variables
C
      K1=2.0*PI/WVL
      K2=K1*CSQRT(EPSR)
      STHI=SIN(THI)
      STHS=SIN(THS)
      CTHI=COS(THI)
      CTHS=COS(THS)
      SPHI=SIN(PHI)
      SPHS=SIN(PHS)
      CPHI=COS(PHI)
      CPHS=COS(PHS)
      QX=K1*(STHI*CPHI-STHS*CPHS)
      QY=K1*(STHI*SPHI-STHS*SPHS)
      QZ=K1*(CTHI-CTHS)
      QT=SQRT(QX*QX+QY*QY)
C
C  Compute APQ and BPQ coefficients
C
      CALL COEF_PO(K1,K2,THI,PHI,THS,PHS,APQ,BPQ,CPQ)
C
C  Check for specular scatter
C
      IF (SPEC) THEN
        TMP=EXP(-QZ*QZ*S*S)
        DO I=1,2
          DO J=1,2
            DO M=1,2
              DO N=1,2
                IC(I,J,M,N)=TMP*APQ(I,J)*CONJG(APQ(M,N))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I=1,2
          DO J=1,2
            DO M=1,2
              DO N=1,2
                IC(I,J,M,N)=CMPLX(0.0,0.0)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C  Compute terms IN (incoherent or Non-coherent) and IS (due to 
C  surface slope).  Because this part of the code involves numerical 
C  integration, it takes the most time.
C
      IF (POTYPE.EQ.1) THEN
        CALL INIS_GAUSS(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,IN,IS)
      ELSE IF (POTYPE.EQ.2) THEN
        CALL INIS_ARB(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,LS,IN,IS)
      ELSE IF (POTYPE.EQ.3) THEN
        CALL INIS_CORR3(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,IN,IS)
      ELSE
        PAUSE 'Improper type for PO model'
      ENDIF
C
C  Compute the W matrix for both incoherent and surface slope
C  contributions.
C
      C0=K1**2/(4.0*PI)**2
C
      WW(1,1)=C0*(IN(1,1,1,1)+IS(1,1,1,1)+IC(1,1,1,1))
      WW(1,2)=C0*(IN(1,2,1,2)+IS(1,2,1,2)+IC(1,2,1,2))
      WW(1,3)=C0*(IN(1,2,1,1)+IS(1,2,1,1)+IC(1,2,1,1))
      WW(1,4)=C0*(IN(1,1,1,2)+IS(1,1,1,2)+IC(1,1,1,2))
      WW(2,1)=C0*(IN(2,1,2,1)+IS(2,1,2,1)+IC(2,1,2,1))
      WW(2,2)=C0*(IN(2,2,2,2)+IS(2,2,2,2)+IC(2,2,2,2))
      WW(2,3)=C0*(IN(2,2,2,1)+IS(2,2,2,1)+IC(2,2,2,1))
      WW(2,4)=C0*(IN(2,1,2,2)+IS(2,1,2,2)+IC(2,1,2,2))
      WW(3,1)=C0*(IN(2,1,1,1)+IS(2,1,1,1)+IC(2,1,1,1))
      WW(3,2)=C0*(IN(2,2,1,2)+IS(2,2,1,2)+IC(2,2,1,2))
      WW(3,3)=C0*(IN(2,2,1,1)+IS(2,2,1,1)+IC(2,2,1,1))
      WW(3,4)=C0*(IN(2,1,1,2)+IS(2,1,1,2)+IC(2,1,1,2))
      WW(4,1)=C0*(IN(1,1,2,1)+IS(1,1,2,1)+IC(1,1,2,1))
      WW(4,2)=C0*(IN(1,2,2,2)+IS(1,2,2,2)+IC(1,2,2,2))
      WW(4,3)=C0*(IN(1,2,2,1)+IS(1,2,2,1)+IC(1,2,2,1))
      WW(4,4)=C0*(IN(1,1,2,2)+IS(1,1,2,2)+IC(1,1,2,2))
C
C  Compute the average modified Mueller matrix for the rough
C  surface.
C
      IF (IMM.EQ.1) THEN
        CALL W_TO_MM(WW,MM)
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE COEF_PO(K1,K2,THI,PHI,THS,PHS,APQ,BPQ,CPQ)
      SAVE
C*********************************************************************
C
C  This subroutine defines the terms a_pq and b_pq for backscatter.
C
C********************************************************************
      REAL*4 K1,THI,PHI,THS,PHS
      COMPLEX K2,APQ(2,2),BPQ(2,2),CPQ(2,2)
C
      COMPLEX STHT,CTHT
      COMPLEX C1,C2,C3,C4,C5,C6,C7
      COMPLEX K2STHI,K2CTHI,K1STHT,K2STHT,K1CTHT,K2CTHT
      COMPLEX RPER0,RPER1,RPAR0,RPAR1   
      COMPLEX ZVV,ZVH,ZHV,ZHH
      REAL*4 K1CTHI,K1STHI,STHI,CTHI,STHS,CTHS
      REAL*4 CPHSI,SPHSI,CPHI,SPHI
      REAL*4 C8,C9,C10
C
      STHI=SIN(THI)
      CTHI=COS(THI)
      STHS=SIN(THS)
      CTHS=COS(THS)
      SPHSI=SIN(PHS-PHI)
      CPHSI=COS(PHS-PHI)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      STHT=K1*STHI/K2
      CTHT=CSQRT(1.0-STHT*STHT)
C
      K1CTHI=K1*CTHI
      K2CTHI=K2*CTHI
      K1STHI=K1*STHI
      K2STHI=K2*STHI
      K2CTHT=K2*CTHT
      K1CTHT=K1*CTHT
      K2STHT=K2*STHT
      K1STHT=K1*STHT
C
      C1=K1CTHI+K2CTHT
      C2=K1CTHI-K2CTHT
      C3=K1STHI+K2STHT
      C4=K2CTHI+K1CTHT
      C5=K2CTHI-K1CTHT
      C6=K2STHI-K1STHT
      C7=K2STHI+K1STHT
C
      RPER0=C1/C2
c      RPER1=RPER0*C3/C2
      RPAR0=C4/C5
c      RPAR1=-(C6-RPAR0*C7)/C5
C
      RPAR1=CMPLX(0.0,0.0)
      RPER1=CMPLX(0.0,0.0)
C       
      C8=CTHI-CTHS
      C9=1.0-CTHI*CTHS
      C10=STHS-STHI*CPHSI
C
      APQ(1,1)= RPAR0*C8*CPHSI
      APQ(1,2)=-RPER0*C9*SPHSI
      APQ(2,1)= RPAR0*C9*SPHSI
      APQ(2,2)= RPER0*C8*CPHSI
C
      ZVV= RPAR0*C10+RPAR1*C8*CPHSI
      ZVH=-SPHSI*(RPER0*STHI*CTHS+RPER1*C9)
      ZHV= SPHSI*(RPAR0*STHI*CTHS+RPAR1*C9)
      ZHH= RPER0*C10+RPER1*C8*CPHSI
C
      BPQ(1,1)=ZVV*CPHI
      BPQ(1,2)=ZVH*CPHI
      BPQ(2,1)=ZHV*CPHI
      BPQ(2,2)=ZHH*CPHI
C
      CPQ(1,1)=ZVV*SPHI
      CPQ(1,2)=ZVH*SPHI
      CPQ(2,1)=ZHV*SPHI
      CPQ(2,2)=ZHH*SPHI
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE INIS_GAUSS(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,IN,IS)
      SAVE
C*********************************************************************
C
C  Subroutine to numerically compute the IN and IS terms.  This 
C  routine is to be used only for Gaussian distributions of the
C  surface height and Gaussian correlation function.  Therefore, 
C  it uses truncated summations instead of numerical integration.
C
C*********************************************************************
      REAL*4 QX,QY,QZ,QT,S,L
      COMPLEX APQ(2,2),BPQ(2,2),CPQ(2,2),IN(2,2,2,2),IS(2,2,2,2)
C
      INTEGER I1,I2,I3,I4,N,NMAX
      REAL*4 SUM,TERM,XINT,ADD,CHK,L2,TN,TS,SUMN,SUMSX,SUMSY
      REAL*4 PREC,PI,XN
      LOGICAL TEST
      DATA NMAX/100/
c      DATA PREC/1.0E-5/
      DATA PREC/1.0E-6/
      DATA PI/3.141592654/
C
C  Do summation for IN and IS (using Gaussian correlation function).
C
      SUM=0.0
      TERM=1.0
      TEST=.FALSE.
      DO N=1,NMAX
        CHK=ABS(SUM)
        XN=FLOAT(N)
        TERM=(TERM/XN)*(QZ*S)**2
        XINT=EXP(-(QT*L)**2/(4.0*XN))/XN
        ADD=TERM*XINT
        SUM=SUM+ADD
        IF (ABS(ADD).LE.(PREC*CHK)) THEN
          IF (TEST) GOTO 10
        ELSE
          TEST=.TRUE.
          IF (N.EQ.NMAX) PAUSE 'Not enough terms in SUM for IN'
        ENDIF  
      ENDDO
C           
   10 L2=L*L
      TN=PI*L2*EXP(-(QZ*S)**2)
      TS=-TN/QZ
      SUMN=TN*SUM
      SUMSX=TS*QX*SUM
      SUMSY=TS*QY*SUM
C
C  Multiply by polarization dependent coefficients and
C  set up matrix.  I1=P, I2=Q, I3=M, I4=N
C
      DO I1=1,2
        DO I2=1,2
          DO I3=1,2
            DO I4=1,2
              IN(I1,I2,I3,I4)=SUMN*APQ(I1,I2)*CONJG(APQ(I3,I4))
              IS(I1,I2,I3,I4)=SUMSX*(BPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(BPQ(I3,I4)))
     &                       +SUMSY*(CPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(CPQ(I3,I4)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE INIS_CORR3(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,IN,IS)
      SAVE
C*********************************************************************
C
C  Subroutine to numerically compute the IN and IS terms.  This 
C  routine is to be used only for Gaussian distributions of the
C  surface height and Gaussian correlation function.  Therefore, 
C  it uses truncated summations instead of numerical integration.
C
C*********************************************************************
      REAL*4 QX,QY,QZ,QT,S,L
      COMPLEX APQ(2,2),BPQ(2,2),CPQ(2,2),IN(2,2,2,2),IS(2,2,2,2)
C
      INTEGER I1,I2,I3,I4,N,NMAX,NNMAX,NE,NO,NCALC
      PARAMETER(NMAX=100,NNMAX=NMAX+5*NMAX/10)
      REAL*4 SUM,TERM,XINT,ADD,CHK,L2,TN,TS,SUMN,SUMSX,SUMSY
      REAL*4 PREC,PI,XN,ARG,OINT(NNMAX),EINT(NNMAX),MU,mGAMMLN
      LOGICAL TEST,EVEN
      DATA PREC/1.0E-5/
      DATA PI/3.141592654/
C
C  Do summation for IN and IS (using corr #3 in Polarimetry book)
C
      SUM=0.0
      TERM=1.0
      TEST=.FALSE.
C
C  First try for N=1,10 modified bessel function evaluations, then
C  if it does not converge go to 20 (then 30, 40, 50, etc.)
C
      ARG=L*QT
      DO N=1,NMAX
        IF (MOD(N,10).EQ.1) THEN
          NE=15*(1+N/10)
          NO=13+15*(N/10)
          CALL RKBESL(ARG,0.0,NE,1,EINT,NCALC)
          NE=NCALC
C          IF (NCALC.NE.NE) STOP 'Error in RKBESL'
          CALL RKBESL(ARG,0.5,NO,1,OINT,NCALC)
          NO=NCALC
C          IF (NCALC.NE.NO) STOP 'Error in RKBESL'
        END IF
        XN=FLOAT(N)
        MU=1.5*XN-1.0
        IF (EVEN(N)) THEN
          IF (INT(MU+1.0).GT.NE) PAUSE 'Error in even RKBESL'
          XINT=EINT(INT(MU+1.0))
        ELSE
          IF (INT(MU+0.5).GT.NO) PAUSE 'Error in odd RKBESL'
          XINT=OINT(INT(MU+0.5))
        ENDIF
        XINT=XINT*((L*QT)**MU)*(2.0**(1.0-MU))/EXP(mGAMMLN(MU+1.0))
        CHK=ABS(SUM)
        TERM=(TERM/XN)*(QZ*S)**2
        ADD=TERM*XINT
        SUM=SUM+ADD
        IF (ABS(ADD).LE.(PREC*CHK)) THEN
          IF (TEST) GOTO 10
        ELSE
          TEST=.TRUE.
          IF (N.EQ.NMAX) PAUSE 'Not enough terms in SUM for IN'
        ENDIF  
      ENDDO
C           
   10 L2=L*L
      TN=PI*L2*EXP(-(QZ*S)**2)
      TS=-TN/QZ
      SUMN=TN*SUM
      SUMSX=TS*QX*SUM
      SUMSY=TS*QY*SUM
C
C  Multiply by polarization dependent coefficients and
C  set up matrix.  I1=P, I2=Q, I3=M, I4=N
C
      DO I1=1,2
        DO I2=1,2
          DO I3=1,2
            DO I4=1,2
              IN(I1,I2,I3,I4)=SUMN*APQ(I1,I2)*CONJG(APQ(I3,I4))
              IS(I1,I2,I3,I4)=SUMSX*(BPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(BPQ(I3,I4)))
     &                       +SUMSY*(CPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(CPQ(I3,I4)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE INIS_ARB(APQ,BPQ,CPQ,QX,QY,QZ,QT,S,L,LS,IN,IS)
      SAVE
C*********************************************************************
C
C  Subroutine to numerically compute the IN and IS terms.  This 
C  routine is to be used for any arbitrary correlation function.
C
C*********************************************************************
      REAL*4 QX,QY,QZ,QT,S,L,LS
      COMPLEX APQ(2,2),BPQ(2,2),CPQ(2,2),IN(2,2,2,2),IS(2,2,2,2)
C
      INTEGER I1,I2,I3,I4,N,NMAX,XN,IFLAG
      REAL*4 ADD,CHK,TN,TS,SUMN,SUMS,PREC,PI,XL,XLS,XQT,XQZ,XS,TMP
      REAL*4 SUMSX,SUMSY,FACTRL,INFINITY,MIDTMP,MID,TMP1
      LOGICAL TEST
      COMMON /CORR1/ XL,XLS,XQT,XQZ,XS,XN
      EXTERNAL CORRI,CORRS,MIDEXP,MIDPNT
      DATA NMAX/100/
      DATA PREC/1.0E-5/
      DATA PI/3.141592654/
      DATA INFINITY/1.0E10/
C
C      WRITE(*,*) 'ENTER INIS_ARB'
      XL=L
      XLS=LS
      XQT=QT
      XQZ=QZ
      XS=S
C
C  Do integration for IN and IS (using arbitrary correlation 
C  function).
C
      SUMN=0.0
      TEST=.FALSE.
      MIDTMP=(0.5*L*L)*(1.0+SQRT(1.0+4.0*(LS/L)**4))
      DO N=1,NMAX
        XN=N
        CHK=ABS(SUMN)
C**********************************************************************
C  Different ways to compute the integral over integrand CORRI        
C**********************************************************************
C        CALL QROMO(CORRI,0.0,INFINITY,TMP,MIDEXP)                    
C----------------------------------------------------------------------
        MID=10.0*SQRT(MIDTMP/FLOAT(XN**2))                             
        CALL QROMO(CORRI,0.0,MID,TMP,MIDPNT)                          
C**********************************************************************
        ADD=(TMP*(QZ*S)**(2*N))/FACTRL(N)
        SUMN=SUMN+ADD
        IF (ABS(ADD).LE.(PREC*CHK)) THEN
          IF (TEST) GOTO 10
        ELSE
          TEST=.TRUE.
          IF (N.EQ.NMAX) PAUSE 'Not enough terms in SUM for IN'
        ENDIF  
      ENDDO
   10 CONTINUE
C**********************************************************************
C  Different ways to compute the integral over integrand CORRS        
C**********************************************************************
C      CALL QROMO(CORRS,0.0,INFINITY,SUMS,MIDEXP)                     
C----------------------------------------------------------------------
      MID=10.0*SQRT(MIDTMP)                                            
      CALL QROMO(CORRS,0.0,MID,SUMS,MIDPNT)                           
C**********************************************************************
      TN=2.0*PI*EXP(-(QZ*S)**2)
      TS=TN*QZ*S*S/QT
      SUMN=TN*SUMN
      SUMSX=TS*QX*SUMS
      SUMSY=TS*QY*SUMS
C
C  Multiply by polarization dependent coefficients and
C  set up matrix.  I1=P, I2=Q, I3=M, I4=N
C
      DO I1=1,2
        DO I2=1,2
          DO I3=1,2
            DO I4=1,2
              IN(I1,I2,I3,I4)=SUMN*APQ(I1,I2)*CONJG(APQ(I3,I4))
              IS(I1,I2,I3,I4)=SUMSX*(BPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(BPQ(I3,I4)))
     &                       +SUMSY*(CPQ(I1,I2)*CONJG(APQ(I3,I4))
     &                              +APQ(I1,I2)*CONJG(CPQ(I3,I4)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C      WRITE(*,*) 'EXIT INIS_ARB'
      RETURN
      END
CC*********************************************************************
CC*********************************************************************
C      SUBROUTINE CORR(X,I,J,RHO,RHOX)
C      SAVE
CC*********************************************************************
CC
CC  This subroutine returns the correlation function and its
CC  derivative for a given rough surface.
CC
CC  INPUTS:
CC
CC    X - Shift from zero
CC    If I=1 then compute RHO
CC    If J=1 then compute RHOX
CC
CC  OUTPUTS:
CC
CC    RHO  - Correlation function
CC    RHOX - Derivative of correlation function
CC
CC  COMMON VARS:
CC
CC    XL  - Large scale correlation length of the surface
CC    XLS - Small scale correlation length of the surface
CC
CC********************************************************************
C      INTEGER I,J,XN
C      REAL*4 XL,XLS,X,RHO,RHOX,XQT,XQZ,XS,TMP
C      COMMON /CORR1/ XL,XLS,XQT,XQZ,XS,XN
CC
C      TMP=XLS**4+(X*XL)**2
C      RHO=EXP(-X*X/TMP**0.5)
C      IF (J.EQ.1) RHOX=-RHO*X*(TMP+XLS**4)/TMP**1.5
CC
C      RETURN
C      END
C*********************************************************************
C*********************************************************************
      FUNCTION CORRI(X)
      SAVE
C*********************************************************************
C
C  This fuction is the integrand for the integral to compute
C  the incoherent intensity for a rough surface with arbitrary
C  correlation function.
C
C  INPUTS:
C
C    X - Shift from zero
C
C  COMMON VARS:
C
C    XL  - Large scale correlation length of the surface
C    XLS - Small scale correlation length of the surface
C
C********************************************************************
      INTEGER XN
      REAL*4 XL,XLS,X,XQT,XQZ,XS,DUM,RHO,CORRI,BESSJ0
      COMMON /CORR1/ XL,XLS,XQT,XQZ,XS,XN
C
      CALL CORR(X,1,0,RHO,DUM)
      CORRI=(RHO**XN)*BESSJ0(XQT*X)*X
C      CORRI=COS(2.0*X)*EXP(-X)
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      FUNCTION CORRS(X)
      SAVE
C*********************************************************************
C
C  This fuction is the integrand for the integral to compute
C  the surface slope intensity for a rough surface with arbitrary
C  correlation function.
C
C  INPUTS:
C
C    X - Shift from zero
C
C  COMMON VARS:
C
C    XL  - Large scale correlation length of the surface
C    XLS - Small scale correlation length of the surface
C
C********************************************************************
      INTEGER XN
      REAL*4 XL,XLS,X,XQT,XQZ,XS,RHO,RHOX,CORRS,BESSJ1
      COMMON /CORR1/ XL,XLS,XQT,XQZ,XS,XN
C
      CALL CORR(X,1,1,RHO,RHOX)
      CORRS=RHOX*BESSJ1(XQT*X)*X*EXP(RHO*(XQZ*XS)**2)
C
      RETURN
      END
C*********************************************************************
C
C  SMALL PERTURBATION SUBROUTINES
C
C*********************************************************************
      SUBROUTINE GND_SP(EPSR,WVL,THI,PHI,THS,PHS,S,L,W,WW,IMM,MM)
      SAVE
C*********************************************************************
C
C  Program to compute the transformation matrix relating the
C  incident and scattered intensity vectors for a rough surface
C  using the SMALL PERTURBATION METHOD.
C  The incident medium is free space.
C
C  INPUTS
C    
C    EPSR    - Complex dielectric constant of the ground surface
C                relative to the medium above the surface.
C    WVL     - Wavelength in free space.
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    S       - Standard deviation of surface heights
C    L       - Surface correlation length
C    W       - Function representing the Fourier transform of the
C                correlation function of the surface.
C                (function of K1X,K1Y,K1XI,K1YI,L)
C    IMM     - IMM=0 means do not compute MM
C              IMM=1 means compute MM
C
C  OUTPUT
C
C    WW  - 4x4 complex W matrix 
C    MM  - 4X4 real modified Mueller matrix relating the incident and 
C            scattered intensity vectors. (without range dependence)
C
C*********************************************************************
C
      INTEGER IMM
      REAL*4 WVL,THI,PHI,THS,PHS,S,L,W,MM(4,4)
      COMPLEX EPSR,WW(4,4)
      EXTERNAL W
C
      REAL*4 K1,PI,F0,WCOR,THIP,CTHS
      COMPLEX FPQ(2,2),K2
      DATA PI/3.141592654/
C
C  Initialize variables
C
      K1=2.0*PI/WVL
      K2=K1*CSQRT(EPSR)
      CTHS=COS(THS)
C
C  Compute the polarization dependent coefficients and 
C  the Fourier transform of the correlation function.
C
      CALL COEF_SP(K1,K2,THI,PHI,THS,PHS,L,S,W,FPQ,WCOR)
C
C  Compute the W matrix.
C
      F0=K1*K1*CTHS*CTHS*WCOR
C
      WW(1,1)=F0*FPQ(1,1)*CONJG(FPQ(1,1))
      WW(1,2)=F0*FPQ(1,2)*CONJG(FPQ(1,2))
      WW(1,3)=F0*FPQ(1,2)*CONJG(FPQ(1,1))
      WW(1,4)=F0*FPQ(1,1)*CONJG(FPQ(1,2))
      WW(2,1)=F0*FPQ(2,1)*CONJG(FPQ(2,1))
      WW(2,2)=F0*FPQ(2,2)*CONJG(FPQ(2,2))
      WW(2,3)=F0*FPQ(2,2)*CONJG(FPQ(2,1))
      WW(2,4)=F0*FPQ(2,1)*CONJG(FPQ(2,2))
      WW(3,1)=F0*FPQ(2,1)*CONJG(FPQ(1,1))
      WW(3,2)=F0*FPQ(2,2)*CONJG(FPQ(1,2))
      WW(3,3)=F0*FPQ(2,2)*CONJG(FPQ(1,1))
      WW(3,4)=F0*FPQ(2,1)*CONJG(FPQ(1,2))
      WW(4,1)=F0*FPQ(1,1)*CONJG(FPQ(2,1))
      WW(4,2)=F0*FPQ(1,2)*CONJG(FPQ(2,2))
      WW(4,3)=F0*FPQ(1,2)*CONJG(FPQ(2,1))
      WW(4,4)=F0*FPQ(1,1)*CONJG(FPQ(2,2))
C
C  Compute the average matrix relating the incident and
C  scattered intensities for the rough surface.
C
      IF (IMM.EQ.1) THEN
        CALL W_TO_MM(WW,MM)
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      SUBROUTINE COEF_SP(K1,K2,THI,PHI,THS,PHS,L,S,W,FPQ,WCOR)
      SAVE
C*********************************************************************
C
C  This subroutine computes the coefficients for FPQ for 
C  polarizations VV, VH, HV, HH.
C
C  INPUTS:
C
C    K1    - wave number in free space
C    K2    - wave number in dielectric
C    THI     - Incident elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHI     - Incident azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    THS     - Scattered elevation angle in radians (measured from z-axis 
C                to the unit vector in the direction of 
C                propagation of the incident wave)
C    PHS     - Scattered azimuth angle in radians (measured from the
C                x-axis to the projection on the xy-plane of the unit
C                vector in the direction of propagation of the incident
C                wave)
C    L       - Surface correlation length
C    S       - Standard deviation of surface heights
C    W       - Function representing the Fourier transform of the
C                correlation function of the surface.
C                (function of K1X,K1Y,K1XI,K1YI,L)
C
C  OUTPUTS:
C
C    FPQ  - Polarization dependent coefficients
C    WCOR - Fourier transform of the correlation function
C
C********************************************************************
      REAL*4 K1,THI,PHI,THS,PHS,L,W,WCOR,S
      COMPLEX K2,FPQ(2,2)
C
      REAL*4 STHI,SPHI,STHS,SPHS,CTHI,CPHI,CTHS,CPHS
      REAL*4 K1XI,K1YI,K1ZI,K1X,K1Y,K1Z,K1P,K1PI
      REAL*4 KK1,KK2
      COMPLEX K2ZI,K2Z,C1,C2,C3,C4,C5,C6
C
      STHI=SIN(THI)
      SPHI=SIN(PHI)
      STHS=SIN(THS)
      SPHS=SIN(PHS)
      CTHI=COS(THI)
      CPHI=COS(PHI)
      CTHS=COS(THS)
      CPHS=COS(PHS)
C
      K1XI=K1*STHI*CPHI
      K1YI=K1*STHI*SPHI
      K1ZI=K1*CTHI
      K2ZI=CSQRT(K2**2-(K1*STHI)**2)
C
      K1X=K1*STHS*CPHS
      K1Y=K1*STHS*SPHS
      K1Z=K1*CTHS
      K2Z=CSQRT(K2**2-(K1*STHS)**2)
C
      K1P=K1*STHS
      K1PI=K1*STHI
C
      KK1=(K1X*K1XI+K1Y*K1YI)/(K1P*K1PI)
      KK2=(K1Y*K1XI-K1X*K1YI)/(K1P*K1PI)
      C1=K2*K2-K1*K1
      C2=K1*K1*K2Z+K2*K2*K1Z
      C3=K1*K1*K2ZI-K2*K2*K1ZI
      C4=K2ZI-K1ZI
      C5=K2Z+K1Z
      C6=K1P*K1PI-(K1*K1/K2**2)*K2Z*K2ZI*KK1
      FPQ(1,1)=2.0*K2*K2*K1ZI*C1*C6/(C2*C3)
      FPQ(1,2)=2.0*K1*K2Z*K1ZI*C1*KK2/(C2*C4)
      FPQ(2,1)=-2.0*K1*K2ZI*K1ZI*C1*KK2/(C5*C3)
      FPQ(2,2)=2.0*K1ZI*C1*KK1/(C5*C4)
C
      WCOR=W(K1X,K1Y,K1XI,K1YI,L,S)
C
      RETURN
      END
CC*********************************************************************
CC*********************************************************************
C      FUNCTION W(K1X,K1Y,K1XI,K1YI,L,S)
C      SAVE
CC*********************************************************************
CC
CC  This fuction is the Fourier transform of the correlation 
CC  function of the rough surface (i.e. the normalized roughness 
CC  spectrum).
CC
CC  THIS PARTICULAR FUNCTIONAL FORM IS FOR A GAUSSIAN CORRELATION
CC  FUNCTION.
CC
CC  INPUTS:
CC 
CC    K1X  - Wave number in x-direction (scattered wave)
CC    K1Y  - Wave number in y-direction (scattered wave)
CC    K1XI - Wave number in x-direction (incident wave)
CC    K1YI - Wave number in y-direction (incident wave)
CC    L    - Correlation length of the surface
CC    S    - Standard deviation of surface heights
CC
CC********************************************************************
C      REAL*4 L,W,K1X,K1Y,K1XI,K1YI,KDP2,S,PI
C      DATA PI/3.141592654/
CC
C      KDP2=(K1X-K1XI)**2+(K1Y-K1YI)**2
C      W=(L*L*S*S/(4.0*PI))*EXP(-0.25*KDP2*L**2)
CC
C      RETURN
C      END
C**********************************************************************
C*********************************************************************
      SUBROUTINE CORR(X,I,J,RHO,RHOX)
      SAVE
C*********************************************************************
C
C  This subroutine returns the correlation function and its
C  derivative for a given rough surface.
C
C  INPUTS:
C
C    X - Shift from zero
C    If I=1 then compute RHO
C    If J=1 then compute RHOX
C
C  OUTPUTS:
C
C    RHO  - Correlation function
C    RHOX - Derivative of correlation function
C
C  COMMON VARS:
C
C    XL  - Large scale correlation length of the surface
C    XLS - Small scale correlation length of the surface
C
C********************************************************************
      INTEGER I,J,XN
      REAL*4 XL,XLS,X,RHO,RHOX,XQT,XQZ,XS,TMP
      COMMON /CORR1/ XL,XLS,XQT,XQZ,XS,XN
C
      TMP=XLS**4+(X*XL)**2
      RHO=EXP(-X*X/TMP**0.5)
      IF (J.EQ.1) RHOX=-RHO*X*(TMP+XLS**4)/TMP**1.5
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
      FUNCTION PROB(HX,HY,M)
      SAVE
C*********************************************************************
C
C  This function returns the value of the surface slope probability
C  density at the X and Y slopes HX and HY, respectively.
C
C  INPUTS:
C
C    HX - Slope in the X-direction
C    HY - Slope in the Y-direction
C    M  - RMS surface slope
C
C********************************************************************
      REAL*4 HX,HY,PROB,PI,M
      DATA PI/3.141592654/
C
      PROB=EXP(-(HX*HX+HY*HY)/(2.0*M*M))
      PROB=PROB/(2.0*PI*M*M)
C
      RETURN
      END
CC*********************************************************************
C**********************************************************************
      SUBROUTINE W_TO_MM(W,MM)
      SAVE
C**********************************************************************
C
C   Subroutine to compute the real 4x4 Modified Mueller Matrix
C   (see page 27 of Radar Polarimetry for Geoscience Applications) 
C   from a given complex 4x4 W Matrix.  Both matrices are assumed
C   to be in the same convention (either BSA or FSA).
C
C   CALL PARAMETERS:
C
C   INPUT:
C
C   W     - Complex 4x4 W Matrix
C
C   OUTPUT:
C
C   MM    - Real 4x4 Modified Mueller Matrix
C
C**********************************************************************
C
      INTEGER I,J
      REAL*4 MM(4,4)
      COMPLEX*8 W(4,4),TMP(4,4),TMP1(4,4),V(4,4),VI(4,4)
      DATA ((V(I,J),J=1,4),I=1,4)/
     &            (1.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(1.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(0.0,0.0),(1.0,0.0),(1.0,0.0),
     &            (0.0,0.0),(0.0,0.0),(0.0,-1.0),(0.0,1.0)/
      DATA ((VI(I,J),J=1,4),I=1,4)/
     &            (1.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(1.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(0.0,0.0),(0.5,0.0),(0.0,0.5),
     &            (0.0,0.0),(0.0,0.0),(0.5,0.0),(0.0,-0.5)/
C
C   Zero out the MM Matrix array.
C
      DO I=1,4
        DO J=1,4
          MM(I,J)=0.0
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C   Compute MM Matrix
C----------------------------------------------------------------------
      CALL CCMULT(4,4,4,V,W,TMP)
      CALL CCMULT(4,4,4,TMP,VI,TMP1)
C
      DO I=1,4
        DO J=1,4
          MM(I,J)=REAL(TMP1(I,J))
        ENDDO
      ENDDO
C
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE MM_TO_MSSO(MM,MSSO)
      SAVE
C**********************************************************************
C
C   Subroutine to compute the real 4x4 Stokes Scattering Operator
C   (see page 30 of Radar Polarimetry for Geoscience Applications) 
C   from a given real 4x4 Modified Mueller Matrix.  Both matrices 
C   are assumed to be in the same convention (either BSA or FSA).
C
C   CALL PARAMETERS:
C
C   MM    - Real 4x4 Modified Mueller Matrix
C   MSSO  - Real 4x4 Stokes Scattering Operator
C
C**********************************************************************
C
      INTEGER I,J
      REAL*4 MM(4,4),MSSO(4,4)
C
C   Zero out the Stokes Scattering Operator array.
C
      DO I=1,4
        DO J=1,4
          MSSO(I,J)=0.0
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C   Compute SSO Matrix
C----------------------------------------------------------------------
      MSSO(1,1)=MM(1,1)
      MSSO(1,2)=MM(1,2)
      MSSO(1,3)=MM(1,3)
      MSSO(1,4)=MM(1,4)
      MSSO(2,1)=MM(2,1)
      MSSO(2,2)=MM(2,2)
      MSSO(2,3)=MM(2,3)
      MSSO(2,4)=MM(2,4)
      MSSO(3,1)=0.5*MM(3,1)
      MSSO(3,2)=0.5*MM(3,2)
      MSSO(3,3)=0.5*MM(3,3)
      MSSO(3,4)=0.5*MM(3,4)
      MSSO(4,1)=-0.5*MM(4,1)
      MSSO(4,2)=-0.5*MM(4,2)
      MSSO(4,3)=-0.5*MM(4,3)
      MSSO(4,4)=-0.5*MM(4,4)
C
      RETURN
      END
C**********************************************************************
C*********************************************************************
      FUNCTION WFTN(K1X,K1Y,K1XI,K1YI,L,S)
      SAVE
C*********************************************************************
C
C  This fuction is the Fourier transform of the correlation 
C  function of the rough surface (i.e. the normalized roughness 
C  spectrum).
C
C  THIS PARTICULAR FUNCTIONAL FORM IS FOR A GAUSSIAN CORRELATION
C  FUNCTION.
C
C  INPUTS:
C 
C    K1X  - Wave number in x-direction (scattered wave)
C    K1Y  - Wave number in y-direction (scattered wave)
C    K1XI - Wave number in x-direction (incident wave)
C    K1YI - Wave number in y-direction (incident wave)
C    L    - Correlation length of the surface
C    S    - Standard deviation of surface heights
C
C********************************************************************
      REAL*4 L,WFTN,K1X,K1Y,K1XI,K1YI,KDP2,S,PI
      DATA PI/3.141592654/
C
      KDP2=(K1X-K1XI)**2+(K1Y-K1YI)**2
      WFTN=(L*L*S*S/(4.0*PI))*EXP(-0.25*KDP2*L**2)
C
      RETURN
      END
C*********************************************************************
C**********************************************************************
C**********************************************************************
      FUNCTION EVEN(N)
C**********************************************************************
      INTEGER N
      LOGICAL EVEN
C
      IF (FLOAT(N/2).EQ.FLOAT(N)/2.0) THEN
        EVEN=.TRUE.
      ELSE
        EVEN=.FALSE.
      ENDIF
C
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE CCMULT(L,M,N,A,B,C)
      SAVE
C***************************************************************************
C
C  Subroutine to perform matrix multiplication for complex
C  matrices.
C
C  INPUTS:
C
C  L - Array size specifier
C  M - Array size specifier
C  R - Array size specifier
C  A - L x M complex array
C  B - M x N complex array
C
C  OUTPUTS:
C
C  C - L x N complex array ( C = A B )
C
C***************************************************************************
      INTEGER L,M,N,I,J,K
      COMPLEX*8 A(L,M),B(M,N),C(L,N)
C
      DO I=1,L
        DO J=1,N
          C(I,J)=CMPLX(0.0,0.0)
          DO K=1,M
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C***************************************************************************
C**********************************************************************
      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE)
      SAVE
C**********************************************************************
      INTEGER JMAX,JMAXP,KM,K,J
      REAL*4 EPS,A,B,SS,DSS
      PARAMETER (EPS=1.E-4,JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      REAL*4 S(JMAXP),H(JMAXP)
      EXTERNAL FUNC
C      CALL INTEG_PLOT(FUNC,A,B)
      H(1)=1.
      DO 11 J=1,JMAX
        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0,SS,DSS)
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      PAUSE 'Too many steps.'
      END
C**********************************************************************
C**********************************************************************
      FUNCTION FACTRL(N)
      SAVE
C**********************************************************************
      INTEGER NTOP,N,J
      REAL*4 A(33),FACTRL,mGAMMLN
      DATA NTOP,A(1)/0,1./
      IF (N.LT.0) THEN
        PAUSE 'negative factorial'
      ELSE IF (N.LE.NTOP) THEN
        FACTRL=A(N+1)
      ELSE IF (N.LE.32) THEN
        DO 11 J=NTOP+1,N
          A(J+1)=J*A(J)
11      CONTINUE
        NTOP=N
        FACTRL=A(N+1)
      ELSE
        FACTRL=EXP(mGAMMLN(N+1.))
      ENDIF
      RETURN
      END
C**********************************************************************
C**********************************************************************
      FUNCTION mGAMMLN(XX)
      SAVE
C**********************************************************************
      INTEGER J
      REAL*4 XX,mGAMMLN
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,mikeMYDBLE
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=mikeMYDBLE(XX)-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      mGAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C**********************************************************************
C**********************************************************************
      FUNCTION BESSJ0(X)
      SAVE
C**********************************************************************
      REAL*4 BESSJ0,X,AX,Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D
     *0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-.785398164
        BESSJ0=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
C**********************************************************************
C**********************************************************************
      FUNCTION BESSJ1(X)
      SAVE
C**********************************************************************
      REAL*4 BESSJ1,X,AX,Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-2.356194491
        BESSJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *SIGN(1.,X)
      ENDIF
      RETURN
      END
C**********************************************************************
C**********************************************************************
      FUNCTION mikeMYDBLE(A)
      SAVE
C**********************************************************************
      REAL*8 mikeMYDBLE
      REAL*4 A
      CHARACTER*30 CHARA
      WRITE(CHARA,10) A
 10   FORMAT(E30.17)
      READ(CHARA,10) mikeMYDBLE
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE RKBESL(X,ALPHA,NB,IZE,BK,NCALC)
C**********************************************************************
C
C  This FORTRAN 77 routine calculates modified Bessel functions
C  of the second kind, K SUB(N+ALPHA) (X), for non-negative
C  argument X, and non-negative order N+ALPHA, with or without
C  exponential scaling.
C
C  Explanation of variables in the calling sequence
C
C  Description of output values ..
C
C X     - Working precision non-negative real argument for which
C         K's or exponentially scaled K's (K*EXP(X))
C         are to be calculated.  If K's are to be calculated,
C         X must not be greater than XMAX (see below).
C ALPHA - Working precision fractional part of order for which 
C         K's or exponentially scaled K's (K*EXP(X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
C         and 2 if exponentially scaled K's are to be calculated.
C BK    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BK
C         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
C         or the corresponding exponentially scaled functions.
C         If (0 .LT. NCALC .LT. NB), BK(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BK, the user should check that 
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = The smallest positive floating-point number such that 
C            1.0+EPS .GT. 1.0
C   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution 
C            to equation:
C               W(X) * (1-1/8X+9/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C   SQXMIN = Square root of beta**minexp
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMIN   = Smallest positive machine number; approximately
C            beta**minexp
C
C
C     Approximate values for some important machines are:
C
C                          beta       minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15
C  Cyber 180/185 
C    under NOS   (S.P.)      2         -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16
C  IBM 3033      (D.P.)     16          -65          63    2.22D-16
C  VAX           (S.P.)      2         -128         127    5.96E-8
C  VAX D-Format  (D.P.)      2         -128         127    1.39D-17
C  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16
C
C
C                         SQXMIN       XINF        XMIN      XMAX
C
C CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858
C Cyber 180/855
C   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342
C IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852
C VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715
C VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715
C VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all K's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, the B-vector is not calculated,
C       and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
C       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
C       the B-vector is not calculated.  Note that again 
C       NCALC .NE. NB.
C
C  0 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BK(I) contains correct function
C       values for I .LE. NCALC, and contains the ratios
C       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by J. B. Campbell
C  (2) that computes values of the Bessel functions K of real
C  argument and real order.  Modifications include the addition
C  of non-scaled functions, parameterization of machine
C  dependencies, and the use of more accurate approximations
C  for SINH and SIN.
C
C References: "On Temme's Algorithm for the Modified Bessel
C              Functions of the Third Kind," Campbell, J. B.,
C              TOMS 6(4), Dec. 1980, pp. 581-586.
C
C             "A FORTRAN IV Subroutine for the Modified Bessel
C              Functions of the Third Kind of Real Order and Real
C              Argument," Campbell, J. B., Report NRC/ERB-925,
C              National Research Council, Canada.
C
C  Latest modification: October 19, 1990
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER I,IEND,ITEMP,IZE,J,K,M,MPLUS1,NB,NCALC
      REAL
     1    A,ALPHA,BLPHA,BK,BK1,BK2,C,D,DM,D1,D2,D3,ENU,EPS,ESTF,ESTM,
     2    EX,FOUR,F0,F1,F2,HALF,ONE,P,P0,Q,Q0,R,RATIO,S,SQXMIN,T,TINYX,
     3    TWO,TWONU,TWOX,T1,T2,WMINF,X,XINF,XMAX,XMIN,X2BY4,ZERO
      DIMENSION BK(NB),P(8),Q(7),R(5),S(4),T(6),ESTM(6),ESTF(7)
C---------------------------------------------------------------------
C  Mathematical constants
C    A = LOG(2.D0) - Euler's constant
C    D = SQRT(2.D0/PI)
C---------------------------------------------------------------------
      DATA HALF,ONE,TWO,ZERO/0.5E0,1.0E0,2.0E0,0.0E0/
      DATA FOUR,TINYX/4.0E0,1.0E-10/
      DATA A/ 0.11593151565841244881E0/,D/0.797884560802865364E0/
CD    DATA HALF,ONE,TWO,ZERO/0.5D0,1.0D0,2.0D0,0.0D0/
CD    DATA FOUR,TINYX/4.0D0,1.0D-10/
CD    DATA A/ 0.11593151565841244881D0/,D/0.797884560802865364D0/
C---------------------------------------------------------------------
C  Machine dependent parameters
C---------------------------------------------------------------------
      DATA EPS/1.19E-7/,SQXMIN/1.08E-19/,XINF/3.40E+38/
      DATA XMIN/1.18E-38/,XMAX/85.337E0/
CD    DATA EPS/2.22D-16/,SQXMIN/1.49D-154/,XINF/1.79D+308/
CD    DATA XMIN/2.23D-308/,XMAX/705.342D0/
C---------------------------------------------------------------------
C  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
C                                         + Euler's constant
C         Coefficients converted from hex to decimal and modified
C         by W. J. Cody, 2/26/82
C  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
C  T    - Approximation for SINH(Y)/Y
C---------------------------------------------------------------------
      DATA P/ 0.805629875690432845E00,    0.204045500205365151E02,
     1        0.157705605106676174E03,    0.536671116469207504E03,
     2        0.900382759291288778E03,    0.730923886650660393E03,
     3        0.229299301509425145E03,    0.822467033424113231E00/
      DATA Q/ 0.294601986247850434E02,    0.277577868510221208E03,
     1        0.120670325591027438E04,    0.276291444159791519E04,
     2        0.344374050506564618E04,    0.221063190113378647E04,
     3        0.572267338359892221E03/
      DATA R/-0.48672575865218401848E+0,  0.13079485869097804016E+2,
     1       -0.10196490580880537526E+3,  0.34765409106507813131E+3,
     2        0.34958981245219347820E-3/
      DATA S/-0.25579105509976461286E+2,  0.21257260432226544008E+3,
     1       -0.61069018684944109624E+3,  0.42269668805777760407E+3/
      DATA T/ 0.16125990452916363814E-9, 0.25051878502858255354E-7,
     1        0.27557319615147964774E-5, 0.19841269840928373686E-3,
     2        0.83333333333334751799E-2, 0.16666666666666666446E+0/
      DATA ESTM/5.20583E1, 5.7607E0, 2.7782E0, 1.44303E1, 1.853004E2,
     1          9.3715E0/
      DATA ESTF/4.18341E1, 7.1075E0, 6.4306E0, 4.25110E1, 1.35633E0,
     1          8.45096E1, 2.0E1/
CD    DATA P/ 0.805629875690432845D00,    0.204045500205365151D02,
CD   1        0.157705605106676174D03,    0.536671116469207504D03,
CD   2        0.900382759291288778D03,    0.730923886650660393D03,
CD   3        0.229299301509425145D03,    0.822467033424113231D00/
CD    DATA Q/ 0.294601986247850434D02,    0.277577868510221208D03,
CD   1        0.120670325591027438D04,    0.276291444159791519D04,
CD   2        0.344374050506564618D04,    0.221063190113378647D04,
CD   3        0.572267338359892221D03/
CD    DATA R/-0.48672575865218401848D+0,  0.13079485869097804016D+2,
CD   1       -0.10196490580880537526D+3,  0.34765409106507813131D+3,
CD   2        0.34958981245219347820D-3/
CD    DATA S/-0.25579105509976461286D+2,  0.21257260432226544008D+3,
CD   1       -0.61069018684944109624D+3,  0.42269668805777760407D+3/
CD    DATA T/ 0.16125990452916363814D-9, 0.25051878502858255354D-7,
CD   1        0.27557319615147964774D-5, 0.19841269840928373686D-3,
CD   2        0.83333333333334751799D-2, 0.16666666666666666446D+0/
CD    DATA ESTM/5.20583D1, 5.7607D0, 2.7782D0, 1.44303D1, 1.853004D2,
CD   1          9.3715D0/
CD    DATA ESTF/4.18341D1, 7.1075D0, 6.4306D0, 4.25110D1, 1.35633D0,
CD   1          8.45096D1, 2.0D1/
C---------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      NCALC = MIN(NB,0)-2
      IF ((NB .GT. 0) .AND. ((ENU .GE. ZERO) .AND. (ENU .LT. ONE))
     1     .AND. ((IZE .GE. 1) .AND. (IZE .LE. 2)) .AND.
     2     ((IZE .NE. 1) .OR. (EX .LE. XMAX)) .AND.
     3     (EX .GT. ZERO))  THEN
            K = 0
            IF (ENU .LT. SQXMIN) ENU = ZERO
            IF (ENU .GT. HALF) THEN
                  K = 1
                  ENU = ENU - ONE
            END IF
            TWONU = ENU+ENU
            IEND = NB+K-1
            C = ENU*ENU
            D3 = -C
            IF (EX .LE. ONE) THEN
C---------------------------------------------------------------------
C  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
C                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
C---------------------------------------------------------------------
                  D1 = ZERO
                  D2 = P(1)
                  T1 = ONE
                  T2 = Q(1)
                  DO 10 I = 2,7,2
                     D1 = C*D1+P(I)
                     D2 = C*D2+P(I+1)
                     T1 = C*T1+Q(I)
                     T2 = C*T2+Q(I+1)
   10             CONTINUE
                  D1 = ENU*D1
                  T1 = ENU*T1
                  F1 = LOG(EX)
                  F0 = A+ENU*(P(8)-ENU*(D1+D2)/(T1+T2))-F1
                  Q0 = EXP(-ENU*(A-ENU*(P(8)+ENU*(D1-D2)/(T1-T2))-F1))
                  F1 = ENU*F0
                  P0 = EXP(F1)
C---------------------------------------------------------------------
C  Calculation of F0 = 
C---------------------------------------------------------------------
                  D1 = R(5)
                  T1 = ONE
                  DO 20 I = 1,4
                     D1 = C*D1+R(I)
                     T1 = C*T1+S(I)
   20             CONTINUE
                  IF (ABS(F1) .LE. HALF) THEN
                        F1 = F1*F1
                        D2 = ZERO
                        DO 30 I = 1,6
                           D2 = F1*D2+T(I)
   30                   CONTINUE
                        D2 = F0+F0*F1*D2
                     ELSE
                        D2 = SINH(F1)/ENU
                  END IF
                  F0 = D2-ENU*D1/(T1*P0)
                  IF (EX .LE. TINYX) THEN
C--------------------------------------------------------------------
C  X.LE.1.0E-10
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                        BK(1) = F0+EX*F0
                        IF (IZE .EQ. 1) BK(1) = BK(1)-EX*BK(1)
                        RATIO = P0/F0
                        C = EX*XINF
                        IF (K .NE. 0) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
C  ALPHA .GE. 1/2
C--------------------------------------------------------------------
                              NCALC = -1
                              IF (BK(1) .GE. C/RATIO) GO TO 500
                              BK(1) = RATIO*BK(1)/EX
                              TWONU = TWONU+TWO
                              RATIO = TWONU
                        END IF
                        NCALC = 1
                        IF (NB .EQ. 1) GO TO 500
C--------------------------------------------------------------------
C  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
C--------------------------------------------------------------------
                        NCALC = -1
                        DO 80 I = 2,NB
                           IF (RATIO .GE. C) GO TO 500
                           BK(I) = RATIO/EX
                           TWONU = TWONU+TWO
                           RATIO = TWONU
   80                   CONTINUE
                        NCALC = 1
                        GO TO 420
                     ELSE
C--------------------------------------------------------------------
C  1.0E-10 .LT. X .LE. 1.0
C--------------------------------------------------------------------
                        C = ONE
                        X2BY4 = EX*EX/FOUR
                        P0 = HALF*P0
                        Q0 = HALF*Q0
                        D1 = -ONE
                        D2 = ZERO
                        BK1 = ZERO
                        BK2 = ZERO
                        F1 = F0
                        F2 = P0
  100                   D1 = D1+TWO
                        D2 = D2+ONE
                        D3 = D1+D3
                        C = X2BY4*C/D2
                        F0 = (D2*F0+P0+Q0)/D3
                        P0 = P0/(D2-ENU)
                        Q0 = Q0/(D2+ENU)
                        T1 = C*F0
                        T2 = C*(P0-D2*F0)
                        BK1 = BK1+T1
                        BK2 = BK2+T2
                        IF ((ABS(T1/(F1+BK1)) .GT. EPS) .OR.
     1                     (ABS(T2/(F2+BK2)) .GT. EPS))  GO TO 100
                        BK1 = F1+BK1
                        BK2 = TWO*(F2+BK2)/EX
                        IF (IZE .EQ. 2) THEN
                              D1 = EXP(EX)
                              BK1 = BK1*D1
                              BK2 = BK2*D1
                        END IF
                        WMINF = ESTF(1)*EX+ESTF(2)
                  END IF
               ELSE IF (EPS*EX .GT. ONE) THEN
C--------------------------------------------------------------------
C  X .GT. ONE/EPS
C--------------------------------------------------------------------
                  NCALC = NB
                  BK1 = ONE / (D*SQRT(EX))
                  DO 110 I = 1, NB
                     BK(I) = BK1
  110             CONTINUE
                  GO TO 500
               ELSE
C--------------------------------------------------------------------
C  X .GT. 1.0
C--------------------------------------------------------------------
                  TWOX = EX+EX
                  BLPHA = ZERO
                  RATIO = ZERO
                  IF (EX .LE. FOUR) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(1)/EX+ESTM(2))
                        M = INT(D2)
                        D1 = D2+D2
                        D2 = D2-HALF
                        D2 = D2*D2
                        DO 120 I = 2,M
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
  120                   CONTINUE
C--------------------------------------------------------------------
C  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
C    recurrence and K(ALPHA,X) from the wronskian
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(3)*EX+ESTM(4))
                        M = INT(D2)
                        C = ABS(ENU)
                        D3 = C+C
                        D1 = D3-ONE
                        F1 = XMIN
                        F0 = (TWO*(C+D2)/EX+HALF*EX/(C+D2+ONE))*XMIN
                        DO 130 I = 3,M
                           D2 = D2-ONE
                           F2 = (D3+D2+D2)*F0
                           BLPHA = (ONE+D1/D2)*(F2+BLPHA)
                           F2 = F2/EX+F1
                           F1 = F0
                           F0 = F2
  130                   CONTINUE
                        F1 = (D3+TWO)*F0/EX+F1
                        D1 = ZERO
                        T1 = ONE
                        DO 140 I = 1,7
                           D1 = C*D1+P(I)
                           T1 = C*T1+Q(I)
  140                   CONTINUE
                        P0 = EXP(C*(A+C*(P(8)-C*D1/T1)-LOG(EX)))/EX
                        F2 = (C+HALF-RATIO)*F1/EX
                        BK1 = P0+(D3*F0-F2+F0+BLPHA)/(F2+F1+F0)*P0
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(3)*EX+ESTF(4)
                     ELSE
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
C  recurrence, for  X .GT. 4.0
C--------------------------------------------------------------------
                        DM = AINT(ESTM(5)/EX+ESTM(6))
                        M = INT(DM)
                        D2 = DM-HALF
                        D2 = D2*D2
                        D1 = DM+DM
                        DO 160 I = 2,M
                           DM = DM-ONE
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
                           BLPHA = (RATIO+RATIO*BLPHA)/DM
  160                   CONTINUE
                        BK1 = ONE/((D+D*BLPHA)*SQRT(EX))
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(5)*(EX-ABS(EX-ESTF(7)))+ESTF(6)
                  END IF
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
C    K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                  BK2 = BK1+BK1*(ENU+HALF-RATIO)/EX
            END IF
C--------------------------------------------------------------------
C  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
C  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
C--------------------------------------------------------------------
            NCALC = NB
            BK(1) = BK1
            IF (IEND .EQ. 0) GO TO 500
            J = 2-K
            IF (J .GT. 0) BK(J) = BK2
            IF (IEND .EQ. 1) GO TO 500
            M = MIN(INT(WMINF-ENU),IEND)
            DO 190 I = 2,M
               T1 = BK1
               BK1 = BK2
               TWONU = TWONU+TWO
               IF (EX .LT. ONE) THEN
                     IF (BK1 .GE. (XINF/TWONU)*EX) GO TO 195
                     GO TO 187
                  ELSE 
                     IF (BK1/EX .GE. XINF/TWONU) GO TO 195
               END IF
  187          CONTINUE
               BK2 = TWONU/EX*BK1+T1
               ITEMP = I
               J = J+1
               IF (J .GT. 0) BK(J) = BK2
  190       CONTINUE
  195       M = ITEMP
            IF (M .EQ. IEND) GO TO 500
            RATIO = BK2/BK1
            MPLUS1 = M+1
            NCALC = -1
            DO 410 I = MPLUS1,IEND
               TWONU = TWONU+TWO
               RATIO = TWONU/EX+ONE/RATIO
               J = J+1
               IF (J .GT. 1) THEN
                     BK(J) = RATIO
                  ELSE
                     IF (BK2 .GE. XINF/RATIO) GO TO 500
                     BK2 = RATIO*BK2
               END IF
  410       CONTINUE
            NCALC = MAX(MPLUS1-K,1)
            IF (NCALC .EQ. 1) BK(1) = BK2
            IF (NB .EQ. 1) GO TO 500
  420       J = NCALC+1
            DO 430 I = J,NB
               IF (BK(NCALC) .GE. XINF/BK(I)) GO TO 500
               BK(I) = BK(NCALC)*BK(I)
               NCALC = I
  430       CONTINUE
      END IF
  500 RETURN
C---------- Last line of RKBESL ----------
      END

C**********************************************************************
      SUBROUTINE MIDPNT(FUNC,A,B,S,N)
      SAVE
C**********************************************************************
      INTEGER N,IT,J
      REAL*4 TNM,A,B,X,FUNC,SUM,DEL,DDEL,S
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      SAVE
C**********************************************************************
      INTEGER NMAX,N,I,M,NS
      PARAMETER (NMAX=10) 
      REAL*4 XA(N),YA(N),C(NMAX),D(NMAX),X,Y,DY,DIF,DIFT,HO,HP,W,DEN
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE W_TO_MSSO(W,MSSO)
      SAVE
C**********************************************************************
C
C   Subroutine to compute the real 4x4 Modified Stokes Scattering
C   Operator 
C   (see page 30 of Radar Polarimetry for Geoscience Applications) 
C   from a given complex 4x4 W Matrix.  Both matrices 
C   are assumed to be in the same convention (either BSA or FSA).
C
C   CALL PARAMETERS:
C
C   W     - Complex 4x4 W Matrix
C   MSSO  - Real 4x4 Modified Stokes Scattering Operator
C
C**********************************************************************
C
      INTEGER I,J
      COMPLEX W(4,4), VTI(4,4), VI(4,4), TMP(4,4), TMP1(4,4)
      REAL*4 MSSO(4,4)
      DATA ((VTI(I,J),J=1,4),I=1,4)/
     &            (1.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(1.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(0.0,0.0),(0.5,0.0),(0.5,0.0),
     &            (0.0,0.0),(0.0,0.0),(0.0,0.5),(0.0,-0.5)/
      DATA ((VI(I,J),J=1,4),I=1,4)/
     &            (1.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(1.0,0.0),(0.0,0.0),(0.0,0.0),
     &            (0.0,0.0),(0.0,0.0),(0.5,0.0),(0.0,0.5),
     &            (0.0,0.0),(0.0,0.0),(0.5,0.0),(0.0,-0.5)/

C
C   Zero out the Stokes Scattering Operator array.
C
      DO I=1,4
        DO J=1,4
          MSSO(I,J)=0.0
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C   Compute SSO Matrix
C----------------------------------------------------------------------

      CALL CCMULT(4,4,4,VTI,W,TMP)
      CALL CCMULT(4,4,4,TMP,VI,TMP1)
C
      DO I=1,4
        DO J=1,4
          MSSO(I,J)=REAL(TMP1(I,J))
        ENDDO
      ENDDO
C
      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE MIDEXP(FUNK,AA,BB,S,N)
      SAVE
C**********************************************************************
      INTEGER N,IT,J
      REAL*4 TNM,A,B,AA,BB,X,FUNC,SUM,DEL,DDEL,S
      EXTERNAL FUNK
      B=EXP(-AA)
      A=0.0
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(FUNK,0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(FUNK,X)
          X=X+DDEL
          SUM=SUM+FUNC(FUNK,X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
C**********************************************************************
C**********************************************************************
      FUNCTION FUNC(FUNK,X)
      SAVE
C**********************************************************************
      EXTERNAL FUNK
      REAL*4 FUNK,X,FUNC
      FUNC=FUNK(-ALOG(X))/X
      RETURN
      END
C**********************************************************************



