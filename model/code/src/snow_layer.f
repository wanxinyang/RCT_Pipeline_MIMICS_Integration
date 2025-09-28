c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c   Subroutine to compute the transmissivity for a layer of
c   snow. 
c***********************************************************************
c
        subroutine  snow_layer(eps_snow,theta_loc,two_way_tau)
        save
c
c***********************************************************************
c
        real two_way_tau, kappa_s, theta_loc
        complex eps_snow

        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL T_SNOW
c
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /R_SNOW/ T_SNOW
c
c***********************************************************************
c
c   slab model -- snow dielectric is given
c              -- neglect scattering in the snow layer
c
c   compute two-way transmissivity of the snow layer
c
        kappa_s = 2.0*k0*abs(aimag(csqrt(eps_snow)))
        two_way_tau = exp(-2.0*kappa_s*t_snow/cos(theta_loc))
c
        return
        end
c
c***********************************************************************

