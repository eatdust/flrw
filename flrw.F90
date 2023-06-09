module flrw
!basic FLRW background integration with possibility to include g*(T)
  use flvars, only : cp, cep, flbgparam
  use flvars, only : c, Mpc,My,LPl
  use flvars, only : tolfl, big
  use flvars, only : flParams
  use flvars, only : correction_rdof, entropy_correction_rdof, equation_of_state_rdof

#ifndef SUNDIALS
  use funcutils, only : easydverk
#else
  use sundials, only : easycvode
#endif  
  
    
  implicit none

  private
       
  public cp, cep, LPl, Mpc, flbgparam
  
  public display_flparams, get_flrw_params
  public set_cosmo_params, set_fiducial_flparams

  public cosmic_time, comoving_distance

  public conformal_time, cosmic_time_normalized
  public conformal_time_normalized

  public comoving_distance_normalized, comoving_distance_today_normalized

  public redshift
  public redshift_equality  
  public redshift_normalized
  public redshift_chioa_normalized

  public redshift_conformal
  public redshift_conformal_normalized
  
  public hubble_today_hz, hubble_today, hubble_normalized, hubble_today_Lpl
  public omegarad_today, omegamat_today, omegalambda_today

  public conformal_deriv_second_a_over_a_normalized
  public hubble_normalized_scalefactor_square
  
#ifdef THERMAL
  public correction_rdof, entropy_correction_rdof
  public nothermal_hubble_normalized_scalefactor_square
#endif

  real(cp), save :: statbuffer
!$omp threadprivate(statbuffer)


contains

    
  subroutine set_fiducial_flparams()
    implicit none
!Planck18 TT+TE+EE+lowE+lensing with r     
     real(cp), parameter :: h = 0.6740
     real(cp), parameter :: OmegaBh2 = 0.02237
     real(cp), parameter :: OmegaCh2 = 0.1199
     real(cp), parameter :: OmegaPhoton = 5.38e-5
!all neutrino relativist -> 1.68
     real(cp), parameter :: OmegaRh2 = 1.68 * OmegaPhoton *h*h
     real(cp), parameter :: OmegaK = 0

     call set_cosmo_params(h,OmegaBh2,OmegaCh2,OmegaRh2,OmegaK)

   end subroutine set_fiducial_flparams



  subroutine set_cosmo_params(h,OmegaBh2,OmegaCh2,OmegaRh2,OmegaK)
#ifdef THERMAL
    use rdof, only : set_splines, energy_rdof_z, entropy_rdof_z
#endif    
    implicit none
    real(cp), intent(in) :: h, OmegaBh2,OmegaCh2,OmegaRh2,OmegaK
         
    flParams%h = h
    flParams%hubbleToday = 100._cp*h/c !Mpc^-1
!    flParams%hubbleToday = 100._cp*h/Mpc*My !My^-1
    flParams%OmegaM = (OmegaBh2 + OmegaCh2)/(h)**2
    flParams%OmegaR = OmegaRh2/h**2
    flParams%OmegaL = 1._cp - flParams%OmegaM - flParams%OmegaR - OmegaK
    flParams%OmegaK = OmegaK
    flParams%zeq = flParams%OmegaM/flParams%OmegaR - 1._cp

#ifdef THERMAL
    call set_splines()      
    flParams%go = energy_rdof_z(z=0._cp)
    flParams%qo = entropy_rdof_z(z=0._cp)
#else
    flParams%go = 1._cp
    flParams%qo = 1._cp
#endif    
    flParams%Hoto = cosmic_time_normalized(0._cp)
    flParams%Hochi = comoving_distance_normalized(big)
    
  end subroutine set_cosmo_params


  function get_flrw_params()
    implicit none
    type(flbgparam) :: get_flrw_params

    get_flrw_params = flParams

  end function get_flrw_params


  

  function hubble_today()
    implicit none
    real(cp) :: hubble_today
    hubble_today = flParams%hubbleToday
    
  end function hubble_today


  

  function hubble_today_hz()
    implicit none
    real(cp) :: hubble_today_hz

    hubble_today_hz = flParams%h * 100._cp/Mpc

  end function hubble_today_hz

!in reduced Planck length  
  function hubble_today_Lpl()
    implicit none
    real(cp) :: hubble_today_Lpl

    hubble_today_Lpl = flParams%h * 100._cp/c * Mpc / LPl

  end function hubble_today_Lpl

  
  
  function omegarad_today()
    implicit none
    real(cp) :: omegarad_today
    
    omegarad_today = flParams%OmegaR

  end function omegarad_today



  function omegamat_today()
    implicit none
    real(cp) :: omegamat_today
    
    omegamat_today = flParams%OmegaM

  end function omegamat_today
  

  function omegalambda_today()
    implicit none
    real(cp) :: omegalambda_today
    
    omegalambda_today = flParams%OmegaL

  end function omegalambda_today

  
  function cosmic_time_today_normalized()
    implicit none
    real(cp) :: cosmic_time_today_normalized

    cosmic_time_today_normalized = flParams%Hoto
    
  end function cosmic_time_today_normalized
    


  

 function comoving_distance_today_normalized()
    implicit none
    real(cp) :: comoving_distance_today_normalized

    comoving_distance_today_normalized = flParams%Hochi

  end function comoving_distance_today_normalized

  
  subroutine display_flparams()
    implicit none
    
    write(*,*)'display_flparams: '
    write(*,*)'H0      = ',flParams%hubbleToday
    write(*,*)'OmegaMo = ',flParams%OmegaM
    write(*,*)'OmegaRo = ',flParams%OmegaR
    write(*,*)'OmegaKo = ',flParams%OmegaK
    write(*,*)'OmegaLo = ',flParams%OmegaL
    write(*,*)'HoChi   = ',flParams%Hochi
    write(*,*)'Hoto    = ',flParams%Hoto
    write(*,*)'zeq     = ',flParams%zeq
    write(*,*)'go      = ',flParams%go
    write(*,*)'qo      = ',flParams%qo
    write(*,*)

  end subroutine display_flparams



  function redshift_equality()
    implicit none
    real(cp) :: redshift_equality

    redshift_equality = flParams%zeq

  end function redshift_equality



  

  function comoving_distance(redshift)   
    implicit none
    real(cp) :: comoving_distance
    real(cp), intent(in) :: redshift

    comoving_distance = comoving_distance_normalized(redshift) &
         /flParams%hubbleToday

  end function comoving_distance



  
  function comoving_distance_normalized(redshift)
    implicit none
    real(cp) :: comoving_distance_normalized
    real(cp), intent(in) :: redshift

    integer, parameter :: neq = 1

    real(cp) :: scaleFactor, scaleStop
    real(cp), dimension(neq) :: chi

    real(cp), parameter :: tol = tolfl
    real(cp), parameter :: ztol = 100 * tol
    
    if (redshift.lt.ztol) then
       comoving_distance_normalized = redshift
       return
    endif
    
!that's the scale factor in unit of a_now    
    scaleFactor = 1._cp/(1._cp + redshift)
    scaleStop = 1._cp
    chi = 0._cp

#ifndef SUNDIALS    
    call easydverk(neq,deriv_conftime_normalized,scaleFactor,chi,scaleStop,tol)
#else
    chi = easycvode(neq,deriv_conftime_normalized,scaleFactor,chi,scaleStop,tol)
#endif
    
    comoving_distance_normalized = chi(neq)

  end function comoving_distance_normalized
  


 

  function conformal_time(redshift)
    implicit none
    real(cp) :: conformal_time
    real(cp), intent(in) :: redshift
   
    conformal_time = conformal_time_normalized(redshift)/flParams%hubbleToday

  end function conformal_time




  function conformal_time_normalized(redshift)
    implicit none
    real(cp) :: conformal_time_normalized
    real(cp), intent(in) :: redshift

    real(cp), parameter :: tol = tolfl
    
    integer, parameter :: neq=1

    real(cp) :: scaleFactor, scaleStart
    real(cp), dimension(neq) :: eta

    !that's the scale factor in unit of a_now    
    scaleFactor = 1._cp/(1._cp + redshift)
    scaleStart = 0._cp
    eta = 0._cp

#ifndef SUNDIALS    
    call easydverk(neq,deriv_conftime_normalized,scaleStart,eta,scaleFactor,tol)
#else
    eta = easycvode(neq,deriv_conftime_normalized,scaleStart,eta,scaleFactor,tol)
#endif    
    !this a0 x eta x Ho
    conformal_time_normalized = eta(neq)

  end function conformal_time_normalized
  

#ifndef SUNDIALS

  subroutine deriv_conftime_normalized(neq,aOveraZero,eta,etaPrime)
    implicit none          
    integer :: neq
    real(cp) :: aOveraZero
    real(cp), dimension(neq) :: eta, etaPrime

    etaPrime = 1._cp/hubble_normalized_scalefactor_square(aOveraZero)
    
  end subroutine deriv_conftime_normalized

#else
  
  function deriv_conftime_normalized(neq,aOveraZero,eta) result(etaPrime)
    implicit none
    integer, intent(in) :: neq
    real(cp), intent(in) :: aOveraZero
    real(cp), dimension(neq), intent(in) :: eta
    real(cp), dimension(neq) :: etaPrime
    
    etaPrime = 1._cp/hubble_normalized_scalefactor_square(aOveraZero)
    
  end function deriv_conftime_normalized

#endif
  

  function cosmic_time(redshift)
    implicit none
    real(cp) :: cosmic_time
    real(cp), intent(in) :: redshift
  
    cosmic_time = cosmic_time_normalized(redshift)/flParams%hubbleToday

  end function cosmic_time


  
  function cosmic_time_normalized(redshift)
    implicit none
    real(cp) :: cosmic_time_normalized
    real(cp), intent(in) :: redshift
    integer, parameter :: neq=1
       
    real(cp) :: scaleFactor, scaleStart
    real(cp), dimension(neq) :: t
   
!that's the scale factor in unit of a_now    
    scaleFactor = 1._cp/(1._cp + redshift)
    scaleStart = 0._cp
    t = 0._cp

#ifndef SUNDIALS    
    call easydverk(neq,deriv_cosmic_time_normalized,scaleStart,t,scaleFactor)
#else
    t = easycvode(neq,deriv_cosmic_time_normalized,scaleStart,t,scaleFactor)
#endif
    
    cosmic_time_normalized = t(neq)

  end function cosmic_time_normalized

  
#ifndef SUNDIALS
  
  subroutine deriv_cosmic_time_normalized(neq,aOveraZero,t,tPrime)
    implicit none          
    integer :: neq
    real(cp) :: aOveraZero
    real(cp), dimension(neq) :: t,tPrime

    tPrime = aOveraZero/hubble_normalized_scalefactor_square(aOveraZero)

  end subroutine deriv_cosmic_time_normalized

#else

  function deriv_cosmic_time_normalized(neq,aOveraZero,t) result(tPrime)
    implicit none
    integer, intent(in) :: neq
    real(cp), intent(in) :: aOveraZero
    real(cp), dimension(neq), intent(in) :: t
    real(cp), dimension(neq) :: tPrime
    
    tPrime = aOveraZero/hubble_normalized_scalefactor_square(aOveraZero)
    
  end function deriv_cosmic_time_normalized

#endif  

  
  function redshift_normalized(cosmicTimeHo)
    use flapprox, only : cosmic_radtime_normalized
    use flapprox, only : redshift_radtime_normalized
    implicit none
    real(cp) :: redshift_normalized
    real(cp), intent(in) :: cosmicTimeHo

    integer, parameter :: neq=1

!should be in the high energy regime for thermal history 2d17
    real(cp), parameter :: zDeepRad = 2d17
!for z>zbig, we invert analytic formulae    
    real(cp), parameter :: zBig = 1d7

    real(cp) :: tsmallHo
    real(cp) :: tnowHo,sqrtomR
    real(cp) :: scaleFactor,tHoDeepRad
    real(cp), dimension(neq) :: a

    logical, parameter :: useHighZApprox = .true.

    if (useHighZApprox) then
    
       tHodeepRad = cosmic_radtime_normalized(zDeepRad)

!rdof is constant in the SM for z>zdeeprad, we do not need to
!integrate anything nor invert thermal history formulae
       if (cosmicTimeHo.lt.tHoDeepRad) then
#ifndef THERMAL    
          sqrtomR = sqrt(flParams%OmegaR)
#else
          sqrtomR = sqrt(flParams%OmegaR * correction_rdof(zDeepRad))
#endif

          redshift_normalized = 1._cp/sqrt(2._cp*cosmicTimeHo &
               *sqrtomR) - 1._cp
          
          return
          
       endif

!zbig is very large redshift for which numerical integration may be
!inaccurate, in the range zbig, zdeeprad, we can skip numerical
!integration but still need to invert thermal history algebraic
!equations
       tsmallHo = cosmic_radtime_normalized(zbig)

       if (cosmicTimeHo.lt.tsmallHo) then

          redshift_normalized = redshift_radtime_normalized(cosmicTimeHo,Q=1._cp)

          return
          
       endif

    endif
       
!Exact integration in between, or all the time is useHighZApprox is false

    tnowHo = cosmic_time_normalized(0._cp)
    a = 1._cp

#ifndef SUNDIALS    
    call easydverk(neq,deriv_scalefactor_normalized,tnowHo,a,cosmicTimeHo)
#else
    a = easycvode(neq,deriv_scalefactor_normalized,tnowHo,a,cosmicTimeHo)
#endif
    
    scaleFactor = a(neq)

    
    redshift_normalized = 1._cp/scaleFactor - 1._cp
    
  end function redshift_normalized


#ifndef SUNDIALS  
   subroutine deriv_scalefactor_normalized(neq,tHo,a,aPrime)
    implicit none
    integer :: neq
    real(cp) :: tHo
    real(cp), dimension(neq) :: a, aPrime

    aPrime = hubble_normalized_scalefactor_square(a(neq))/a(neq)
    
  end subroutine deriv_scalefactor_normalized

#else
 function deriv_scalefactor_normalized(neq,tHo,a) result(aPrime)
    implicit none
    integer, intent(in) :: neq
    real(cp), intent(in) :: tHo
    real(cp), dimension(neq), intent(in) :: a
    real(cp), dimension(neq) :: aPrime
    
    aPrime = hubble_normalized_scalefactor_square(a(neq))/a(neq)
    
  end function deriv_scalefactor_normalized
#endif  

  
  function redshift(cosmicTime)
    implicit none
    real(cp) :: redshift
    real(cp), intent(in) :: cosmicTime
   
    redshift = redshift_normalized(cosmicTime*flParams%hubbleToday)

  end function redshift




  function redshift_conformal_normalized(etaHo)
    use flapprox, only : conformal_radmattime_normalized
    use flapprox, only : redshift_conformal_radmattime_normalized
    implicit none
    real(cp) :: redshift_conformal_normalized
    real(cp), intent(in) :: etaHo
    
    integer, parameter :: neq=1
    real(cp), dimension(neq) :: a

!should be in the high energy regime for thermal history 2d17
    real(cp), parameter :: zDeepRad = 2d17
!for z>zbig, we invert analytic formulae    
    real(cp), parameter :: zBig = 1d7

    real(cp) :: etasmallHo, etaHoDeepRad
    real(cp) :: etanowHo, scaleFactor, sqrtomR

    real(cp), parameter :: tol = tolfl


    logical, parameter :: useHighZApprox = .true.

    
    if (useHighZApprox) then

       etaHoDeepRad = conformal_radmattime_normalized(zDeepRad)

!rdof is constant in the SM for z>zdeeprad, we do not need to
!integrate anything nor invert thermal history formulae
       if (etaHo.lt.etaHoDeepRad) then
#ifndef THERMAL    
          sqrtomR = sqrt(flParams%OmegaR)
#else
          sqrtomR = sqrt(flParams%OmegaR * correction_rdof(zDeepRad))
#endif

          redshift_conformal_normalized = 1._cp &
               /(0.25_cp*flParams%OmegaM * etaHo*etaHo + sqrtomR * etaHo) &
               - 1._cp
          
          return
          
       endif
       
       etasmallHo = conformal_radmattime_normalized(zBig)

       if (etaHo.lt.etasmallHo) then

          redshift_conformal_normalized = redshift_conformal_radmattime_normalized(etaHo,Q=1._cp)

          return
          
       endif
       
    endif

    
    
!Exact integration in between, or all the time is useHighZApprox is false
    etanowHo = conformal_time_normalized(0._cp)
    a = 1._cp

#ifndef SUNDIALS    
    call easydverk(neq,deriv_conformal_scalefactor_normalized,etanowHo,a,etaHo,tol)
#else
    a = easycvode(neq,deriv_conformal_scalefactor_normalized,etanowHo,a,etaHo,tol)
#endif    
    scaleFactor = a(neq)

    redshift_conformal_normalized = 1._cp/scaleFactor - 1._cp
        
  end function redshift_conformal_normalized

  
#ifndef SUNDIALS
  subroutine deriv_conformal_scalefactor_normalized(neq,etaHo,a,aPrime)
    implicit none
    integer :: neq
    real(cp) :: etaHo
    real(cp), dimension(neq) :: a, aPrime

    aPrime = hubble_normalized_scalefactor_square(a(neq))
    
  end subroutine deriv_conformal_scalefactor_normalized

#else

   function deriv_conformal_scalefactor_normalized(neq,etaHo,a) result(aPrime)
    implicit none
    integer, intent(in) :: neq
    real(cp), intent(in) :: etaHo
    real(cp), dimension(neq), intent(in) :: a
    real(cp), dimension(neq) :: aPrime
    
    aPrime = hubble_normalized_scalefactor_square(a(neq))
    
  end function deriv_conformal_scalefactor_normalized

#endif

  
  function redshift_conformal(conformalTime)
    implicit none
    real(cp) :: redshift_conformal
    real(cp), intent(in) :: conformalTime

    redshift_conformal = redshift_conformal_normalized(conformalTime*flParams%hubbleToday)

  end function redshift_conformal
  

  
!a^2 H
  function hubble_scalefactor_square(scaleFactor)
    implicit none
    real(cp) :: hubble_scalefactor_square
    real(cp), intent(in) :: scaleFactor

    hubble_scalefactor_square = flParams%hubbleToday &
         * hubble_normalized_scalefactor_square(scaleFactor)
  
  end function hubble_scalefactor_square




  
   
!a^2 H/Ho  
#ifndef THERMAL

  function hubble_normalized_scalefactor_square(scaleFactor)
    implicit none
    real(cp) :: hubble_normalized_scalefactor_square
    real(cp), intent(in) :: scaleFactor

    hubble_normalized_scalefactor_square = nothermal_hubble_normalized_scalefactor_square(scaleFactor)
    
  end function hubble_normalized_scalefactor_square


  
  function conformal_deriv_second_a_over_a_normalized(scaleFactor)
    implicit none
    real(cp) :: conformal_deriv_second_a_over_a_normalized
    real(cp), intent(in) :: scaleFactor

    conformal_deriv_second_a_over_a_normalized &
         = nothermal_hubble_normalized_scalefactor_square(scaleFactor)
    
  end function conformal_deriv_second_a_over_a_normalized
    
#else  
    
  function hubble_normalized_scalefactor_square(scaleFactor)
    implicit none
    real(cp) :: hubble_normalized_scalefactor_square
    real(cp), intent(in) :: scaleFactor

    real(cp) :: z

    z = 1._cp/scaleFactor - 1._cp
    
    hubble_normalized_scalefactor_square &
         = sqrt(flParams%OmegaR*correction_rdof(z) &
         + flParams%OmegaM*scaleFactor &
         + flParams%OmegaL*scaleFactor**4 &
         + flParams%OmegaK* scaleFactor**2)
    
  end function hubble_normalized_scalefactor_square


  
!this is (a''/a) * 1/(ao^2 Ho^2) = (calH^2 + calH')/(calHo^2)
  function conformal_deriv_second_a_over_a_normalized(scaleFactor)
    implicit none
    real(cp) :: conformal_deriv_second_a_over_a_normalized
    real(cp), intent(in) :: scaleFactor

    real(cp) :: z

    z = 1._cp/scaleFactor - 1._cp
    
    conformal_deriv_second_a_over_a_normalized = flParams%OmegaK &
         + 0.5_cp * flParams%OmegaM/scaleFactor &
         + 2._cp*flParams%OmegaL*scaleFactor*scaleFactor &
         + 0.5_cp*flParams%OmegaR*(1._cp - 3._cp*equation_of_state_rdof(z)) &
         /scaleFactor**2

  end function conformal_deriv_second_a_over_a_normalized

  
#endif
  
  
  function nothermal_hubble_normalized_scalefactor_square(scaleFactor)
    implicit none
    real(cp) :: nothermal_hubble_normalized_scalefactor_square
    real(cp), intent(in) :: scaleFactor

    nothermal_hubble_normalized_scalefactor_square &
         = sqrt(flParams%OmegaR + flParams%OmegaM*scaleFactor &
         + flParams%OmegaL*scaleFactor**4 &
         + flParams%OmegaK* scaleFactor**2)

  end function nothermal_hubble_normalized_scalefactor_square

 

  
  function hubble_normalized(z)
    implicit none
    real(cp) :: hubble_normalized
    real(cp), intent(in) :: z
    real(cp) :: scaleFactor, a2HoverHo

    scaleFactor = 1._cp/(1._cp + z)

    a2HoverHo = hubble_normalized_scalefactor_square(scaleFactor)
    
    hubble_normalized = a2HoverHo/scaleFactor/ScaleFactor
    
  end function hubble_normalized


!this is (a''/a) * 1/(ao^2 Ho^2) = (calH^2 + calH')/(calHo^2) with wrad = 1/3
  function nothermal_conformal_deriv_second_a_over_a_normalized(z)
    implicit none
    real(cp) :: nothermal_conformal_deriv_second_a_over_a_normalized
    real(cp), intent(in) :: z

    real(cp) :: scaleFactor

    scaleFactor = 1._cp/(1._cp + z)

    
    nothermal_conformal_deriv_second_a_over_a_normalized = flParams%OmegaK &
         + 0.5_cp * flParams%OmegaM/scaleFactor &
         + 2._cp*flParams%OmegaL*scaleFactor*scaleFactor        

  end function nothermal_conformal_deriv_second_a_over_a_normalized
  
  
  
  !returns z such that Ho Chi(z) (1+z) = Hochioa
  function redshift_chioa_normalized(Hochioa)
    use funcutils, only : brent
    implicit none
    real(cp) :: redshift_chioa_normalized
    real(cp), intent(in) :: Hochioa

    real(cp) :: lnzp1min, lnzp1max
    real(cp), parameter :: tol = tolfl
    
    lnzp1min = 0._cp
    lnzp1max = log(big)/4._cp
    statbuffer = Hochioa
    
    redshift_chioa_normalized &
         = exp(brent(lnzp1min,lnzp1max,find_redshift_chioa_normalized,tol)) - 1._cp
         
  end function redshift_chioa_normalized

  function find_redshift_chioa_normalized(lnzp1)
    implicit none
    real(cp) :: find_redshift_chioa_normalized
    real(cp), intent(in) :: lnzp1

    real(cp) :: zp1,zp1Hochi
    
    zp1 = exp(lnzp1)
    zp1Hochi = zp1*comoving_distance_normalized(zp1-1._cp)
    
    find_redshift_chioa_normalized = zp1Hochi - statbuffer
    
  end function find_redshift_chioa_normalized

  

  
end module flrw
