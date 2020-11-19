module flrw
!basic FLRW background integration with possibility to include g*(T)
  implicit none

  integer, parameter :: cp = kind(1._8)
  integer, parameter :: cep = kind(1._16)
  
  private


  interface cosmic_scalingtime
     module procedure cp_cosmic_scalingtime, ep_cosmic_scalingtime
  end interface

  interface cosmic_scalingtime_normalized
     module procedure cp_cosmic_scalingtime_normalized, ep_cosmic_scalingtime_normalized
  end interface cosmic_scalingtime_normalized
  
  interface cosmic_mattime
     module procedure cp_cosmic_mattime, ep_cosmic_mattime
  end interface cosmic_mattime

  interface cosmic_mattime_normalized
     module procedure cp_cosmic_mattime_normalized, ep_cosmic_mattime_normalized
  end interface cosmic_mattime_normalized

  interface cosmic_radtime
     module procedure cp_cosmic_radtime, ep_cosmic_radtime
  end interface cosmic_radtime

  interface cosmic_radtime_normalized
     module procedure cp_cosmic_radtime_normalized, ep_cosmic_radtime_normalized
  end interface cosmic_radtime_normalized

  interface redshift_scalingtime
     module procedure cp_redshift_scalingtime, ep_redshift_scalingtime
  end interface redshift_scalingtime

  interface redshift_scalingtime_normalized
     module procedure cp_redshift_scalingtime_normalized, ep_redshift_scalingtime_normalized
  end interface redshift_scalingtime_normalized

  interface redshift_radtime_normalized
     module procedure cp_redshift_radtime_normalized, ep_redshift_radtime_normalized
  end interface redshift_radtime_normalized
  
  interface redshift_toa_scalingtime_normalized
     module procedure cp_redshift_toa_scalingtime_normalized, ep_redshift_toa_scalingtime_normalized
  end interface redshift_toa_scalingtime_normalized
  
  interface redshift_toa_radtime_normalized
     module procedure cp_redshift_toa_radtime_normalized
  end interface redshift_toa_radtime_normalized

  
  type flbgparam
     real(cp) :: h, hubbleToday !Mpc
     real(cp) :: OmegaM, OmegaR, OmegaL
 !for convenience
     real(cp) :: OmegaK
     real(cp) :: zeq
     real(cp) :: Hoto, Hochi
!for thermal history
     real(cp) :: go, qo
     
  end type flbgparam

  real(cp), parameter :: c = 299792.458 !km/s
  real(cp), parameter :: Mpc = 3.0857d19 !km 
  real(cp), parameter :: My = 365.25*24*3600*1d6 !s
  real(cp), parameter :: LPl= 8.1026d-38 !km
  
  type(flbgparam), save :: flParams
 
  real(cp), save :: statbuffer, statpower
!$omp threadprivate(statbuffer,statpower)  

  real(cp), dimension(2), save :: partbuffer
!$omp threadprivate(partbuffer)  
  
  real(cp), parameter :: tolfl = 100*epsilon(1._cp)
  real(cp), parameter :: big = epsilon(1._cp)*huge(1._cp)



  public cp, cep, LPl, Mpc
  
  public display_flparams
  public set_cosmo_params, set_fiducial_flparams
  public cosmic_time, cosmic_scalingtime
  public conformal_time, redshift, redshift_scalingtime, redshift_equality
  public redshift_crossing, comoving_distance

  public cosmic_scalingtime_normalized, redshift_scalingtime_normalized
  public cosmic_mattime_normalized, cosmic_radtime_normalized
  public redshift_toa_scalingtime_normalized, comoving_distance_normalized
  public redshift_radtime_normalized, redshift_toa_radtime_normalized
  public conformal_time_normalized, cosmic_time_normalized, redshift_normalized
  public redshift_chioa_normalized, cosmic_scalingtime_today_normalized
  public cosmic_time_today_normalized, comoving_distance_today_normalized
  public redshift_tchipower_scalingtime_normalized
  
  public hubble_today_hz, hubble_today, hubble_normalized, hubble_today_Lpl
  public omegarad_today, omegamat_today, omegalambda_today
  
#ifdef THERMAL
  public correction_rdof, entropy_correction_rdof
  public nothermal_hubble_normalized_scalefactor_square
#endif

  
contains

  subroutine set_fiducial_flparams()
    implicit none
     
     real(cp), parameter :: h = 0.679
     real(cp), parameter :: OmegaBh2 = 0.02227
     real(cp), parameter :: OmegaCh2 = 0.1184
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
    


  function cosmic_scalingtime_today_normalized()
    implicit none
    real(cp) :: cosmic_scalingtime_today_normalized
    
    cosmic_scalingtime_today_normalized = cosmic_scalingtime_normalized(0._cp)

  end function cosmic_scalingtime_today_normalized



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


  
  function redshift_crossing()
    implicit none
    real(cp) :: redshift_crossing

    redshift_crossing = 9._cp*(1._cp + flParams%zeq)/16._cp - 1._cp

  end function redshift_crossing


  

  function comoving_distance(redshift)   
    implicit none
    real(cp) :: comoving_distance
    real(cp), intent(in) :: redshift

    comoving_distance = comoving_distance_normalized(redshift) &
         /flParams%hubbleToday

  end function comoving_distance



  
  function comoving_distance_normalized(redshift)
    use functools, only : easydverk    
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
    
    call easydverk(neq,deriv_conftime_normalized,scaleFactor,chi,scaleStop,tol)

    comoving_distance_normalized = chi(neq)

  end function comoving_distance_normalized
  


 

  function conformal_time(redshift)
    implicit none
    real(cp) :: conformal_time
    real(cp), intent(in) :: redshift
   
    conformal_time = conformal_time_normalized(redshift)/flParams%hubbleToday

  end function conformal_time




  function conformal_time_normalized(redshift)
    use functools, only : easydverk
    implicit none
    real(cp) :: conformal_time_normalized
    real(cp), intent(in) :: redshift

    integer, parameter :: neq=1

    real(cp) :: scaleFactor, scaleStart
    real(cp), dimension(neq) :: eta

    !that's the scale factor in unit of a_now    
    scaleFactor = 1._cp/(1._cp + redshift)
    scaleStart = 0._cp
    eta = 0._cp

    call easydverk(neq,deriv_conftime_normalized,scaleStart,eta,scaleFactor)

    !this a0 x eta x Ho
    conformal_time_normalized = eta(neq)

  end function conformal_time_normalized
  



  subroutine deriv_conftime_normalized(neq,aOveraZero,eta,etaPrime)
    implicit none          
    integer :: neq
    real(cp) :: aOveraZero
    real(cp), dimension(neq) :: eta, etaPrime

    etaPrime = 1._cp/hubble_normalized_scalefactor_square(aOveraZero)

  end subroutine deriv_conftime_normalized

  

  function cosmic_time(redshift)
    implicit none
    real(cp) :: cosmic_time
    real(cp), intent(in) :: redshift
  
    cosmic_time = cosmic_time_normalized(redshift)/flParams%hubbleToday

  end function cosmic_time


  
  function cosmic_time_normalized(redshift)
    use functools, only : easydverk
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

    call easydverk(neq,deriv_cosmic_time_normalized,scaleStart,t,scaleFactor)
    
    cosmic_time_normalized = t(neq)

  end function cosmic_time_normalized

  

  subroutine deriv_cosmic_time_normalized(neq,aOveraZero,t,tPrime)
    implicit none          
    integer :: neq
    real(cp) :: aOveraZero
    real(cp), dimension(neq) :: t,tPrime

!    real(cp) :: scaleFactor
!    scaleFactor = aOveraZero

    tPrime = aOveraZero/hubble_normalized_scalefactor_square(aOveraZero)

  end subroutine deriv_cosmic_time_normalized


  
  function redshift_normalized(cosmicTimeHo)
    use functools, only : easydverk
    implicit none
    real(cp) :: redshift_normalized
    real(cp), intent(in) :: cosmicTimeHo

    integer, parameter :: neq=1
!should be in the high energy regime for thermal history 2d17
    real(cp), parameter :: zDeepRad = 2d17
    real(cp), parameter :: aDeepRad = 1._cp/(1._cp + zDeepRad)

!for z>zbig, we invert analytic formulae    
    real(cp), parameter :: zBig = 1d7
    real(cp) :: tsmallHo
    
    real(cp) :: tnowHo,sqrtomR
    real(cp) :: scaleFactor,tHoDeepRad
    real(cp), dimension(neq) :: a,lna


    
    
    tHodeepRad = cosmic_radtime_normalized(zDeepRad)

!rdof is constant in the SM for z>zdeeprad    
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

!that is very large redshift for which numerical integration may be inaccurate    
    tsmallHo = cosmic_radtime_normalized(zbig)
   
    if (cosmicTimeHo.lt.tsmallHo) then
       redshift_normalized = redshift_radtime_normalized(cosmicTimeHo,Q=1._cp)
       return
    endif
    
    
!let's integrate    

    tnowHo = cosmic_time_normalized(0._cp)
    a = 1._cp

    call easydverk(neq,deriv_scalefactor_normalized,tnowHo,a,cosmicTimeHo)
    scaleFactor = a(neq)

    
    redshift_normalized = 1._cp/scaleFactor - 1._cp
    
  end function redshift_normalized


  
   subroutine deriv_scalefactor_normalized(neq,tHo,a,aPrime)
    implicit none
    integer :: neq
    real(cp) :: tHo
    real(cp), dimension(neq) :: a, aPrime

    aPrime = hubble_normalized_scalefactor_square(a(neq))/a(neq)
    
  end subroutine deriv_scalefactor_normalized

  

  
  
  function redshift(cosmicTime)
    implicit none
    real(cp) :: redshift
    real(cp), intent(in) :: cosmicTime
   
    redshift = redshift_normalized(cosmicTime*flParams%hubbleToday)

  end function redshift



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

    hubble_normalized_scalefactor_square &
         = sqrt(flParams%OmegaR + flParams%OmegaM*scaleFactor &
         + flParams%OmegaL*scaleFactor**4 &
         + flParams%OmegaK* scaleFactor**2)
    
  end function hubble_normalized_scalefactor_square

!to potentially implement a step like minimal thing
  function correction_rdof(z)
    implicit none
    real(cp), intent(in) :: z
    real(cp) :: correction_rdof
    
    correction_rdof = 1._cp
    
  end function correction_rdof

  
  function entropy_correction_rdof(z)
    real(cp), intent(in) :: z
    real(cp) :: entropy_correction_rdof


    entropy_correction_rdof = 1._cp
    
  end function entropy_correction_rdof

  
#else  


  function nothermal_hubble_normalized_scalefactor_square(scaleFactor)
    implicit none
    real(cp) :: nothermal_hubble_normalized_scalefactor_square
    real(cp), intent(in) :: scaleFactor

    nothermal_hubble_normalized_scalefactor_square &
         = sqrt(flParams%OmegaR + flParams%OmegaM*scaleFactor &
         + flParams%OmegaL*scaleFactor**4 &
         + flParams%OmegaK* scaleFactor**2)

  end function nothermal_hubble_normalized_scalefactor_square

 

  
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




!one spline evaluation
  function correction_rdof(z)
    use rdof, only : dp, correction_rdof_z
    implicit none
    real(cp), intent(in) :: z
    real(cp) :: correction_rdof
    
    correction_rdof = correction_rdof_z(real(z,dp))
    
  end function correction_rdof
  
  
  
!requires two spline evaluations
  function old_correction_rdof(z)
    use rdof, only : dp, energy_rdof_z, entropy_rdof_z
    implicit none
    real(cp), intent(in) :: z
    real(cp) :: old_correction_rdof
    real(cp), parameter :: fourthird = 4._cp/3._cp
    real(cp) :: g,q,go,qo
    real(dp) :: redshift

    redshift = real(z,dp)
    
    g = real(energy_rdof_z(redshift),cp)
    q = real(entropy_rdof_z(redshift),cp)

    go = real(energy_rdof_z(0._dp),cp)
    qo = real(entropy_rdof_z(0._dp),cp)

    old_correction_rdof = g/go * (qo/q)**fourthird
    
  end function old_correction_rdof


  function entropy_correction_rdof(z)
    use rdof, only : dp, entropy_rdof_z
    real(cp) :: entropy_correction_rdof
    real(cp), intent(in) :: z
    real(dp) :: redshift

    redshift = real(z,dp)

    entropy_correction_rdof = (real(entropy_rdof_z(redshift),cp)/flParams%qo)
    
  end function entropy_correction_rdof

  
#endif
  

   

  
  

  function hubble_normalized(z)
    implicit none
    real(cp) :: hubble_normalized
    real(cp), intent(in) :: z
    real(cp) :: scaleFactor, a2HoverHo

    scaleFactor = 1._cp/(1._cp + z)

    a2HoverHo = hubble_normalized_scalefactor_square(scaleFactor)
    
    hubble_normalized = a2HoverHo/scaleFactor/ScaleFactor
    
  end function hubble_normalized
  
  

  
  function cp_cosmic_scalingtime(redshift)
    implicit none
    real(cp) :: cp_cosmic_scalingtime
    real(cp), intent(in) :: redshift
    
    cp_cosmic_scalingtime = cp_cosmic_scalingtime_normalized(redshift)/flParams%hubbleToday

  end function cp_cosmic_scalingtime


  
  function ep_cosmic_scalingtime(redshift)
    implicit none
    real(cep) :: ep_cosmic_scalingtime
    real(cep), intent(in) :: redshift
    
    ep_cosmic_scalingtime = ep_cosmic_scalingtime_normalized(redshift)/real(flParams%hubbleToday,cep)

  end function ep_cosmic_scalingtime


  
  function cp_cosmic_scalingtime_normalized(redshift)
    implicit none
    real(cp) :: cp_cosmic_scalingtime_normalized
    real(cp), intent(in) :: redshift

    real(cp) :: zcross
    
    zcross = redshift_crossing()

    if (redshift.gt.zcross) then
       cp_cosmic_scalingtime_normalized = cp_cosmic_radtime_normalized(redshift)
    elseif (redshift.le.zcross) then
       cp_cosmic_scalingtime_normalized = cp_cosmic_mattime_normalized(redshift)
    else
       write(*,*)'redshift= ',redshift
       stop 'cp_cosmic_scalingtime_normalized: internal error!'
    endif

  end function cp_cosmic_scalingtime_normalized
  


  
  function ep_cosmic_scalingtime_normalized(redshift)
    implicit none
    real(cep) :: ep_cosmic_scalingtime_normalized
    real(cep), intent(in) :: redshift

    real(cep) :: zcross
    
    zcross = real(redshift_crossing(),cep)

    if (redshift.gt.zcross) then
       ep_cosmic_scalingtime_normalized = ep_cosmic_radtime_normalized(redshift)
    elseif (redshift.le.zcross) then
       ep_cosmic_scalingtime_normalized = ep_cosmic_mattime_normalized(redshift)
    else
       write(*,*)'redshift= ',redshift
       stop 'ep_cosmic_scalingtime_normalized: internal error!'
    endif

  end function ep_cosmic_scalingtime_normalized
  

  

  
  function cp_cosmic_mattime(redshift)
    implicit none
!pure matter era
    real(cp) :: cp_cosmic_mattime
    real(cp), intent(in) :: redshift

    cp_cosmic_mattime = cp_cosmic_mattime_normalized(redshift)/flParams%hubbleToday
    
  end function cp_cosmic_mattime


  

  function ep_cosmic_mattime(redshift)
    implicit none
!pure matter era
    real(cep) :: ep_cosmic_mattime
    real(cep), intent(in) :: redshift

    ep_cosmic_mattime = ep_cosmic_mattime_normalized(redshift)/real(flParams%hubbleToday,cep)
    
  end function ep_cosmic_mattime

  
  
  
  function cp_cosmic_mattime_normalized(redshift)
    implicit none
!pure matter era
    real(cp) :: cp_cosmic_mattime_normalized
    real(cp), intent(in) :: redshift

    real(cp) :: sqrtomM

    sqrtomM = sqrt(flParams%OmegaM)
    
    cp_cosmic_mattime_normalized = 2._cp/(3._cp &
         * sqrtomM)/(1._cp + redshift)**(3._cp/2._cp)
    
  end function cp_cosmic_mattime_normalized



  function ep_cosmic_mattime_normalized(redshift)
    implicit none
!pure matter era
    real(cep) :: ep_cosmic_mattime_normalized
    real(cep), intent(in) :: redshift

    real(cep) :: sqrtomM
    
    sqrtomM = sqrt(real(flParams%OmegaM,cep))
    
    ep_cosmic_mattime_normalized = 2._cep/(3._cep &
         * sqrtomM)/(1._cep + redshift)**(3._cep/2._cep)
    
  end function ep_cosmic_mattime_normalized
  

  
  function cp_cosmic_radtime(redshift)
    implicit none
!pure radiation era or THERMAL leading order
    real(cp) :: cp_cosmic_radtime
    real(cp), intent(in) :: redshift
    
    cp_cosmic_radtime = cp_cosmic_radtime_normalized(redshift)/flParams%hubbleToday
       
  end function cp_cosmic_radtime



  
  function ep_cosmic_radtime(redshift)
    implicit none
!pure radiation era or THERMAL leading order
    real(cep) :: ep_cosmic_radtime
    real(cep), intent(in) :: redshift
    
    ep_cosmic_radtime = ep_cosmic_radtime_normalized(redshift)/real(flParams%hubbleToday,cep)
       
  end function ep_cosmic_radtime
  
  
  
  function cp_cosmic_radtime_normalized(redshift)
    implicit none
    !pure radiation era
    real(cp) :: cp_cosmic_radtime_normalized
    real(cp), intent(in) :: redshift

    real(cp) :: sqrtomR

    
#ifndef THERMAL    
    sqrtomR = sqrt(flParams%OmegaR)
#else
!this is only the leading order term of an expansion in correction_rdof derivatives
    sqrtomR = sqrt(flParams%OmegaR * correction_rdof(redshift))
#endif
    
    cp_cosmic_radtime_normalized = 1._cp/(2._cp &
         * sqrtomR)/(1._cp + redshift)**2

  end function cp_cosmic_radtime_normalized



  
  function ep_cosmic_radtime_normalized(redshift)
    implicit none
!pure radiation era
    real(cep) :: ep_cosmic_radtime_normalized
    real(cep), intent(in) :: redshift

    real(cep) :: sqrtomR

#ifndef THERMAL
    sqrtomR = sqrt(real(flParams%OmegaR,cep))
#else
    sqrtomR = sqrt(real(flParams%OmegaR*correction_rdof(real(redshift,cp)),cep))
#endif
    
    ep_cosmic_radtime_normalized = 1._cep/(2._cep &
         * sqrtomR )/(1._cep + redshift)**2
       
  end function ep_cosmic_radtime_normalized




  
  function cp_redshift_scalingtime(cosmicTime)
    implicit none
    real(cp) :: cp_redshift_scalingtime
    real(cp), intent(in) :: cosmicTime

    cp_redshift_scalingtime = cp_redshift_scalingtime_normalized(cosmicTime*flParams%hubbleToday)
    
  end function cp_redshift_scalingtime




  function ep_redshift_scalingtime(cosmicTime)
    implicit none
    real(cep) :: ep_redshift_scalingtime
    real(cep), intent(in) :: cosmicTime

    ep_redshift_scalingtime = ep_redshift_scalingtime_normalized(cosmicTime*real(flParams%hubbleToday,cep))
    
  end function ep_redshift_scalingtime
  


  
  
  function cp_redshift_scalingtime_normalized(cosmicTimeHo)
    implicit none
    real(cp) :: cp_redshift_scalingtime_normalized
    real(cp), intent(in) :: cosmicTimeHo

    real(cp):: zcross, tcrossHo
    real(cp) :: zeq, sqrtomM
    
    sqrtomM = sqrt(flParams%OmegaM)
    zcross = redshift_crossing()
    tcrossHo = cp_cosmic_radtime_normalized(zcross)
   
    if (cosmicTimeHo.lt.tcrossHo) then
       cp_redshift_scalingtime_normalized = cp_redshift_radtime_normalized(cosmicTimeHo,Q=1._cp)
    elseif (cosmicTimeHo.ge.tcrossHo) then
       cp_redshift_scalingtime_normalized = (2._cp/(3._cp*cosmicTimeHo &
            * sqrtomM))**(2._cp/3._cp) - 1._cp
    else
       write(*,*)'tHo= ',cosmicTimeHo
       stop 'cp_redshift_scalingtime_normalized: t screwed!'
    endif
  end function cp_redshift_scalingtime_normalized


 

  function ep_redshift_scalingtime_normalized(cosmicTimeHo)
    implicit none
    real(cep) :: ep_redshift_scalingtime_normalized
    real(cep), intent(in) :: cosmicTimeHo

    real(cep):: zcross, tcrossHo, sqrtomM

    sqrtomM = sqrt(real(flParams%OmegaM,cep))    
    zcross = real(redshift_crossing(),cep)
    tcrossHo = ep_cosmic_radtime_normalized(zcross)
   
    if (cosmicTimeHo.lt.tcrossHo) then
       ep_redshift_scalingtime_normalized = ep_redshift_radtime_normalized(cosmicTimeHo,Q=1._cep)
    elseif (cosmicTimeHo.ge.tcrossHo) then
       ep_redshift_scalingtime_normalized = (2._cep/(3._cep*cosmicTimeHo &
            * sqrtomM))**(2._cep/3._cep) - 1._cep
    else
       stop 'ep_redshift_scalingtime_normalized: t screwed!'
    endif
  end function ep_redshift_scalingtime_normalized



  function brent_redshift_radtime_normalized(cosmicTimeHo)
    use functools, only : brent
    implicit none
    real(cp), intent(in) :: cosmicTimeHo
    real(cp) :: brent_redshift_radtime_normalized

    real(cp) :: lnzp1min, lnzp1max
    real(cp), parameter :: tol = tolfl
    
    lnzp1min = 0._cp
    lnzp1max = log(big)/4._cp
    statbuffer = cosmicTimeHo
    
     brent_redshift_radtime_normalized &
         = exp(brent(lnzp1min,lnzp1max,find_brent_redshift_radtime_normalized,tol)) - 1._cp
    

  end function brent_redshift_radtime_normalized


  function find_brent_redshift_radtime_normalized(lnzp1)
    implicit none
    real(cp) :: find_brent_redshift_radtime_normalized
    real(cp), intent(in) :: lnzp1

    real(cp) :: zp1,tHo
    
    zp1 = exp(lnzp1)
    tHo = cosmic_radtime_normalized(zp1-1._cp)
    
    find_brent_redshift_radtime_normalized = tHo - statbuffer
    
  end function find_brent_redshift_radtime_normalized


  
  
  
  recursive function cp_redshift_radtime_normalized(cosmicTimeHo,Q,y) result(z)
    implicit none
    real(cp), intent(in) :: cosmicTimeHo
    real(cp), intent(in), optional :: Q,y
    real(cp) :: z, dz

    real(cp), parameter :: tol = tolfl
    
    real(cp) :: sqrtomR

    if (present(Q)) then
       sqrtomR = sqrt(flParams%OmegaR*Q)
    else
       sqrtomR = sqrt(flParams%OmegaR)
    endif

    z = 1._cp/sqrt(2._cp*sqrtomR*cosmicTimeHo) - 1._cp
    
#ifdef THERMAL
    if (.not.present(Q)) return
    if (present(y)) then
       dz = 2*(y - z)/(y + z)
       if (abs(dz).lt.tol) return
    endif
    z = cp_redshift_radtime_normalized(cosmicTimeHo,correction_rdof(z),z)
#endif    
    
  end function cp_redshift_radtime_normalized
  



  recursive function ep_redshift_radtime_normalized(cosmicTimeHo,Q,y) result(z)
    implicit none
    real(cep), intent(in) :: cosmicTimeHo
    real(cep), intent(in), optional :: Q,y
    real(cep) :: z, dz

    real(cep), parameter :: tol = tolfl
    
    real(cep) :: sqrtomR

    if (present(Q)) then
       sqrtomR = sqrt(real(flParams%OmegaR,cep)*Q)
    else
       sqrtomR = sqrt(real(flParams%OmegaR,cep))
    endif

    z = 1._cep/sqrt(2._cep*sqrtomR*cosmicTimeHo) - 1._cep
    
#ifdef THERMAL
    if (.not.present(Q)) return
    if (present(y)) then
       dz = 2*(y - z)/(y + z)
       if (abs(dz).lt.tol) return
    endif
    z = ep_redshift_radtime_normalized(cosmicTimeHo,real(correction_rdof(real(z,cp)),cep),z)
#endif    
    
  end function ep_redshift_radtime_normalized
  
  
  
  
!returns z such that Ho t(z)(1+z) = Hotoa  
  function cp_redshift_toa_scalingtime_normalized(Hotoa)
    implicit none
    real(cp) :: cp_redshift_toa_scalingtime_normalized
    real(cp), intent(in) :: Hotoa

    real(cp) :: zcross, Hotoacross, omM

    omM = flParams%OmegaM
    
    zcross = redshift_crossing()
    Hotoacross = cp_cosmic_radtime_normalized(zcross)*(1._cp + zcross)
    
    if (Hotoa.lt.Hotoacross) then
       cp_redshift_toa_scalingtime_normalized = cp_redshift_toa_radtime_normalized(Hotoa,Q=1._cp)
    elseif (Hotoa.ge.Hotoacross) then
       cp_redshift_toa_scalingtime_normalized = 4._cp/(9._cp*Hotoa*Hotoa*omM) - 1._cp
    else
       stop 'cp_redshift_toa_scalingtime_normalized: screwed!'
    endif
    
  end function cp_redshift_toa_scalingtime_normalized



  
  function ep_redshift_toa_scalingtime_normalized(Hotoa)
    implicit none
    real(cep) :: ep_redshift_toa_scalingtime_normalized
    real(cep), intent(in) :: Hotoa

    real(cep) :: zcross, Hotoacross, omM

    omM = real(flParams%OmegaM,cep)
    
    zcross = real(redshift_crossing(),cep)
    Hotoacross = ep_cosmic_radtime_normalized(zcross) * (1._cep + zcross)

    if (Hotoa.lt.Hotoacross) then
       ep_redshift_toa_scalingtime_normalized =  ep_redshift_toa_radtime_normalized(Hotoa,Q=1._cep)
    elseif (Hotoa.ge.Hotoacross) then
       ep_redshift_toa_scalingtime_normalized = 4._cep/(9._cep*Hotoa*Hotoa*omM) - 1._cep
    else
       stop 'ep_redshift_toa_scalingtime_normalized: screwed!'
    endif

  end function ep_redshift_toa_scalingtime_normalized
  



  recursive function cp_redshift_toa_radtime_normalized(Hotoa,Q,y) result(z)
    implicit none
    real(cp), intent(in) :: Hotoa
    real(cp), intent(in), optional :: Q, y
    real(cp) :: z

    real(cp), parameter :: tol = tolfl
    
    real(cp) :: sqrtomR, dz

    if (present(Q)) then
       sqrtomR = sqrt(flParams%OmegaR*Q)
    else
       sqrtomR = sqrt(flParams%OmegaR)
    endif

    z = 1._cp/(2._cp*Hotoa*sqrtomR) - 1._cp
               
#ifdef THERMAL
    if (.not.present(Q)) return
    if (present(y)) then
       dz = 2*(y - z)/(y + z)
       if (abs(dz).lt.tol) return
    endif
    z = cp_redshift_toa_radtime_normalized(Hotoa,correction_rdof(z),z)
#endif    
    
  end function cp_redshift_toa_radtime_normalized

  

  recursive function ep_redshift_toa_radtime_normalized(Hotoa,Q,y) result(z)
    implicit none
    real(cep), intent(in) :: Hotoa
    real(cep), intent(in) , optional :: Q,y
    real(cep) :: z
    
    real(cep), parameter :: tol = tolfl

    real(cep) :: sqrtomR, dz   

    if (present(Q)) then
       sqrtomR = sqrt(real(flParams%OmegaR,cep)*Q)
    else
       sqrtomR = sqrt(flParams%OmegaR)
    endif
        
    z = 1._cep/(2._cep*Hotoa*sqrtomR) - 1._cep
    
#ifdef THERMAL
    if (.not.present(Q)) return
    if (present(y)) then
       dz = 2*(y - z)/(y + z)
       if (abs(dz).lt.tol) return
    endif
    z = ep_redshift_toa_radtime_normalized(Hotoa,real(correction_rdof(real(z,cp)),cep),z)
#endif
    
  end function ep_redshift_toa_radtime_normalized


  
  
  
!returns z such that Ho Chi(z) (1+z) = Hochioa
  function redshift_chioa_normalized(Hochioa)
    use functools, only : brent
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


  

!returns z such that  Ho t(z) (1+z) [Chi(z)/t(z)]^alpha = numval
  function redshift_tchipower_scalingtime_normalized(alpha,numval)
    use functools, only : brent
    implicit none
    real(cp) :: redshift_tchipower_scalingtime_normalized
    real(cp), intent(in) :: alpha,numval

    real(cp) :: lnzp1min, lnzp1max, lnzp1
    real(cp), parameter :: tol = tolfl

    
    lnzp1min = 0._cp
    lnzp1max = log(big)/2._cp
    statbuffer = numval
    statpower = alpha
    
    lnzp1 = brent(lnzp1min,lnzp1max,find_redshift_tchipower_scalingtime_normalized,tol)
    
    redshift_tchipower_scalingtime_normalized = exp(lnzp1) - 1._cp
         
  end function redshift_tchipower_scalingtime_normalized

  function find_redshift_tchipower_scalingtime_normalized(lnzp1)
    implicit none
    real(cp) :: find_redshift_tchipower_scalingtime_normalized
    real(cp), intent(in) :: lnzp1

    real(cp) :: zp1,Hochi,Hot

    if (lnzp1.eq.0._cp) then
       find_redshift_tchipower_scalingtime_normalized = -statbuffer
       return
    endif
    
    zp1 = exp(lnzp1)
    
    Hochi = comoving_distance_normalized(zp1-1._cp)
    Hot = cosmic_scalingtime_normalized(zp1-1._cp)
    
    find_redshift_tchipower_scalingtime_normalized = zp1*Hot*(Hochi/Hot)**statpower - statbuffer
    
  end function find_redshift_tchipower_scalingtime_normalized
  
  
end module flrw
