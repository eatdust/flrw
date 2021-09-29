module flapprox
  use flvars, only : cp, cep, flbgparam
  use flvars, only : flParams
  use flvars, only : c, Mpc,My,LPl
  use flvars, only : tolfl, big
  use flvars, only : flParams, correction_rdof, entropy_correction_rdof
  
  implicit none
  
  private

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

  interface redshift_radtime_normalized
     module procedure cp_redshift_radtime_normalized, ep_redshift_radtime_normalized
  end interface redshift_radtime_normalized

  interface redshift_toa_radtime_normalized
     module procedure cp_redshift_toa_radtime_normalized
  end interface redshift_toa_radtime_normalized


  real(cp), save :: statbuffer
!$omp threadprivate(statbuffer)  

  
  
  public conformal_radmattime_normalized
  public conformal_radtime_equality_normalized

  public redshift_radtime_normalized
 
  public redshift_toa_radtime_normalized
  public redshift_conformal_radmattime_normalized
  
  public cosmic_mattime_normalized, cosmic_radtime_normalized


contains


!returns the radiation-era calH0 * conformal time at which rhomat=rhorad or
!z=zeq. This is the shift for matter era solutions needed by
!instantaneous transitions radiation matter.
  function conformal_radtime_equality_normalized()
    implicit none
    real(cp) :: conformal_radtime_equality_normalized

    conformal_radtime_equality_normalized = sqrt(flParams%OmegaR)/flParams%OmegaM
    

  end function conformal_radtime_equality_normalized

  
  
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
  
  



  

  recursive function redshift_conformal_radmattime_normalized(etaHo,Q,y) result(z)
    implicit none
    real(cp) :: z
    real(cp), intent(in) :: etaHo
    real(cp), intent(in), optional :: Q,y

    real(cp), parameter :: tol = tolfl

    real(cp) :: sqrtomR, omRoM
    real(cp) :: scaleFactor, dz
    
    if (present(Q)) then
       sqrtomR = sqrt(flParams%OmegaR*Q)
    else
       sqrtomR = sqrt(flParams%OmegaR)
    endif

    
    scaleFactor = 0.25_cp*flParams%OmegaM * etaHo*etaHo + sqrtomR * etaHo
    
    z = 1._cp/scaleFactor - 1._cp
    
#ifdef THERMAL
    if (.not.present(Q)) return
    if (present(y)) then
       dz = 2*(y - z)/(y + z)
       if (abs(dz).lt.tol) return
    endif
    z = redshift_conformal_radmattime_normalized(etaHo,correction_rdof(z),z)
#endif    

    
  end function redshift_conformal_radmattime_normalized






!returns etaHo from redshift, approximated if thermal is on
  function conformal_radmattime_normalized(z)
    implicit none
    real(cp) :: conformal_radmattime_normalized
    real(cp), intent(in) :: z

    real(cp) ::  scaleFactor
    real(cp) :: Q

    Q = 1._cp
    scaleFactor = 1._cp / (1._cp + z)

#ifdef THERMAL
    Q = correction_rdof(z)
#endif
    
    conformal_radmattime_normalized = 2._cp * scaleFactor &
         / ( sqrt(flParams%OmegaR*Q + flParams%OmegaM*scaleFactor) &
         + sqrt(flParams%OmegaR*Q) ) 
        
  end function conformal_radmattime_normalized









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


  
  
 
  
end module flapprox
