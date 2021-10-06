module flapprox
!exact, but particular, solutions of the FLRW equations (pure rad,
!pure mat, mixture rad/mat, instantaneous transition)
  use flvars, only : cp, cep, flbgparam
  use flvars, only : flParams
  use flvars, only : c, Mpc,My,LPl
  use flvars, only : tolfl, big
  use flvars, only : flParams, correction_rdof, entropy_correction_rdof
  
  implicit none
  
  private

      
  interface cosmic_puremattime
     module procedure cp_cosmic_puremattime, ep_cosmic_puremattime
  end interface cosmic_puremattime

  interface cosmic_puremattime_normalized
     module procedure cp_cosmic_puremattime_normalized, ep_cosmic_puremattime_normalized
  end interface cosmic_puremattime_normalized

  interface redshift_puremattime_normalized
     module procedure cp_redshift_puremattime_normalized, ep_redshift_puremattime_normalized
  end interface redshift_puremattime_normalized
  
  interface redshift_toa_puremattime_normalized
     module procedure cp_redshift_toa_puremattime_normalized, ep_redshift_toa_puremattime_normalized
  end interface redshift_toa_puremattime_normalized
  
  interface cosmic_instmattime
     module procedure cp_cosmic_instmattime, ep_cosmic_instmattime
  end interface cosmic_instmattime

  interface cosmic_instmattime_normalized
     module procedure cp_cosmic_instmattime_normalized, ep_cosmic_instmattime_normalized
  end interface cosmic_instmattime_normalized

  interface redshift_instmattime_normalized
     module procedure cp_redshift_instmattime_normalized, ep_redshift_instmattime_normalized
  end interface redshift_instmattime_normalized
  
  interface redshift_toa_instmattime_normalized
     module procedure cp_redshift_toa_instmattime_normalized, ep_redshift_toa_instmattime_normalized
  end interface redshift_toa_instmattime_normalized

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
     module procedure cp_redshift_toa_radtime_normalized, ep_redshift_toa_radtime_normalized
  end interface redshift_toa_radtime_normalized

!INSTMAT stands for matter era with a smooth transition RAD->MAT at zeq
#ifdef INSTMAT
  
  interface cosmic_mattime
     module procedure cp_cosmic_instmattime, ep_cosmic_instmattime
  end interface cosmic_mattime

  interface cosmic_mattime_normalized
     module procedure cp_cosmic_instmattime_normalized, ep_cosmic_instmattime_normalized
  end interface cosmic_mattime_normalized

  interface redshift_mattime_normalized
     module procedure cp_redshift_instmattime_normalized, ep_redshift_instmattime_normalized 
  end interface redshift_mattime_normalized

  interface redshift_toa_mattime_normalized
     module procedure cp_redshift_toa_instmattime_normalized, ep_redshift_toa_instmattime_normalized
  end interface redshift_toa_mattime_normalized

!PUREMAT default behaviour  
#else
  
  interface cosmic_mattime
     module procedure cp_cosmic_puremattime, ep_cosmic_puremattime
  end interface cosmic_mattime

   interface cosmic_mattime_normalized
     module procedure cp_cosmic_puremattime_normalized, ep_cosmic_puremattime_normalized
  end interface cosmic_mattime_normalized

  interface redshift_mattime_normalized
     module procedure cp_redshift_puremattime_normalized, ep_redshift_puremattime_normalized 
  end interface redshift_mattime_normalized

  interface redshift_toa_mattime_normalized
     module procedure cp_redshift_toa_puremattime_normalized, ep_redshift_toa_puremattime_normalized
  end interface redshift_toa_mattime_normalized

#endif  
  
  real(cp), save :: statbuffer
!$omp threadprivate(statbuffer)  

!pure matter era
  public cosmic_puremattime, cosmic_puremattime_normalized
  public redshift_puremattime_normalized
  public redshift_toa_puremattime_normalized
  public conformal_puremattime_normalized
  
!matter era after instantaneous transition
  public cosmic_instmattime, cosmic_instmattime_normalized
  public conformal_instmattime_normalized
  public redshift_instmattime_normalized
  public redshift_toa_instmattime_normalized

!alias of either one of the other
  public cosmic_mattime, cosmic_mattime_normalized
  public redshift_mattime_normalized, redshift_toa_mattime_normalized
  
!pure radiation era  
  public cosmic_radtime, cosmic_radtime_normalized  
  public conformal_radtime_equality_normalized
  public redshift_radtime_normalized 
  public redshift_toa_radtime_normalized


  
!mixture matter and radiation
  public conformal_radmattime_normalized
  public redshift_conformal_radmattime_normalized
  


  
contains


!pure matter era
  function cp_cosmic_puremattime(redshift)
    implicit none
    real(cp) :: cp_cosmic_puremattime
    real(cp), intent(in) :: redshift

    cp_cosmic_puremattime = cp_cosmic_puremattime_normalized(redshift)/flParams%hubbleToday
    
  end function cp_cosmic_puremattime


  
!pure matter era
  function ep_cosmic_puremattime(redshift)
    implicit none
    real(cep) :: ep_cosmic_puremattime
    real(cep), intent(in) :: redshift

    ep_cosmic_puremattime = ep_cosmic_puremattime_normalized(redshift)/real(flParams%hubbleToday,cep)
    
  end function ep_cosmic_puremattime

  
  
!pure matter era  
  function cp_cosmic_puremattime_normalized(redshift)
    implicit none

    real(cp) :: cp_cosmic_puremattime_normalized
    real(cp), intent(in) :: redshift

    real(cp) :: sqrtomM

    sqrtomM = sqrt(flParams%OmegaM)
    
    cp_cosmic_puremattime_normalized = 2._cp/(3._cp &
         * sqrtomM)/(1._cp + redshift)**(3._cp/2._cp)
    
  end function cp_cosmic_puremattime_normalized


!pure matter era
  function ep_cosmic_puremattime_normalized(redshift)
    implicit none

    real(cep) :: ep_cosmic_puremattime_normalized
    real(cep), intent(in) :: redshift

    real(cep) :: sqrtomM
    
    sqrtomM = sqrt(real(flParams%OmegaM,cep))
    
    ep_cosmic_puremattime_normalized = 2._cep/(3._cep &
         * sqrtomM)/(1._cep + redshift)**(3._cep/2._cep)
    
  end function ep_cosmic_puremattime_normalized
  

  function cp_redshift_puremattime_normalized(cosmicTimeHo)
    implicit none
    real(cp) :: cp_redshift_puremattime_normalized
    real(cp), intent(in) :: cosmicTimeHo

    real(cp) :: sqrtomM

    sqrtomM = sqrt(flParams%OmegaM)

    cp_redshift_puremattime_normalized = (2._cp/(3._cp*cosmicTimeHo &
         * sqrtomM))**(2._cp/3._cp) - 1._cp

  end function cp_redshift_puremattime_normalized


  function ep_redshift_puremattime_normalized(cosmicTimeHo)
    implicit none
    real(cep) :: ep_redshift_puremattime_normalized
    real(cep), intent(in) :: cosmicTimeHo

    real(cep) :: sqrtomM

    sqrtomM = sqrt(real(flParams%OmegaM,cep))

    ep_redshift_puremattime_normalized = (2._cep/(3._cep*cosmicTimeHo &
         * sqrtomM))**(2._cep/3._cep) - 1._cep

  end function ep_redshift_puremattime_normalized


!returns z such that Ho t(z)(1+z) = Hotoa  
  function cp_redshift_toa_puremattime_normalized(Hotoa)
    implicit none
    real(cp) :: cp_redshift_toa_puremattime_normalized
    real(cp), intent(in) :: Hotoa

    real(cp) :: omM

    omM = flParams%OmegaM
        
    cp_redshift_toa_puremattime_normalized = 4._cp/(9._cp*Hotoa*Hotoa*omM) - 1._cp
   
    
  end function cp_redshift_toa_puremattime_normalized

  function ep_redshift_toa_puremattime_normalized(Hotoa)
    implicit none
    real(cep) :: ep_redshift_toa_puremattime_normalized
    real(cep), intent(in) :: Hotoa

    real(cep) :: omM

    omM = real(flParams%OmegaM,cep)
        
    ep_redshift_toa_puremattime_normalized = 4._cep/(9._cep*Hotoa*Hotoa*omM) - 1._cep
   
  end function ep_redshift_toa_puremattime_normalized
  



  function conformal_puremattime_normalized(z)
    implicit none
    real(cp) :: conformal_puremattime_normalized
    real(cp), intent(in) :: z

    conformal_puremattime_normalized = 2._cp / sqrt(flParams%OmegaM * (1._cp + z))

  end function conformal_puremattime_normalized





  
  
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




  
!for test, inversion of redshift from time using brent
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


  
  
!the clever way, using recursion  
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


     



  
!returns the radiation-era calHo * conformal time at which rhomat=rhorad or
!z=zeq. This is the shift for matter era solutions needed by
!instantaneous transition radiation to matter.
  function conformal_radtime_equality_normalized()
    implicit none
    real(cp) :: conformal_radtime_equality_normalized

    conformal_radtime_equality_normalized = sqrt(flParams%OmegaR)/flParams%OmegaM
    

  end function conformal_radtime_equality_normalized

  

!returns the radiation-era Ho t at which rhomat=rhorad or z=zeq. This is the shift for matter era solutions with instantaenous transition radiation to matter
  function cosmic_radtime_equality_normalized()
    implicit none
    real(cp) :: cosmic_radtime_equality_normalized

    cosmic_radtime_equality_normalized = 0.5_cp * flParams%OmegaR**1.5_cp/flParams%OmegaM**2
    

  end function cosmic_radtime_equality_normalized



!matter era after instantaneous transition t -> t + 1/3 teq


  function cp_cosmic_instmattime(redshift)
    implicit none

    real(cp) :: cp_cosmic_instmattime
    real(cp), intent(in) :: redshift

    cp_cosmic_instmattime = cp_cosmic_instmattime_normalized(redshift) /flParams%hubbleToday
   
  end function cp_cosmic_instmattime


  function ep_cosmic_instmattime(redshift)
    implicit none

    real(cep) :: ep_cosmic_instmattime
    real(cep), intent(in) :: redshift

    ep_cosmic_instmattime = ep_cosmic_instmattime_normalized(redshift) /flParams%hubbleToday
   
  end function ep_cosmic_instmattime

  
  

  function cp_cosmic_instmattime_normalized(redshift)
    implicit none

    real(cp) :: cp_cosmic_instmattime_normalized
    real(cp), intent(in) :: redshift

    cp_cosmic_instmattime_normalized = cp_cosmic_puremattime_normalized(redshift) &
         - cosmic_radtime_equality_normalized()/3._cp
   
  end function cp_cosmic_instmattime_normalized

  
  function ep_cosmic_instmattime_normalized(redshift)
    implicit none

    real(cep) :: ep_cosmic_instmattime_normalized
    real(cep), intent(in) :: redshift

    ep_cosmic_instmattime_normalized = ep_cosmic_puremattime_normalized(redshift) &
         - real(cosmic_radtime_equality_normalized(),cep)/3._cep
   
  end function ep_cosmic_instmattime_normalized
  


  function cp_redshift_instmattime_normalized(cosmicTimeHo)
    implicit none
    real(cp) :: cp_redshift_instmattime_normalized
    real(cp), intent(in) :: cosmicTimeHo

    real(cp) :: sqrtomM, teqHo

    sqrtomM = sqrt(flParams%OmegaM)
    teqHo = cosmic_radtime_equality_normalized()

    cp_redshift_instmattime_normalized = (2._cp/((3._cp*cosmicTimeHo + teqHo) &
         * sqrtomM))**(2._cp/3._cp) - 1._cp

  end function cp_redshift_instmattime_normalized


  function ep_redshift_instmattime_normalized(cosmicTimeHo)
    implicit none
    real(cep) :: ep_redshift_instmattime_normalized
    real(cep), intent(in) :: cosmicTimeHo

    real(cep) :: sqrtomM, teqHo

    sqrtomM = sqrt(real(flParams%OmegaM,cep))
    teqHo = real(cosmic_radtime_equality_normalized(),cep)
    
    ep_redshift_instmattime_normalized = (2._cep/((3._cep*cosmicTimeHo + teqHo) &
         * sqrtomM))**(2._cep/3._cep) - 1._cep

  end function ep_redshift_instmattime_normalized


!returns z such that Ho t(z)(1+z) = Hotoa  
  function cp_redshift_toa_instmattime_normalized(Hotoa)
    implicit none
    real(cp) :: cp_redshift_toa_instmattime_normalized
    real(cp), intent(in) :: Hotoa

    real(cp) :: a,b,x,det,teqHo
    real(cp), parameter :: onethird = 1._cp/3._cp
    real(cp), parameter :: twothird = 2._cp/3._cp

    teqHo = cosmic_radtime_equality_normalized()
    
    a = 3._cp*Hotoa/teqHo
    b = 2._cp/(teqHo*sqrt(flParams%OmegaM))

!solution of x^3 + ax -b =0

    det = 9._cp*b+sqrt(12._cp*a**3 + 81._cp*b*b)
    
    x = 2._cp**onethird* (det**twothird - 2._cp*3._cp**onethird*a) &
         / (6._cp**twothird * det**onethird)
    
    cp_redshift_toa_instmattime_normalized = x*x - 1._cp
       
  end function cp_redshift_toa_instmattime_normalized

  
  function ep_redshift_toa_instmattime_normalized(Hotoa)
    implicit none
    real(cep) :: ep_redshift_toa_instmattime_normalized
    real(cep), intent(in) :: Hotoa

    real(cep) :: a,b,x,det,teqHo
    real(cep), parameter :: onethird = 1._cep/3._cep
    real(cep), parameter :: twothird = 2._cep/3._cep

    teqHo = real(cosmic_radtime_equality_normalized(),cep)
    
    a = 3._cep*Hotoa/teqHo
    b = 2._cep/(teqHo*sqrt(flParams%OmegaM))

!solution of x^3 + ax -b =0

    det = 9._cep*b+sqrt(12._cep*a**3 + 81._cep*b*b)
    
    x = 2._cep**onethird* (det**twothird - 2._cep*3._cep**onethird*a) &
         / (6._cep**twothird * det**onethird)
    
    ep_redshift_toa_instmattime_normalized = x*x - 1._cep
       
  end function ep_redshift_toa_instmattime_normalized
  


  
  
  function conformal_instmattime_normalized(z)
    implicit none
    real(cp) :: conformal_instmattime_normalized
    real(cp), intent(in) :: z
   
    conformal_instmattime_normalized = 2._cp / sqrt(flParams%OmegaM * (1._cp + z)) &
         - conformal_radtime_equality_normalized()

  end function conformal_instmattime_normalized

  
  
  
  

!mixture of radiation and matter

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














 
  
end module flapprox
