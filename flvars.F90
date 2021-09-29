module flvars
  implicit none

  private

  integer, parameter :: cp = kind(1._8)
  integer, parameter :: cep = kind(1._16)

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

  real(cp), parameter :: tolfl = 100*epsilon(1._cp)
  real(cp), parameter :: big = epsilon(1._cp)*huge(1._cp)

  public cp, cep, flbgparam
  public c, Mpc, My, LPl
  public tolfl, big
  
  public flParams, correction_rdof, entropy_correction_rdof

contains

  
#ifndef THERMAL
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
!for hiding possible type differences

  function correction_rdof(z)
    use rdof, only : dp, correction_rdof_z
    implicit none
    real(cp), intent(in) :: z
    real(cp) :: correction_rdof
    
    correction_rdof = correction_rdof_z(real(z,dp))
    
  end function correction_rdof
  

  
  function entropy_correction_rdof(z)
    use rdof, only : dp, entropy_rdof_z
    real(cp) :: entropy_correction_rdof
    real(cp), intent(in) :: z
    real(dp) :: redshift

    redshift = real(z,dp)

    entropy_correction_rdof = (real(entropy_rdof_z(redshift),cp)/flParams%qo)
    
  end function entropy_correction_rdof

#endif  

end module flvars
