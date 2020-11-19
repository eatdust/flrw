module cosmospline
  use prec, only : dp
  use bspline, only : dbsnak, dbsint, dbsval
  use flrw
  implicit none


  private


  integer, parameter :: order = 5
  integer, save :: ndata = 0

  real(dp), save :: lnzmin, lnzmax
  real(dp), dimension(:), allocatable :: knots, bcoef


  public check_flrw_splines, initialize_flrw_splines, free_flrw_splines
  public spline_comoving_distance_normalized
  
contains

  
  function check_flrw_splines()
    implicit none
    logical :: check_flrw_splines

    check_flrw_splines = allocated(knots)

    if (check_flrw_splines) then
       if (.not.allocated(bcoef)) stop 'check_flrw_splines: bcoef not found!'
    endif
    
  end function check_flrw_splines

  

  subroutine free_flrw_splines()
    implicit none

    if (check_flrw_splines()) deallocate(knots,bcoef)

  end subroutine free_flrw_splines


  
  subroutine initialize_flrw_splines(zmin,zmax,nevals)
    implicit none
    real(dp), intent(in) :: zmin, zmax
    integer, intent(in) :: nevals

    real(dp), dimension(:), allocatable :: lnz, lnchi
    integer :: i

    
    if (zmin.eq.0._dp) then
       lnzmin = log(epsilon(1._dp))
    else
       lnzmin = log(zmin)
    endif

    lnzmax = log(zmax)*(1._dp+epsilon(1._dp))

    ndata = nevals
    
    allocate(knots(ndata+order))
    allocate(bcoef(ndata))

    allocate(lnz(ndata))
    allocate(lnchi(ndata))

    do i=1,nevals
       lnz(i) = lnzmin + real(i-1,dp)*(lnzmax-lnzmin)/real(nevals-1,dp)
       lnchi(i) = log(comoving_distance_normalized(exp(lnz(i))))
!       print *,'lnz,lnchi',lnz(i),lnchi(i)
    enddo

    call dbsnak(ndata,lnz,order,knots)
    call dbsint(ndata,lnz,lnchi,order,knots,bcoef)

    deallocate(lnz,lnchi)
    
  end subroutine initialize_flrw_splines



  
  function spline_comoving_distance_normalized(z)
    implicit none
    real(dp) :: spline_comoving_distance_normalized
    real(dp), intent(in) :: z

    real(dp) :: lnz

    if (z.lt.exp(lnzmin)) then
       spline_comoving_distance_normalized = z
       return
    endif

    if (z.gt.exp(lnzmax)) then
       write(*,*)'spline_comoving_distance_normalized:'
       write(*,*)'z= zmax= ',z,exp(lnzmax)
       stop
    endif
    
    lnz = log(z)
    
    spline_comoving_distance_normalized = exp(dbsval(lnz,order,knots,ndata,bcoef))
    
  end function spline_comoving_distance_normalized
  

  
end module cosmospline
