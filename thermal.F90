program thermal
  use prec, only : dp
  use iotools
  use rdof
  use flrw

  implicit none

  real(dp) :: x
  real(dp) :: lnxmin, lnxmax
  real(dp), dimension(:), allocatable :: xdata, qdata, gdata, adata, cdata

  real(dp) :: zp1,z,a
  real(dp) :: lnzp1min, lnzp1max
  
  real(dp) :: qstar, gstar, zcross

  real(dp) :: tHo, tscalHo, tradHo, tmatHo
  
  integer :: npoints
  integer :: i

  
!  call readump_data('HP_B_thg.sav',xdata, qdata, gdata, adata, cdata)
!  stop

  call preprocessed_data(xdata, qdata, gdata, adata, cdata)

!  call set_splines()

  
#ifdef THERMAL
  call set_fiducial_flparams()
#else
  call set_splines()
#endif
  
  call delete_file('gandqstar_x.dat')

  npoints = 1000
  lnxmin = 0._dp
  lnxmax = 50._dp
  
  do i = 1, npoints
     x = exp(lnxmin + real(i-1,dp)*(lnxMax - lnxMin) &
          /real(npoints-1,dp))

     qstar = entropy_rdof_x(x)
     gstar = energy_rdof_x(x)
     
     call livewrite('gandqstar_x.dat',x*To/GeV,qstar,gstar)
  enddo
  
  call delete_file('gandqstar_z.dat')
  call delete_file('hubble_z.dat')
  call delete_file('rdofcorr.dat')
  call delete_file('times.dat')
  
  lnzp1min = 0.1_dp
  lnzp1max = 40._dp

  zcross = redshift_crossing()
  print *,'zcross',zcross,redshift_equality()

  read(*,*)
  !  print *,'test',cosmic_mattime_normalized(zcross),cosmic_radtime_normalized(zcross) &
!, cosmic_scalingtime_normalized(zcross)
  
  do i = 1, npoints
     zp1 = exp(lnzp1min + real(i-1,dp)*(lnzp1Max - lnzp1Min) &
          /real(npoints-1,dp))
     z = zp1 - 1._dp
     gstar = energy_rdof_z(z)
     qstar = entropy_rdof_z(z)
     x = x_rdof(z)

     tHo = cosmic_time_normalized(z)
     tscalHo = cosmic_scalingtime_normalized(z)
     tmatHo = cosmic_mattime_normalized(z)
     tradHo = cosmic_radtime_normalized(z)

     a = 1._dp/(1._dp + z)

!     x = redshift_radtime_normalized(tradHo,Q=1._dp)

!     x = redshift_toa_radtime_normalized(tradHo/a,Q=1._dp)
!     print *,'scaling z= x-z= ',z,(x-z)/(x+z)*2
     
     
     call livewrite('times.dat',z,tHo,tscalHo,tmatHo,tradHo)
     
     call livewrite('gandqstar_z.dat',x*To/GeV,gstar,qstar)

     call livewrite('hubble_z.dat',z,hubble_normalized(z)/zp1/zp1)
#ifdef THERMAL     
     call livewrite('rdofcorr.dat',z,correction_rdof(z))
#endif     
     
  enddo
     

     

end program thermal
