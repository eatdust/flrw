program thermalmain
  use iotools
  use rdof
  use flvars
  use flrw
  use flapprox

  implicit none
  
  real(dp) :: x
  real(dp) :: lnxmin, lnxmax
  real(dp), dimension(:), allocatable :: xdata, qdata, gdata, adata, cdata

  real(dp) :: zp1,z,a
  real(dp) :: lnzp1min, lnzp1max
  
  real(dp) :: qstar, gstar, zeq

  real(dp) :: tHo, tscalHo, tradHo, tpmatHo,timatHo

  real(dp) :: etaHoEx, etaHoApp
  
  integer :: npoints
  integer :: i

  
!  call readump_data('HP_B_thg.sav',xdata, qdata, gdata, adata, cdata)
!  stop

  call preprocessed_data(xdata, qdata, gdata, adata, cdata)

!  call set_splines()

  

  call set_fiducial_flparams()

#ifndef THERMAL
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

  zeq = redshift_equality()
  
  print *,'zeq=',zeq

!  read(*,*)
  print *,'test',cosmic_puremattime_normalized(zeq) &
       ,cosmic_instmattime_normalized(zeq) &
       ,cosmic_radtime_normalized(zeq)
  
  do i = 1, npoints
     zp1 = exp(lnzp1min + real(i-1,dp)*(lnzp1Max - lnzp1Min) &
          /real(npoints-1,dp))
     z = zp1 - 1._dp
     gstar = energy_rdof_z(z)
     qstar = entropy_rdof_z(z)
     x = x_rdof(z)

     tHo = cosmic_time_normalized(z)
     timatHo = cosmic_instmattime_normalized(z)
     tpmatHo = cosmic_puremattime_normalized(z)
     tradHo = cosmic_radtime_normalized(z)

     a = 1._dp/(1._dp + z)

     x = redshift_radtime_normalized(tradHo,Q=1._dp)

!     x = redshift_toa_radtime_normalized(tradHo/a,Q=1._dp)
!     print *,'scaling z= x-z= ',z,(x-z)/(x+z)*2

     etaHoEx = conformal_time_normalized(z)


     call livewrite('eta.dat',z,etaHoEx)
     
     call livewrite('times.dat',z,tHo,tpmatHo,timatHo,tradHo)
     
     call livewrite('gandqstar_z.dat',x*To/GeV,gstar,qstar)

     call livewrite('hubble_z.dat',z,hubble_normalized(z)/zp1/zp1)
#ifdef THERMAL     
     call livewrite('rdofcorr.dat',z,correction_rdof(z))
#endif     
     
  enddo
     

     

end program thermalmain
