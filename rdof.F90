!From Micomegas tabulated values of g* and q*. We convert T,
!originally in GeV into the dimensionless value x = T/Tcmb

module rdof
  use precision, only : fdp
  use bspline, only : dbsnak, dbsint, dbsval, dbsder

  implicit none

  private

  integer, parameter :: dp = fdp
  
  integer, parameter :: order = 3
  integer, save :: ndata = 0
  real(dp), save :: xmin = 0._dp, xmax = 0._dp
  real(dp), save :: zmin = 0._dp, zmax = 0._dp
  
  real(dp), dimension(:), allocatable :: xknots, zknots
  real(dp), dimension(:), allocatable :: gxbcoef, qxbcoef
  real(dp), dimension(:), allocatable :: gzbcoef, qzbcoef
  real(dp), dimension(:), allocatable :: xbcoef, zbcoef
  real(dp), dimension(:), allocatable :: czbcoef, wp1zbcoef
  
!max number of records in precomputed data file
  integer, parameter :: nrecmax = 300
  integer, parameter :: ncols = 3

  real(dp), parameter :: To = 2.725 !K

  real(dp), parameter :: eV = 11604.505 !K
  real(dp), parameter :: GeV = eV*1d9 !K

  real(dp), parameter :: MpInGeV = 2.435323186423786d18 !GeV
  real(dp), parameter :: MpInK = MpInGeV*GeV

  
  logical, parameter :: display = .false.

!define CREATEPP to generate binary preprocessed data in rdof.pp
!from data file

  
  public To, GeV, MpInK, dp
  
  public readump_data, preprocessed_data
  public set_splines, check_splines, free_splines
  public energy_rdof_x, entropy_rdof_x
  public energy_rdof_z, entropy_rdof_z, correction_rdof_z, eos_rdof_z
  public x_rdof, z_rdof
contains




  

  subroutine readump_data(filename,xdata, qdata, gdata, zdata, cdata, wp1data)
    implicit none

    character(len=*), intent(in) :: filename
    real(dp), dimension(:), allocatable, intent(inout) :: xdata, qdata, gdata, zdata, cdata, wp1data

    integer, parameter :: nunit = 110

    integer :: ioerr, nrec,i

    real(dp), dimension(ncols), save :: statbuffer
    real(dp), dimension(nrecmax,ncols) :: buffer

    real(dp) :: qo
    
    write(*,*)'reading g* and q* values table...'


    open(unit=nunit,file=filename,status='old')   

    do i=1,nrecmax
       read(nunit,*,iostat=ioerr) statbuffer         

       buffer(i,1) = statbuffer(1) * GeV/To
       buffer(i,2:3) = statbuffer(2:3)

       if (ioerr.ne.0) exit

    enddo

    close(nunit)

    write(*,*)
    write(*,*)'readump_data:'
    write(*,*)'number of records: ',i-1

    nrec = i-1

    allocate(xdata(nrec),qdata(nrec),gdata(nrec),zdata(nrec), cdata(nrec), wp1data(nrec))
    xdata = buffer(:,1)
    qdata = buffer(:,2)
    gdata = buffer(:,3)

    
    
    ndata = nrec
    xmin = xdata(ndata)
    xmax = xdata(1)

    zdata = redshift_from_entropy(xdata,qdata/qdata(ndata))
    zmin = zdata(ndata)
    zmax = zdata(1)

    cdata = correction_rdof_from_all(xdata,qdata/qdata(ndata),gdata/gdata(ndata))
    wp1data = eosp1_rdof_from_all(xdata,qdata,gdata)
    
#ifdef CREATEPP    
    write(*,*)'readump_data: overwritting rdof.pp, press ENTER to proceed!'
    read(*,*)
    open(unit=nunit,file='rdof.pp',status='unknown')
    do i=1,nrec
       write(nunit,*)'PPDATA(',i,',',xdata(i),',',qdata(i),',',gdata(i) &
            ,',',zdata(i),',',cdata(i),',',wp1data(i),')'
    enddo
    close(nunit)
    write(*,*)'readump_data: new rdof.pp dumped!'
#endif       

  end subroutine readump_data

 
  
  subroutine preprocessed_data(xdata,qdata,gdata,zdata,cdata,wp1data)
    implicit none
    real(dp), dimension(nrecmax,ncols+3) :: buffer
    real(dp), dimension(:), allocatable, intent(inout) :: xdata,qdata,gdata,zdata,cdata,wp1data
    integer :: i,nrec

    buffer = -1._dp

#define PPDATA(fooidx,foox,fooq,foog,fooa,fooc,foow) \
    buffer(fooidx,1) = foox; \
    buffer(fooidx,2) = fooq; \
    buffer(fooidx,3) = foog ; \
    buffer(fooidx,4) = fooa ; \
    buffer(fooidx,5) = fooc ; \
    buffer(fooidx,6) = foow
#ifndef CREATEPP
#include "rdof.pp"
#endif    
#undef PPDATA
 
    do i=1,nrecmax
       if (buffer(i,1).eq.-1._dp) exit
    enddo

    if (display) then
       write(*,*)
       write(*,*)'preprocessed_data:'
       write(*,*)'number of records: ',i-1
    endif

    nrec = i-1

    allocate(xdata(nrec),qdata(nrec),gdata(nrec),zdata(nrec),cdata(nrec),wp1data(nrec))
    xdata = buffer(:,1)
    qdata = buffer(:,2)
    gdata = buffer(:,3)
    zdata = buffer(:,4)
    cdata = buffer(:,5)
    wp1data = buffer(:,6)
    
    ndata = nrec
    xmin = xdata(ndata)
    xmax = xdata(1)
    zmin = zdata(ndata)
    zmax = zdata(1)

    if (display) then
       write(*,*)'zmin= zmax= ',zmin,zmax
       write(*,*)'xmin= xmax= ',xmin,xmax
    endif
       
  end subroutine preprocessed_data



  
  function check_splines()
    implicit none
    logical :: check_splines

    check_splines = allocated(xknots)
  
!sanity checks
    if (check_splines) then      
       if ((.not.allocated(gxbcoef)).or.(.not.allocated(qxbcoef)) &
            .or.(.not.allocated(xbcoef)).or.(.not.allocated(zknots)) &
            .or.(.not.allocated(gzbcoef)).or.(.not.allocated(qzbcoef)) &
            .or.(.not.allocated(czbcoef)).or.(.not.allocated(wp1zbcoef)) &
            .or.(.not.allocated(zbcoef))) &
            stop 'check_splines: bcoefs not found!'
    endif

  end function  check_splines
  


  subroutine free_splines()
    implicit none

    if (check_splines()) then
       deallocate(xknots,gxbcoef,qxbcoef)
       deallocate(zknots,xbcoef,zbcoef,gzbcoef,qzbcoef,czbcoef,wp1zbcoef)
    endif

  end subroutine free_splines



  
  subroutine set_splines()
    implicit none
    character(len=*), parameter :: tablefile = 'HP_B_thg.sav'
    real(dp), dimension(:), allocatable :: xdata, gdata, qdata, zdata, cdata, wp1data

!set it to true to use preprocessed data
    logical, parameter :: usePP = .true.
    
    if (check_splines()) stop 'set_splines: splines already set!'

    if (usePP) then
       call preprocessed_data(xdata, qdata, gdata, zdata, cdata, wp1data)
    else
       call readump_data(tablefile, xdata, qdata, gdata, zdata, cdata, wp1data)
    endif

    allocate(xknots(ndata+order),zknots(ndata+order))
    allocate(gxbcoef(ndata),qxbcoef(ndata),xbcoef(ndata),zbcoef(ndata))
    allocate(gzbcoef(ndata),qzbcoef(ndata),czbcoef(ndata),wp1zbcoef(ndata))

    call dbsnak(ndata,xdata(ndata:1:-1),order,xknots)
    call dbsint(ndata,xdata(ndata:1:-1),gdata(ndata:1:-1),order,xknots,gxbcoef)
    call dbsint(ndata,xdata(ndata:1:-1),qdata(ndata:1:-1),order,xknots,qxbcoef)
    call dbsint(ndata,xdata(ndata:1:-1),zdata(ndata:1:-1),order,xknots,zbcoef)
    
    call dbsnak(ndata,zdata(ndata:1:-1),order,zknots)
    call dbsint(ndata,zdata(ndata:1:-1),xdata(ndata:1:-1),order,zknots,xbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),gdata(ndata:1:-1),order,zknots,gzbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),qdata(ndata:1:-1),order,zknots,qzbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),cdata(ndata:1:-1),order,zknots,czbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),wp1data(ndata:1:-1),order,zknots,wp1zbcoef)
    
    deallocate(xdata,gdata,qdata,zdata,cdata,wp1data)

  end subroutine set_splines



  
!returns a/ao assuming entropy conservation
  elemental function redshift_from_entropy(x,qoverqo)
    implicit none
    real(dp), intent(in) :: x,qoverqo
    real(dp) :: redshift_from_entropy
    real(dp), parameter :: onethird = 1._dp/3._dp

    redshift_from_entropy = (qoverqo)**(onethird) * x - 1._dp
        
  end function redshift_from_entropy



  
  elemental function correction_rdof_from_all(x,qoverqo,govergo)
    implicit none
    real(dp), intent(in) :: x, qoverqo, govergo
    real(dp) :: correction_rdof_from_all
    real(dp), parameter :: mfourthird = -4._dp/3._dp
    
    correction_rdof_from_all = govergo * (qoverqo)**mfourthird

  end function correction_rdof_from_all
  


  elemental function eosp1_rdof_from_all(x,q,g)
    implicit none
    real(dp), intent(in) :: x, q, g
    real(dp) :: eosp1_rdof_from_all
    real(dp), parameter :: fourthird = 4._dp/3._dp
!w+1
    eosp1_rdof_from_all = fourthird * q/g

  end function eosp1_rdof_from_all

  
  
  recursive function energy_rdof_x(x) result(g)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: g
    
    if (x.lt.xmin) then
       g = energy_rdof_x(xmin)
    elseif (x.gt.xmax) then
       g = energy_rdof_x(xmax)
    else
       g = dbsval(x,order,xknots,ndata,gxbcoef)
    end if
    
       
  end function energy_rdof_x
  


  

  
  recursive function entropy_rdof_x(x) result(q)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: q

    if (x.lt.xmin) then
       q = entropy_rdof_x(xmin)
    elseif (x.gt.xmax) then
       q = entropy_rdof_x(xmax)
    else
       q = dbsval(x,order,xknots,ndata,qxbcoef)
    end if

  end function entropy_rdof_x



  
  
  recursive function x_rdof(z) result(x)
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: x

    if (z.lt.zmin) then
       x = x_rdof(zmin) * (1+z)/(1+zmin)
    elseif (z.gt.zmax) then
       x = x_rdof(zmax) * (1+z)/(1+zmax)
    else
       x = dbsval(z, order, zknots, ndata, xbcoef)
    endif

  end function x_rdof



  recursive function z_rdof(x) result(z)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: z

    if (x.lt.xmin) then
       z = z_rdof(xmin) * (1+x)/(1+xmin)
    elseif (x.gt.xmax) then
       z = z_rdof(xmax) * (1+x)/(1+xmax)
    else
       z = dbsval(x, order, xknots, ndata, zbcoef)
    endif

  end function z_rdof

  
  

  recursive function energy_rdof_z(z) result(g)
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: g
    
    if (z.lt.zmin) then
       g = energy_rdof_z(zmin)
    elseif (z.gt.zmax) then
       g = energy_rdof_z(xmax)
    else
       g = dbsval(z,order,zknots,ndata,gzbcoef)
    end if
    
       
  end function energy_rdof_z
  

  
  
  
  recursive function entropy_rdof_z(z) result(q)
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: q

    if (z.lt.zmin) then
       q = entropy_rdof_z(zmin)
    elseif (z.gt.zmax) then
       q = entropy_rdof_z(zmax)
    else
       q = dbsval(z, order, zknots, ndata, qzbcoef)
    end if
    
  end function entropy_rdof_z
  


  recursive function correction_rdof_z(z) result(c)
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: c

    if (z.lt.zmin) then
       c = correction_rdof_z(zmin)
    elseif (z.gt.zmax) then
       c = correction_rdof_z(zmax)
    else
       c = dbsval(z, order, zknots, ndata, czbcoef)
    end if

  end function correction_rdof_z
  

  
  recursive function eos_rdof_z(z) result(w)
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: w

!w becomes crazy under this value, possible due to e+e-
!annihilation. It is computed from s=(rho+P)/T, the chemical potential
!is not included. Non-conservation of the number of particles breaks
!this calculation
    real(dp) :: zee = 4.7955d9
    
    if (z.lt.zee) then
       w = eos_rdof_z(zee)
    elseif (z.gt.zmax) then
       w = eos_rdof_z(zmax)
    else
       w = dbsval(z, order, zknots, ndata, wp1zbcoef) - 1._dp
    end if

  end function eos_rdof_z
  
end module rdof
