!From Micomegas tabulated values of g* and q*. We convert T,
!originally in GeV into the dimensionless value x = T/Tcmb

module rdof
  use prec, only : dp
  use bspline, only : dbsnak, dbsint, dbsval, dbsder

  implicit none

  private
  
  integer, parameter :: order = 3
  integer, save :: ndata = 0
  real(dp), save :: xmin = 0._dp, xmax = 0._dp
  real(dp), save :: zmin = 0._dp, zmax = 0._dp
  
  real(dp), dimension(:), allocatable :: xknots, zknots
  real(dp), dimension(:), allocatable :: gxbcoef, qxbcoef
  real(dp), dimension(:), allocatable :: gzbcoef, qzbcoef, czbcoef, xbcoef

!max number of records in precomputed data file
  integer, parameter :: nrecmax = 300
  integer, parameter :: ncols = 3

  real(dp), parameter :: To = 2.725 !K
  real(dp), parameter :: eV = 11604.505 !K
  real(dp), parameter :: GeV = eV*1d9 !K

  logical, parameter :: display = .false.
!set it to true to generate binary preprocessed data in rdof.pp
!from data file
  logical, parameter :: createPP = .false.  
  
  public To, GeV, dp
  
  public readump_data, preprocessed_data
  public set_splines, check_splines, free_splines
  public energy_rdof_x, entropy_rdof_x
  public x_rdof, energy_rdof_z, entropy_rdof_z, correction_rdof_z
  
contains




  

  subroutine readump_data(filename,xdata, qdata, gdata, zdata, cdata)
    implicit none

    character(len=*), intent(in) :: filename
    real(dp), dimension(:), allocatable, intent(inout) :: xdata, qdata, gdata, zdata, cdata

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

    allocate(xdata(nrec),qdata(nrec),gdata(nrec),zdata(nrec), cdata(nrec))
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
    
    
    if (createPP) then
       open(unit=nunit,file='rdof.pp',status='new')
       do i=1,nrec
          write(nunit,*)'PPDATA(',i,',',xdata(i),',',qdata(i),',',gdata(i) &
               ,',',zdata(i),',',cdata(i),')'
       enddo
       close(nunit)
    endif

  end subroutine readump_data

 
  
  subroutine preprocessed_data(xdata,qdata,gdata,zdata,cdata)
    implicit none
    real(dp), dimension(nrecmax,ncols+2) :: buffer
    real(dp), dimension(:), allocatable, intent(inout) :: xdata,qdata,gdata,zdata,cdata
    integer :: i,nrec

    buffer = -1._dp
    
#define PPDATA(fooidx,foox,fooq,foog,fooa,fooc) \
    buffer(fooidx,1) = foox; \
    buffer(fooidx,2) = fooq; \
    buffer(fooidx,3) = foog ; \
    buffer(fooidx,4) = fooa ; \
    buffer(fooidx,5) = fooc
#include "rdof.pp"
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

    allocate(xdata(nrec),qdata(nrec),gdata(nrec),zdata(nrec),cdata(nrec))
    xdata = buffer(:,1)
    qdata = buffer(:,2)
    gdata = buffer(:,3)
    zdata = buffer(:,4)
    cdata = buffer(:,5)
    
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
            .or.(.not.allocated(czbcoef))) stop 'check_splines: bcoefs not found!'
    endif

  end function  check_splines
  


  subroutine free_splines()
    implicit none

    if (check_splines()) then
       deallocate(xknots,gxbcoef,qxbcoef)
       deallocate(zknots,xbcoef,gzbcoef,qzbcoef,czbcoef)
    endif

  end subroutine free_splines



  
  subroutine set_splines()
    implicit none
    character(len=*), parameter :: tablefile = 'HP_B_thg.sav'
    real(dp), dimension(:), allocatable :: xdata, gdata, qdata, zdata, cdata

!set it to true to use preprocessed data
    logical, parameter :: usePP = .true.
    
    if (check_splines()) stop 'set_splines: splines already set!'

    if (usePP) then
       call preprocessed_data(xdata, qdata, gdata, zdata, cdata)
    else
       call readump_data(tablefile, xdata, qdata, gdata, zdata, cdata)
    endif

    allocate(xknots(ndata+order),zknots(ndata+order))
    allocate(gxbcoef(ndata),qxbcoef(ndata),xbcoef(ndata))
    allocate(gzbcoef(ndata),qzbcoef(ndata),czbcoef(ndata))

    call dbsnak(ndata,xdata(ndata:1:-1),order,xknots)
    call dbsint(ndata,xdata(ndata:1:-1),gdata(ndata:1:-1),order,xknots,gxbcoef)
    call dbsint(ndata,xdata(ndata:1:-1),qdata(ndata:1:-1),order,xknots,qxbcoef)

    call dbsnak(ndata,zdata(ndata:1:-1),order,zknots)
    call dbsint(ndata,zdata(ndata:1:-1),xdata(ndata:1:-1),order,zknots,xbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),gdata(ndata:1:-1),order,zknots,gzbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),qdata(ndata:1:-1),order,zknots,qzbcoef)
    call dbsint(ndata,zdata(ndata:1:-1),cdata(ndata:1:-1),order,zknots,czbcoef)
    
    deallocate(xdata,gdata,qdata,zdata,cdata)

  end subroutine set_splines



  
!returns a/ao assuming entropy conservation
  elemental function redshift_from_entropy(x,qoverqo)
    implicit none
    real(dp), intent(in) :: x,qoverqo
    real(dp) :: redshift_from_entropy

    redshift_from_entropy = (qoverqo)**(1._dp/3._dp) * x - 1._dp
        
  end function redshift_from_entropy



  
  elemental function correction_rdof_from_all(x,qoverqo,govergo)
    implicit none
    real(dp), intent(in) :: x, qoverqo, govergo
    real(dp) :: correction_rdof_from_all
    real(dp), parameter :: mfourthird = -4._dp/3._dp
    
    correction_rdof_from_all = govergo * (qoverqo)**mfourthird

  end function correction_rdof_from_all
  


  
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
  
  
end module rdof
