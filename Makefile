# >>> DESIGNED FOR GMAKE <<<


ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=gfortran
FFLAGS= -O -fopenmp -ffree-line-length-none -D$(ext) -DTHERMAL 
#-DTHERMAL -DCREATEPP
endif

OBJTH=precision.o bspline.o rdof.o iotools.o functools.o flvars.o flapprox.o flrw.o

thermalmain.$(ext) : $(OBJTH) thermalmain.o
	$(FC) $(FFLAGS) $(OBJTH) thermalmain.o -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(FFLAGS) -c $*.F90

%.o: %.f03
	$(FC) $(FFLAGS) -c $*.f03

%.o: %.F08
	$(FC) $(FFLAGS) -c $*.F08

clean:
	rm *.$(ext) *.o *.mod


