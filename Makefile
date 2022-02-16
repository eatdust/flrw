# >>> DESIGNED FOR GMAKE <<<


ext=$(shell uname | cut -c1-3)

FC=gfortran
FFLAGS= -O -fopenmp -ffree-line-length-none -D$(ext) -DTHERMAL
#-DTHERMAL -DCREATEPP
SUNDIALS=yes
LDFLAGS=


ifneq ($(SUNDIALS),)
FFLAGS += -DSUNDIALS
LDFLAGS += -lsundials_cvode -lsundials_fcvode_mod
OBJSUN+=sundials.o
endif

OBJTH=precision.o bspline.o rdof.o iotools.o functools.o $(OBJSUN) flvars.o flapprox.o flrw.o

thermalmain.$(ext) : $(OBJTH) thermalmain.o
	$(FC) $(FFLAGS) $(OBJTH) thermalmain.o -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(FFLAGS) -c $*.F90

%.o: %.f03
	$(FC) $(FFLAGS) -c $*.f03

%.o: %.F03
	$(FC) $(FFLAGS) -c $*.F03

%.o: %.F08
	$(FC) $(FFLAGS) -c $*.F08

clean:
	rm *.$(ext) *.o *.mod


