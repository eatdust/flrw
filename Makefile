# >>> DESIGNED FOR GMAKE <<<

FC=f95

ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=mpif90
FFLAGS= -O3 -fopenmp -march=native -ffree-line-length-none -D$(ext) -DTHERMAL -DMPI -DMPISCHED
# -Wopenmp-simd -DMPI -DMPISCHED -DTHERMAL
#bugbuster
#FFLAGS = -g -debug variable_locations -inline_debug_info -CB -check all 
#LDFLAGS= -lgsl -lgslcblas
endif
#hyper2F1.o

OBJTH=prec.o bspline.o rdof.o iotools.o functools.o flrw.o

thermal.$(ext) : $(OBJTH) thermal.o
	$(FC) $(FFLAGS) $(OBJTH) thermal.o -o $@ $(LDFLAGS)

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


