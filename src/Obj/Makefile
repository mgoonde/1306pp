# -*- Makefile -*-
SHELL=/bin/sh

EXE=../KMC.x

F90=gfortran-mp-7

FFLAGS= -Wall -O -fdefault-real-8 $(FDEBUG) $(LSHARED)
FDEBUG= -g -fbounds-check
LSHARED= -cpp -ldl #-DSHARED_LIB 

all: $(EXE)

VPATH=../
S90 = lev00_tools.f90 \
      lev00_new_types.f90 \
      lev00_lib_dlopen.f90 \
      sys02_adsrbdesrb.f90 \
      lev00_interface_lib.f90 \
      lev01_algorithms.f90 \
      lev01_kernel.f90 \
      lev02_system_init.f90 \
      lev02_system_rate.f90 \
      lev02_neig_list.f90 

O90=$(S90:.f90=.o)

$(EXE): $(O90)
	$(F90) $(FFLAGS) $(O90) -o $(EXE)

%.o: %.f90
	$(F90) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod $(EXE)
