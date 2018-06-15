compiler = gfortran
objects = routines.o maintest.o
flags = -O -fbounds-check
double = -fdefault-real-8

#all: maintest.x f90nautyinterf.mod hashtest.x event_hash.x cycles.x hashing_sites.x
#all: f90nautyinterf.mod hashtest.x event_hash.x cycles.x hashing_sites.x
#all: f90nautyinterf.mod event_hash.x 
#all: f90nautyinterf.mod one_event.x parse_syst.x 
#all: f90nautyinterf.mod one_event.x p11.x 
all: f90nautyinterf.mod one_event22.x p22.x

ffnautyex1_sestic.o:
	gcc -c ./ctofinterf/ffnautyex1_sestic.c -I./nauty/nauty26r11/

#maintest.x : $(objects)
#	$(compiler) $(double) -o maintest.x $(objects)

%.o: %.f90
	$(compiler) $(double) ${flags} -c $<

f90nautyinterf.mod: ./ctofinterf/interface.f90
	$(compiler) $(double) ${flags} -c $<

sort_index.o: ./ctofinterf/sort_index.f90
	$(compiler) $(double) ${flags} -c $<

hashtest.x : routines.o hashtest.o f90nautyinterf.mod
	$(compiler) $(double) -o hashtest.x routines.o hashtest.o interface.o

event_hash.x: routines.o f90nautyinterf.mod sort_index.o event_hash.o ffnautyex1_sestic.o
	$(compiler) $(double) -o event_hash.x routines.o event_hash.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a

one_event.x: routines.o f90nautyinterf.mod sort_index.o one_event.o ffnautyex1_sestic.o
	$(compiler) $(double) -o one_event.x routines.o one_event.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a


one_event22.x: routines.o f90nautyinterf.mod sort_index.o one_event22.o ffnautyex1_sestic.o
	$(compiler) $(double) -o one_event22.x routines.o one_event22.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a


#cycles.x: routines.o f90nautyinterf.mod sort_index.o cycles.o ffnautyex1_sestic.o
#	$(compiler) $(double) -o cycles.x routines.o cycles.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a

hashing_sites.x: routines.o f90nautyinterf.mod sort_index.o ffnautyex1_sestic.o hashing_sites.o
	$(compiler) $(double) -o hashing_sites.x routines.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a hashing_sites.o

parse_syst.x: routines.o f90nautyinterf.mod sort_index.o ffnautyex1_sestic.o parse_syst.o
	$(compiler) $(double) -o parse_syst.x routines.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a parse_syst.o

p11.x: routines.o f90nautyinterf.mod sort_index.o ffnautyex1_sestic.o p11.o
	$(compiler) $(double) -o p11.x routines.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a p11.o

p22.x: routines.o f90nautyinterf.mod sort_index.o ffnautyex1_sestic.o p22.o
	$(compiler) $(double) -o p22.x routines.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a p22.o


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f *.x
