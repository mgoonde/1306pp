compiler = gfortran
objects = routines.o maintest.o
flags = -O
double = -fdefault-real-8

#all: maintest.x f90nautyinterf.mod hashtest.x event_hash.x cycles.x hashing_sites.x
all: f90nautyinterf.mod hashtest.x event_hash.x cycles.x hashing_sites.x

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

cycles.x: routines.o f90nautyinterf.mod sort_index.o cycles.o ffnautyex1_sestic.o
	$(compiler) $(double) -o cycles.x routines.o cycles.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a

hashing_sites.x: routines.o f90nautyinterf.mod sort_index.o ffnautyex1_sestic.o hashing_sites.o
	$(compiler) $(double) -o hashing_sites.x routines.o interface.o sort_index.o ffnautyex1_sestic.o ./nauty/nauty26r11/nauty.a hashing_sites.o


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f *.x
