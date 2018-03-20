compiler = gfortran
objects = routines.o maintest.o
flags = -O
double = -fdefault-real-8

all: maintest.x hashtest.x

hashtest.x : routines.o hashtest.o
	$(compiler) $(double) -o hashtest.x routines.mod routines.o hashtest.o

maintest.x : $(objects)
	$(compiler) $(double) -o maintest.x routines.mod $(objects)

%.o: %.f90
	$(compiler) $(double) ${flags} -c $<

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f *.x
