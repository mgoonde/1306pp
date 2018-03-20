# -*- Makefile -*-
SHELL=/bin/sh

############################################
# derived makefile variables
OBJ=$(SRC:%.f90=Obj/%.o)
############################################

default:
	$(MAKE) $(MFLAGS) -C Obj

clean:
	$(MAKE) $(MFLAGS) -C Obj clean
