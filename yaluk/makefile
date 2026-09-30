# File: makefile
# Purpose: compile and link atp
#
# Change this line (and directories if necessary)
# Change paramet.inc, if necessary
#
SYSTEM = LINUX
FOR = gfortran
AR = ar
# U: store each member's real timestamp in the .a (Debian's ar defaults to D,
# which zeroes timestamps and would make make rebuild everything every time)
ARFLAGS = rvU
#
# OBJECTS is a customized list of modules that are being worked on
# and are present in the current directory

OBJECTS = yaluk_exec.o \
	   		 	yaluk_init.o \
		   		em_field_funct.o \
		   		basic_functions.o

# ------------------------

FFLAGS = -Ofast -ffree-form -ffree-line-length-512 -fopenmp -std=gnu -w -march=native --std=legacy

LIBRARY = libyaluk.a

.f.o:
	$(FOR) -c $(FFLAGS) $<
# The .o files are archive members: make compiles them, adds them to the library
# via the implicit (%): % rule and, being intermediate, deletes them afterwards
$(LIBRARY) : $(LIBRARY)($(OBJECTS))

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(LIBRARY) *.mod
