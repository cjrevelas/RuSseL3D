#######################################################################
#                  Makefile of SCF-FEM code                           #
#######################################################################

MAKE_MPI_RUN=1
MAKE_PRODUCTION_RUN=0

BOTH_OPTIONS=
PROD_OPTIONS=
DEBUG_OPTIONS=-DDEBUG_OUTPUTS

# MPI/MUMPS SECTION
ifeq ($(MAKE_MPI_RUN),0)
CPPFLAGS = -Wno-unused-dummy-argument
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_SERIAL
else
CPPFLAGS = -Wno-unused-dummy-argument -DUSE_MPI
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR
endif

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

# Flags of production and debug runs
PROD=-O3 -cpp $(CPPFLAGS) $(PROD_OPTIONS) $(BOTH_OPTIONS)
DEBUG=-O0 -g -fcheck=all -Wall -cpp $(CPPFLAGS) $(DEBUG_OPTIONS) $(BOTH_OPTIONS)

# Choose between PROD and DEBUG run
ifeq ($(MAKE_PRODUCTION_RUN),0)
FCFLAGS=$(DEBUG)
else
FCFLAGS=$(PROD)
endif

LIBFS=#-lstdc++ #-lm

MODULES =  mdata_mod.o xdata_mod.o constants_mod.o kcw_mod.o fhash_mod.o mpistuff_mod.o
OBJECTS =  matrix_assemble.o  part_fun_phi.o\
	   scfinout.o  simpsonkoef.o spat_3d.o surf_pot.o tetshp.o qprint.o \
           mesh_io_3d.o gauss_3d.o edwards_free_film_fem.o adh_ten.o main.o mumps_sub.o

.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBFS)  $*.f90

CMD=fem_3d.exe
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS) 
	   $(FL) -o $(CMD) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

mumps_sub.o :
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/src -c $*.f90 $(OUTF)$*.o

fhash_modules: fhash.f90 fhash_modules.f90
	$(FC) $(FFLAGS) -c fhash_modules.f90

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) *.o *.mod

cleaner:
	$(RM) *.o *.mod *.exe *.out.txt fort.*

test:
	./TEST_INTEGRITY/test_integrity.sh
