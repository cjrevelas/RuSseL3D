#######################################################################
#                  Makefile of SCF-FEM code                           #
#######################################################################

MAKE_MPI_RUN=1
MAKE_PRODUCTION_RUN=1

# MPI/MUMPS SECTION
ifeq ($(MAKE_MPI_RUN),1)
CPPFLAGS = -DUSE_MPI -Wno-unused-dummy-argument
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR
else
CPPFLAGS = -Wno-unused-dummy-argument
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_SERIAL
endif

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

# Flags of production and debug runs
PROD=-O3 -cpp $(CPPFLAGS)
DEBUG=-O0 -g -fcheck=all -Wall -cpp $(CPPFLAGS)

# Choose between PROD and DEBUG run
ifeq ($(MAKE_PRODUCTION_RUN),1)
FCFLAGS=$(PROD)
else
FCFLAGS=$(DEBUG)
endif

LIBFS=#-lstdc++ #-lm

MODULES =  mdata_mod.o xdata_mod.o constants_mod.o kcw_mod.o fhash_mod.o mpistuff_mod.o
OBJECTS =  matrix_assemble.o  part_fun_phi.o\
	   scfinout.o  simpsonkoef.o spat_3d.o surf_pot.o tetshp.o qprint.o \
           mesh_io_3d.o gauss_3d.o edwards_free_film_fem.o adh_ten.o main.o mumps_sub.o

.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBFS)  $*.f90

CMD=FEM_3d.exe
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS) 
	   $(FL) -o $(CMD) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

mumps_sub.o :
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/src -c $*.f90 $(OUTF)$*.o

fhash_modules: fhash.f90 fhash_modules.f90
	$(FC) $(FFLAGS) -c fhash_modules.f90

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) *.o *.mod *.x *.exe

cleaner:
	$(RM) *.o *.mod *.x *.exe *.out.txt

test:
	./TEST_INTEGRITY/test_integrity.sh
