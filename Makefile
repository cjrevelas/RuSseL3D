#######################################################################
#                  Makefile of SCF-FEM code                           #
#######################################################################

MAKE_MPI_RUN=1
MAKE_PRODUCTION_RUN=0

PROD_OPTIONS=
DEBUG_OPTIONS= -DDEBUG_OUTPUTS #-DPRINT_AFULL
BOTH_OPTIONS= #-DVARIABLE_DS_SCHEME# -DMUMPS_REPORT

# MPI/MUMPS SECTION
ifeq ($(MAKE_MPI_RUN),0)
MPI_OPTIONS=
topdir = /home/cjrevelas/MUMPS/mumps_serial
else
MPI_OPTIONS = -DUSE_MPI
#topdir = /home/cjrevelas/MUMPS/mumps_par
topdir = /home/asgouros/TrioStountzes/MUMPS/TEMP_25Aug2019_MUMPS_5.2.1_PAR
endif


libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

# Flags of production and debug runs
FC_PROD=-O3
FC_DEBUG=-O0 -g -fcheck=all -Wall
#Other debug flags: -Wextra -Wno-unused-dummy-argument

# Choose between PROD and DEBUG run
ifeq ($(MAKE_PRODUCTION_RUN),0)
CPPFLAGS=$(DEBUG_OPTIONS) $(BOTH_OPTIONS) $(MPI_OPTIONS)
FCFLAGS=$(FC_DEBUG) -cpp $(CPPFLAGS)
else
CPPFLAGS=$(PROD_OPTIONS) $(BOTH_OPTIONS) $(MPI_OPTIONS)
FCFLAGS=$(FC_PROD) -cpp $(CPPFLAGS)
endif

LIBFS=#-lstdc++ #-lm

MODULES = xdata_mod.o constants_mod.o kcw_mod.o fhash_mod.o \
          mpistuff_mod.o error_handing_mod.o write_helper_mod.o
OBJECTS = matrix_assemble.o part_fun_phi.o \
	  scfinout.o simpsonkoef.o quadinterp_koef.o spat_3d.o \
          surf_pot.o tetshp.o qprint.o mesh_io_3d.o gauss_3d.o \
          edwards_film_fem.o adh_ten.o main.o mumps_sub.o convolution.o  

.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBFS)  $*.f90

CMD=fem_3d.exe
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS)
	   $(FL) -o $(CMD) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

mumps_sub.o :
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/src -cpp $(CPPFLAGS) -c $*.f90 $(OUTF)$*.o

fhash_modules: fhash.f90 fhash_modules.f90
	$(FC) $(FFLAGS) -c fhash_modules.f90

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) *.o *.mod

cleaner:
	$(RM) *.o *.mod *.exe *.out.txt *.bin fort.* 

test:
	./TEST_INTEGRITY/test_integrity.sh
