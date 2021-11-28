#######################################################################
#                  Makefile of SCF-FEM code                           #
#######################################################################
MAKE_MPI_RUN=0
MAKE_PRODUCTION_RUN=0

PROD_OPTIONS=
DEBUG_OPTIONS= -DDEBUG_OUTPUTS #-DPRINT_AFULL
BOTH_OPTIONS= #-DVARIABLE_DS_SCHEME# -DMUMPS_REPORT

# MPI/MUMPS SECTION
ifeq ($(MAKE_MPI_RUN),0)
MPI_OPTIONS=
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_SERIAL
else
MPI_OPTIONS = -DUSE_MPI
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR
endif

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc                                #defines PLAT, LIBEXT, LIBMUMPS_COMMON, FC   
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

OBJDIR=obj
SRCDIR=src
RUNDIR=run

MODULES = $(join $(OBJDIR),/xdata_mod.o)\
	  $(join $(OBJDIR),/constants_mod.o)\
	  $(join $(OBJDIR),/kcw_mod.o)\
	  $(join $(OBJDIR),/fhash_mod.o)\
          $(join $(OBJDIR),/mpistuff_mod.o)\
	  $(join $(OBJDIR),/error_handing_mod.o)\
  	  $(join $(OBJDIR),/write_helper_mod.o)\

OBJECTS = $(join $(OBJDIR),/matrix_assemble.o)\
          $(join $(OBJDIR),/part_fun_phi.o)\
	  $(join $(OBJDIR),/scfinout.o)\
          $(join $(OBJDIR),/simpsonkoef.o)\
          $(join $(OBJDIR),/quadinterp_koef.o)\
          $(join $(OBJDIR),/spat_3d.o)\
          $(join $(OBJDIR),/surf_pot.o)\
	  $(join $(OBJDIR),/tetshp.o)\
	  $(join $(OBJDIR),/qprint.o)\
	  $(join $(OBJDIR),/mesh_io_3d.o)\
	  $(join $(OBJDIR),/gauss_3d.o)\
          $(join $(OBJDIR),/edwards_film_fem.o)\
	  $(join $(OBJDIR),/adh_ten.o)\
	  $(join $(OBJDIR),/main.o)\
	  $(join $(OBJDIR),/mumps_sub.o)\
	  $(join $(OBJDIR),/convolution.o)

#First rule has been changed from implicit to explicit
OBJTEMP=$(join $(OBJDIR),/%.o)
SRCTEMP=$(join $(SRCDIR),/%.f90)
$(OBJTEMP): $(SRCTEMP)
	$(FC) $(FCFLAGS) $(LIBFS) -Jobj/ -c -o $@ $?

CMD=$(join $(RUNDIR),/fem_3d.exe)
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS)
	   $(FL) -o $(CMD) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

MUMPSTEMP=$(join $(OBJDIR),/mumps_sub.o)
$(MUMPSTEMP):
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/src -Iobj/ -cpp $(CPPFLAGS) -c $(join $(SRCDIR),/mumps_sub.f90) $(OUTF)$*.o

FHASHTEMP=$(join $(OBJDIR),/fhash_modules)
$(FHASHTEMP): $(join $(SRCDIR),/fhash.f90) $(join $(SRCDIR),/fhash_modules.f90)
	$(FC) $(FFLAGS) -c $(join $(SRCDIR),/fhash_modules.f90)

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) $(join $(OBJDIR),/*.o) $(join $(OBJDIR),/*.mod)

cleaner:
	$(RM) $(join $(OBJDIR),/*.o) $(join $(OBJDIR),/*.mod) $(join $(RUNDIR),/*.exe) $(join $(RUNDIR),/*.out.txt) $(join $(RUNDIR),/*.bin) $(join $(RUNDIR),/fort.*)

test:
	./TEST_INTEGRITY/test_integrity.sh
