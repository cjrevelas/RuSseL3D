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
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_SERIAL
else
MPI_OPTIONS = -DUSE_MPI
#topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR
topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR_3
#topdir = /home/asgouros/TrioStountzes/MUMPS/MUMPS_5.2.1_PAR_EXPERIMENTAL
endif

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc                                #defines PLAT, LIBEXT, LIBMUMPS_COMMON, FC
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

# Flags of production and debug runs
FC_PROD=-O3
FC_DEBUG=-O0 -g -fcheck=all -Wall
#FC_DEBUG=-O0 -g -fcheck=all -Wall -Wextra -Wno-unused-dummy-argument

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
MODDIR=$(SRCDIR)/mod
RUNDIR=run

MODULES=$(OBJDIR)/parser_vars_mod.o\
	$(OBJDIR)/constants_mod.o\
	$(OBJDIR)/kcw_mod.o\
	$(OBJDIR)/fhash_mod.o\
	$(OBJDIR)/mpistuff_mod.o\
	$(OBJDIR)/error_handing_mod.o\
	$(OBJDIR)/write_helper_mod.o\

OBJECTS=$(OBJDIR)/matrix_assemble.o\
	$(OBJDIR)/part_fun.o\
	$(OBJDIR)/grafted_chains.o\
	$(OBJDIR)/periodic_dumper.o\
	$(OBJDIR)/export_field.o\
	$(OBJDIR)/parser.o\
	$(OBJDIR)/grafted_init_cond.o\
	$(OBJDIR)/init_field.o\
	$(OBJDIR)/calc_scf_params.o\
	$(OBJDIR)/init_time.o\
	$(OBJDIR)/simpsonkoef.o\
	$(OBJDIR)/quadinterp_koef.o\
	$(OBJDIR)/spat_3d.o\
	$(OBJDIR)/surf_pot.o\
	$(OBJDIR)/tetshp.o\
	$(OBJDIR)/qprint.o\
	$(OBJDIR)/mesh_io_3d.o\
	$(OBJDIR)/gauss_3d.o\
	$(OBJDIR)/edwards.o\
	$(OBJDIR)/adh_ten.o\
	$(OBJDIR)/main.o\
	$(OBJDIR)/mumps_sub.o\
	$(OBJDIR)/get_sys_time.o\
	$(OBJDIR)/dirichlet.o\
	$(OBJDIR)/convolution.o

$(OBJDIR)/%.o : $(MODDIR)/%.f90 #$(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(LIBFS) -J$(OBJDIR) -c -o $@ $?

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(LIBFS) -J$(OBJDIR) -c -o $@ $?

EXEC=$(RUNDIR)/fem_3d.exe

$(EXEC):$(LIBDMUMPS) $(MODULES) $(OBJECTS)
	$(FL) -o $(EXEC) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

$(OBJDIR)/mumps_sub.o:
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/$(SCRDIR) -I$(OBJDIR) -cpp $(CPPFLAGS) -c $(SRCDIR)/mumps_sub.f90 $(OUTF)$*.o

$(OBJDIR)/fhash_modules : $(SRCDIR)/fhash.f90 $(SRCDIR)/fhash_modules.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/fhash_modules.f90

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod

cleaner:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/*.exe $(RUNDIR)/*.out.txt $(RUNDIR)/*.bin $(RUNDIR)/fort.*

test:
	./test_integrity/test_integrity.sh
