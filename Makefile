###################################################################################################################################################
#                                                            MAKEFILE OF SCF-FEM CODE                                                             #
###################################################################################################################################################
#--------------------------------------------------SET COMPILER FLAGS AND LIBRARIES TO BE LINKED--------------------------------------------------#
MAKE_MPI_RUN        = 0
MAKE_PRODUCTION_RUN = 0
PROD_OPTIONS        =
DEBUG_OPTIONS       = -DDEBUG_OUTPUTS #-DPRINT_AFULL
BOTH_OPTIONS        = -DVARIABLE_DS_SCHEME# -DMUMPS_REPORT

ifeq ($(MAKE_MPI_RUN),0)
MPI_OPTIONS =
topdir      = /home/cjrevelas/programs/mumps-5.2.1/mumps_serial/
else
MPI_OPTIONS = -DUSE_MPI
topdir      = /usr/share/mumps-5.2.1_par_opt
endif

libdir = $(topdir)/lib

.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
LIBDMUMPS       = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

FC_PROD  = -O3
FC_DEBUG = -O0 -g -fcheck=all -Wall # -Wextra -Wno-unused-dummy-argument

ifeq ($(MAKE_PRODUCTION_RUN),0)
CPPFLAGS = $(DEBUG_OPTIONS) $(BOTH_OPTIONS) $(MPI_OPTIONS)
FCFLAGS  = $(FC_DEBUG) -cpp $(CPPFLAGS)
else
CPPFLAGS = $(PROD_OPTIONS) $(BOTH_OPTIONS) $(MPI_OPTIONS)
FCFLAGS  = $(FC_PROD) -cpp $(CPPFLAGS)
endif

LIBFS=#-lstdc++ #-lm
#--------------------------------------------------------SET FILE PATHS AND EXECUTABLE NAME-------------------------------------------------------#
OBJDIR=obj
SRCDIR=src
MODDIR=$(SRCDIR)/mod
RUNDIR=run

MODULES=$(OBJDIR)/parser_vars.o\
	$(OBJDIR)/constants.o\
	$(OBJDIR)/kcw.o\
	$(OBJDIR)/geometry.o\
	$(OBJDIR)/fhash_mod.o\
	$(OBJDIR)/mpistuff.o\
	$(OBJDIR)/write_helper.o\
	$(OBJDIR)/interp.o\
	$(OBJDIR)/iofiles.o\
	$(OBJDIR)/error_handing.o\
	$(OBJDIR)/force_fields.o\
	$(OBJDIR)/delta.o\

OBJECTS=$(OBJDIR)/matrix_assemble.o\
	$(OBJDIR)/part_fun.o\
	$(OBJDIR)/grafted_chains.o\
	$(OBJDIR)/periodic_dumper.o\
	$(OBJDIR)/export_field.o\
	$(OBJDIR)/parser.o\
	$(OBJDIR)/init_field.o\
	$(OBJDIR)/init_arrays.o\
	$(OBJDIR)/calc_scf_params.o\
	$(OBJDIR)/init_time.o\
	$(OBJDIR)/simpsonkoef.o\
	$(OBJDIR)/quadinterp_koef.o\
	$(OBJDIR)/spat_3d.o\
	$(OBJDIR)/tetshp.o\
	$(OBJDIR)/find_delta.o\
	$(OBJDIR)/qprint.o\
	$(OBJDIR)/mesh.o\
	$(OBJDIR)/gausspoints.o\
	$(OBJDIR)/edwards.o\
	$(OBJDIR)/energies.o\
	$(OBJDIR)/main.o\
	$(OBJDIR)/mumps_sub.o\
	$(OBJDIR)/get_sys_time.o\
	$(OBJDIR)/dirichlet.o\
	$(OBJDIR)/convolution.o

EXEC=$(RUNDIR)/fem_3d.exe
#-----------------------------------------------------------------COMPILE AND LINK----------------------------------------------------------------#
#modules compilation
$(OBJDIR)/%.o : $(MODDIR)/%.f90
	$(FC) $(FCFLAGS) $(LIBFS) -J$(OBJDIR) -c -o $@ $?

#routines compilation
$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(LIBFS) -J$(OBJDIR) -c -o $@ $?

#link object files to create the executable code
$(EXEC):$(LIBDMUMPS) $(MODULES) $(OBJECTS)
	$(FL) -o $(EXEC) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

$(OBJDIR)/mumps_sub.o:
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/$(SCRDIR) -I$(OBJDIR) -cpp $(CPPFLAGS) -c $(SRCDIR)/mumps_sub.f90 $(OUTF)$*.o

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p
#-------------------------------------------------------------CLEAN WORKING DIRECTORY-------------------------------------------------------------#
clean:
	$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod

cleaner:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/*.exe $(RUNDIR)/*.out.txt $(RUNDIR)/*out.bin $(RUNDIR)/fort.*
#-----------------------------------------------------------RUN TESTS TO VERIFY CHANGES-----------------------------------------------------------#
test:
	./test_integrity/test_integrity.sh
###################################################################################################################################################
