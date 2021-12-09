###################################################################################################################################################
#                                                            MAKEFILE OF SCF-FEM CODE                                                             #
###################################################################################################################################################
#--------------------------------------------------SET COMPILER FLAGS AND LIBRARIES TO BE LINKED--------------------------------------------------#
MAKE_MPI_RUN        = 0
MAKE_PRODUCTION_RUN = 0
PROD_OPTIONS        =
DEBUG_OPTIONS       = -DDEBUG_OUTPUTS #-DPRINT_AFULL
BOTH_OPTIONS        = -ffree-line-length-512# -DMUMPS_REPORT

ifeq ($(MAKE_MPI_RUN),0)
MPI_OPTIONS =
topdir      = /home/cjrevelas/programs/mumps/mumps_serial
else
MPI_OPTIONS = -DUSE_MPI
topdir      = /home/cjrevelas/programs/mumps/mumps_par
endif

libdir = $(topdir)/lib

.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
LIBDMUMPS       = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

FC_PROD  = -O3
FC_DEBUG = -O0 -g -fcheck=all -Wall -Wextra \
           -m64 -g -O0 -pedantic-errors -frepack-arrays -fdump-core -fbounds-check\
           -fimplicit-none -fbacktrace -ffree-line-length-none -frange-check\
           -Wall -Waliasing -Wampersand\
           -Wsurprising -Wunderflow -W #-Wno-unused-dummy-argument

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
        $(OBJDIR)/flags.o\
        $(OBJDIR)/eos.o\
	$(OBJDIR)/kcw.o\
	$(OBJDIR)/geometry.o\
	$(OBJDIR)/fhash_mod.o\
	$(OBJDIR)/mpistuff.o\
	$(OBJDIR)/write_helper.o\
	$(OBJDIR)/arrays.o\
	$(OBJDIR)/iofiles.o\
	$(OBJDIR)/error_handing.o\
	$(OBJDIR)/force_fields.o\
	$(OBJDIR)/delta.o\
	$(OBJDIR)/hist.o\

OBJECTS=$(OBJDIR)/parser.o\
	$(OBJDIR)/import_mesh.o\
	$(OBJDIR)/tools_histogram.o\
	$(OBJDIR)/tools_matinv4.o\
	$(OBJDIR)/solver_mumps.o\
	$(OBJDIR)/solver_dirichlet.o\
	$(OBJDIR)/solver_edwards.o\
	$(OBJDIR)/contour_convolution.o\
        $(OBJDIR)/contour_interp.o\
	$(OBJDIR)/is_node_inside_el.o\
	$(OBJDIR)/get_node_volume.o\
	$(OBJDIR)/get_contour_coeffs.o\
	$(OBJDIR)/get_sys_time.o\
        $(OBJDIR)/get_contour_step.o\
        $(OBJDIR)/get_gradient.o\
	$(OBJDIR)/get_delta_numer.o\
	$(OBJDIR)/init_field.o\
	$(OBJDIR)/init_arrays.o\
	$(OBJDIR)/init_delta.o\
	$(OBJDIR)/init_scf_params.o\
	$(OBJDIR)/init_chain_contour.o\
	$(OBJDIR)/fem_gausspoints.o\
	$(OBJDIR)/fem_integration.o\
        $(OBJDIR)/fem_matrix_assemble.o\
	$(OBJDIR)/fem_interpolation.o\
	$(OBJDIR)/fem_tetshpfun.o\
	$(OBJDIR)/export_computes.o\
	$(OBJDIR)/export_phi_indiv.o\
	$(OBJDIR)/export_q.o\
	$(OBJDIR)/export_delta.o\
	$(OBJDIR)/export_field_bin.o\
	$(OBJDIR)/export_phi_nodal.o\
        $(OBJDIR)/export_phi_smear.o\
	$(OBJDIR)/export_field.o\
	$(OBJDIR)/export_brush.o\
	$(OBJDIR)/export_brush99.o\
	$(OBJDIR)/export_chains_area.o\
        $(OBJDIR)/compute_phi_indiv.o\
        $(OBJDIR)/compute_phi_end_middle_nodal.o\
	$(OBJDIR)/compute_part_func_mx.o\
	$(OBJDIR)/compute_number_of_chains.o\
	$(OBJDIR)/compute_energies.o\
        $(OBJDIR)/compute_stretching_energy.o\
	$(OBJDIR)/main.o

EXEC=$(RUNDIR)/RuSseL
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

$(OBJDIR)/solver_mumps.o:
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/$(SCRDIR) -I$(OBJDIR) -cpp $(CPPFLAGS) -c $(SRCDIR)/solver_mumps.f90 $(OUTF)$*.o

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p
#-------------------------------------------------------------CLEAN WORKING DIRECTORY-------------------------------------------------------------#
clean:
	$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod

cleaner:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/RuSseL* $(RUNDIR)/*.out.txt $(RUNDIR)/*out.bin $(RUNDIR)/fort.*

#-----------------------------------------------------------RUN TESTS TO VERIFY CHANGES-----------------------------------------------------------#
test:
	./test_integrity/test_integrity.sh

cleantest:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/RuSseL* $(RUNDIR)/*.out.txt $(RUNDIR)/*out.bin $(RUNDIR)/fort.* $(RUNDIR)/field.in.bin $(RUNDIR)/LOG*
###################################################################################################################################################
