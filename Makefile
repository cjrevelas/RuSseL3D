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
	$(OBJDIR)/interp.o\
	$(OBJDIR)/arrays.o\
	$(OBJDIR)/iofiles.o\
	$(OBJDIR)/error_handing.o\
	$(OBJDIR)/force_fields.o\
	$(OBJDIR)/delta.o\
	$(OBJDIR)/hist.o\

OBJECTS=$(OBJDIR)/matrix_assemble.o\
	$(OBJDIR)/get_part_func.o\
	$(OBJDIR)/get_nchains.o\
	$(OBJDIR)/get_volnp.o\
	$(OBJDIR)/get_interp_quad_coeff.o\
	$(OBJDIR)/get_sys_time.o\
        $(OBJDIR)/get_contour_step.o\
	$(OBJDIR)/parser.o\
	$(OBJDIR)/init_field.o\
	$(OBJDIR)/init_arrays.o\
	$(OBJDIR)/init_delta.o\
	$(OBJDIR)/init_scf_params.o\
	$(OBJDIR)/init_chain_contour.o\
	$(OBJDIR)/alloc_hist.o\
	$(OBJDIR)/interp_fem.o\
	$(OBJDIR)/matinv4.o\
	$(OBJDIR)/is_node_inside_el.o\
	$(OBJDIR)/spat_3d.o\
	$(OBJDIR)/tetshp.o\
	$(OBJDIR)/compute_delta_numer.o\
	$(OBJDIR)/import_mesh.o\
	$(OBJDIR)/gausspoints.o\
	$(OBJDIR)/edwards.o\
	$(OBJDIR)/energies.o\
	$(OBJDIR)/mumps_sub.o\
	$(OBJDIR)/dirichlet.o\
	$(OBJDIR)/convolution.o\
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
        $(OBJDIR)/compute_gradient.o\
        $(OBJDIR)/compute_phi_end_middle_nodal.o\
	$(OBJDIR)/main.o\

EXEC=$(RUNDIR)/fem3d.exe
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

cleantest:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/*.exe $(RUNDIR)/*.out.txt $(RUNDIR)/*out.bin $(RUNDIR)/fort.* $(RUNDIR)/field.in.bin $(RUNDIR)/LOG*
###################################################################################################################################################
