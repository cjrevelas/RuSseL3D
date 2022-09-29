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
topdir      = /home/cjrevelas/libraries/mumps_serial
else
MPI_OPTIONS = -DUSE_MPI
topdir      = /home/cjrevelas/libraries/mumps_serial
endif

libdir = $(topdir)/lib

.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
LIBDMUMPS       = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

FC_PROD  = -O3
FC_DEBUG = -O1 -g -fcheck=all -Wall -Wextra \
           -m64 -g -pedantic-errors -frepack-arrays -fdump-core -fbounds-check\
           -fimplicit-none -fbacktrace -ffree-line-length-none -frange-check\
           -Wall -Waliasing -Wampersand\
           -Wsurprising -Wunderflow -W -Wno-unused-dummy-argument

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

MODULES=$(OBJDIR)/parser_vars_mod.o\
	$(OBJDIR)/constants_mod.o\
        $(OBJDIR)/flags_mod.o\
	$(OBJDIR)/defaults_mod.o\
        $(OBJDIR)/eos_mod.o\
	$(OBJDIR)/kcw_mod.o\
	$(OBJDIR)/fhash_mod.o\
	$(OBJDIR)/geometry_mod.o\
	$(OBJDIR)/mpistuff_mod.o\
	$(OBJDIR)/write_helper_mod.o\
	$(OBJDIR)/arrays_mod.o\
	$(OBJDIR)/iofiles_mod.o\
	$(OBJDIR)/error_handing_mod.o\
	$(OBJDIR)/force_fields_mod.o\
	$(OBJDIR)/delta_mod.o\
	$(OBJDIR)/hist_mod.o\

OBJECTS=$(OBJDIR)/ComputeContourCoeffs.o\
        $(OBJDIR)/ComputeContourStep.o\
        $(OBJDIR)/ComputeDeltaNumerical.o\
        $(OBJDIR)/ComputeGradient.o\
	$(OBJDIR)/ComputeNodeVolume.o\
	$(OBJDIR)/ComputeNumberOfChains.o\
	$(OBJDIR)/ComputePartitionMatrix.o\
        $(OBJDIR)/ComputeIndivProfile.o\
        $(OBJDIR)/ComputeStretchingEnergy.o\
	$(OBJDIR)/ContourConvolution.o\
        $(OBJDIR)/ContourInterpolation.o\
        $(OBJDIR)/ExportAdsorbed.o\
	$(OBJDIR)/ExportBrush99.o\
	$(OBJDIR)/ExportBrush.o\
	$(OBJDIR)/ExportChainsArea.o\
	$(OBJDIR)/ExportComputes.o\
	$(OBJDIR)/ExportDelta.o\
	$(OBJDIR)/ExportEnergies.o\
	$(OBJDIR)/ExportFieldBinary.o\
	$(OBJDIR)/ExportFieldAscii.o\
        $(OBJDIR)/ExportEndMiddleProfile.o\
	$(OBJDIR)/ExportIndivProfile.o\
	$(OBJDIR)/ExportNodalProfile.o\
        $(OBJDIR)/ExportSmearedProfile.o\
	$(OBJDIR)/ExportPropagator.o\
	$(OBJDIR)/ExportVtuProfiles.o\
	$(OBJDIR)/ExportVtuIndiv.o\
	$(OBJDIR)/FemApplyPeriodicity.o\
	$(OBJDIR)/FemNonZeroEntries.o\
	$(OBJDIR)/FemPeriodicEdges.o\
	$(OBJDIR)/FemPeriodicCorners.o\
	$(OBJDIR)/FemGaussPoints.o\
	$(OBJDIR)/FemIntegration.o\
	$(OBJDIR)/FemInterpolation.o\
	$(OBJDIR)/FemIsPointInElement.o\
        $(OBJDIR)/FemMatrixAssemble.o\
	$(OBJDIR)/FemShapeFunctions.o\
	$(OBJDIR)/InitArrays.o\
	$(OBJDIR)/InitChainContour.o\
	$(OBJDIR)/InitDelta.o\
	$(OBJDIR)/InitField.o\
	$(OBJDIR)/InitScfParams.o\
	$(OBJDIR)/MeshAppendPeriodicPairs.o\
	$(OBJDIR)/MeshPeriodicNodePairs.o\
	$(OBJDIR)/MeshPeriodicDestNeighbors.o\
	$(OBJDIR)/MeshBulkNodePairs.o\
	$(OBJDIR)/MeshDirichletFaces.o\
	$(OBJDIR)/MeshElementsOfNode.o\
	$(OBJDIR)/MeshFaceEntities.o\
	$(OBJDIR)/MeshPeriodicFaces.o\
	$(OBJDIR)/MeshPeriodicEdges.o\
	$(OBJDIR)/MeshProfile.o\
        $(OBJDIR)/ParserInput.o\
	$(OBJDIR)/ParserMesh.o\
	$(OBJDIR)/SolverEdwards.o\
	$(OBJDIR)/SolverMumps.o\
	$(OBJDIR)/ToolsHistogram.o\
	$(OBJDIR)/ToolsMatrixInversion.o\
	$(OBJDIR)/ToolsSystemTime.o\
	$(OBJDIR)/Main.o

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

$(OBJDIR)/SolverMumps.o:
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/$(SCRDIR) -I$(OBJDIR) -cpp $(CPPFLAGS) -c $(SRCDIR)/SolverMumps.f90 $(OUTF)$*.o

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p
#-------------------------------------------------------------CLEAN WORKING DIRECTORY-------------------------------------------------------------#
clean:
	$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod

cleaner:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/RuSseL* $(RUNDIR)/o.* $(RUNDIR)/d.* $(RUNDIR)/e.* $(RUNDIR)/*out.bin $(RUNDIR)/fort.*

#-----------------------------------------------------------RUN TESTS TO VERIFY CHANGES-----------------------------------------------------------#
test:
	./test_integrity/test_integrity.sh

cleantest:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/RuSseL* $(RUNDIR)/o.* $(RUNDIR)/d.* $(RUNDIR)/e.* $(RUNDIR)/*out.bin $(RUNDIR)/fort.* $(RUNDIR)/in.field.bin $(RUNDIR)/LOG*
###################################################################################################################################################
