#
#  Makefile of the FEM code
#
#######################################################################
# MUMPS SECTION                                                       #
#######################################################################
#tolis
# Set the path of the MUMPS directory
#topdir = /home/asgouros/FEM/MUMPS_5.2.1
topdir = /home/parallels/Documents/ARIS_MUMPS

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)
FCFLAGS=-O3
LIBFS=#-lstdc++ #-lm

MODULES =  mdata.o xdata.o kcw.o
OBJECTS =  interpolation.o  matrix_assemble.o  part_fun_phi.o\
	       scfinout.o  simpsonkoef.o spat_3d.o surf_pot.o tetshp.o tetshp10.o qprint.o \
           mesh_io_3d.o gauss_3d.o edwards_free_film_fem.o adh_ten_alt.o  main.o mumps_sub_tolis_2.o


.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBFS)  $*.f90

CMD=fem_3d.x
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS)
#main: $(LIBDMUMPS) $(MODULES) $(OBJECTS) $(MAIN_OBJECT)
	   $(FL) -o $(CMD) $(OPTL) $(MODULES) $(OBJECTS)  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

mumps_sub_tolis_2.o :
	$(FC) $(OPTF) $(INCS) -I. -I$(topdir)/include -I$(topdir)/src -c $*.f90 $(OUTF)$*.o

.SUFFIXES: (.SUFFIXES) .F .f90 .h .p

clean:
	$(RM) *.o *.mod *.x
