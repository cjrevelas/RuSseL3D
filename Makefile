#######################################################################
#                                                                     #          
#                  Makefile of SCF-FEM code                           #
#                                                                     #
#######################################################################
#                      MUMPS SECTION                                  #
#######################################################################
#Set the path of the MUMPS directory
topdir = /home/cjrevelas/Documents/mumps_5.2.1

libdir = $(topdir)/lib
.SECONDEXPANSION:
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

PROD=-O3 -cpp -Wno-unused-dummy-argument
DEBUG=-O0 -g -fcheck=all -Wall -cpp -Wno-unused-dummy-argument

FCFLAGS=$(PROD)
#FCFLAGS=$(DEBUG)

LIBFS=#-lstdc++ #-lm

MODULES =  mdata.o xdata.o kcw.o fhash_modules.o
OBJECTS =  matrix_assemble.o  part_fun_phi.o\
	   scfinout.o  simpsonkoef.o spat_3d.o surf_pot.o tetshp.o qprint.o \
           mesh_io_3d.o gauss_3d.o edwards_free_film_fem.o adh_ten_alt.o main.o mumps_sub.o


.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBFS)  $*.f90

CMD=bench_3d.exe
$(CMD):$(LIBDMUMPS) $(MODULES) $(OBJECTS) 
#main: $(LIBDMUMPS) $(MODULES) $(OBJECTS) $(MAIN_OBJECT)
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
