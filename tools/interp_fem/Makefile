#######################################################################
#                  Makefile of FEM-INTERPOLATION code                 #
#######################################################################
MAKE_PRODUCTION_RUN = 1
FC                  = /usr/bin/gfortran
PROD_OPTIONS        =
DEBUG_OPTIONS       = -DDEBUG_OUTPUTS
BOTH_OPTIONS        =

.SECONDEXPANSION:
FC_PROD  = -O3
#FC_DEBUG = -O0 -g -fcheck=all # -Wall -Wextra -Wno-unused-dummy-argument


FC_DEBUG = -m64 -g -O0 -pedantic-errors -frepack-arrays -fdump-core -fbounds-check\
           -fimplicit-none -fbacktrace -ffree-line-length-none -frange-check\
           -Wall -Waliasing -Wampersand\
           -Wsurprising -Wunderflow -W

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

MODULES=

OBJECTS=$(OBJDIR)/isInside.o\
	$(OBJDIR)/matinv4.o\
	$(OBJDIR)/interp.o\
	$(OBJDIR)/tetshp.o\
	$(OBJDIR)/gausspoints.o\
	$(OBJDIR)/main.o

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(LIBFS) -J$(OBJDIR) -c -o $@ $?

EXEC=$(RUNDIR)/interp.exe

$(EXEC):$(MODULES) $(OBJECTS)
	$(FC) -o $(EXEC) $(MODULES) $(OBJECTS)

.SUFFIXES: .F .f90 .h .p

clean:
	$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod

cleaner:
	$(RM)  $(OBJDIR)/*.o $(OBJDIR)/*.mod $(RUNDIR)/*.exe $(RUNDIR)/*.out.txt $(RUNDIR)/*.log
