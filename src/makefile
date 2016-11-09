# Makefile for Main program
FC= # mpif90/gfortran/...

# MUMPS Library
LIBMUMPSTOPDIR = # PROVIDE THE PATH TO MUMPS LIBRARY
# (1) ORDERINGS
LMETISDIR = /usr/local/lib
LMETIS    = -L$(LMETISDIR) -lmetis
LPORDDIR = $(LIBMUMPSTOPDIR)/PORD/lib/
LPORD    = -L$(LPORDDIR) -lpord
LSCOTCHDIR = /usr/local/lib
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

# (2) LIBS
SCALAP  = -lscalapack-openmpi -lblacs-openmpi  -lblacsF77init-openmpi -lblacsCinit-openmpi
LIBPAR = $(SCALAP)   -lmpi #-lmpi_f77
LIBS = $(LIBPAR)
#LIBBLAS = -lblas
#LIBBLAS = -L/usr/lib/openblas-base -lopenblas
LIBBLAS = -L/usr/lib/atlas-base/atlas -lblas
LIBOTHERS = -lpthread

# (3) MUMPS
LIBMUMPS_COMMON = $(LIBMUMPSTOPDIR)/lib/libmumps_common.a
LIBMUMPS_D = $(LIBMUMPSTOPDIR)/lib/libdmumps.a
LIBMUMPS_Z = $(LIBMUMPSTOPDIR)/lib/libzmumps.a
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
INCLUDEMUMPS = -I$(LIBMUMPSTOPDIR)/include
INCLUDEMPI = -I/usr/lib/openmpi/include

# COMPILATION FLAGS
#      -Wall -fbounds-check -g -std=f95 -ftrace=full
#FFLAGS= -O0 -Wfatal-errors -z muldefs -fbounds-check -ggdb -pg -g
FFLAGS= -Ofast -z muldefs

# COMPILATION
radau_main:	units.o utilities.o matrixrw.o radau_common.o radau_fx_seq.o radau_sa_seq.o radau_ad_seq.o radau_main.o
	$(FC) -o radau_main radau_main.o units.o utilities.o matrixrw.o radau_common.o radau_fx_seq.o radau_sa_seq.o radau_ad_seq.o radau_main.o  $(LIBS) $(FFLAGS) $(LIBMUMPS_Z) $(LIBMUMPS_D) $(LIBMUMPS_COMMON) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

# DEPENDENCIES
units.o:	units.f90
	$(FC) -c $(FFLAGS) units.f90 

utilities.o:	utilities.f90
	$(FC) -c $(FFLAGS) utilities.f90

matrixrw.o:	matrixrw.f90
	$(FC) -c $(FFLAGS) matrixrw.f90

radau_common.o:	radau_common.f90
	$(FC) -c $(FFLAGS) $(INCLUDEMUMPS) $(INCLUDEMPI) radau_common.f90

radau_fx_seq.o:	radau_fx_seq.f90
	$(FC) -c $(FFLAGS) $(INCLUDEMUMPS) $(INCLUDEMPI) -DALLOW_NON_INIT radau_fx_seq.f90

radau_sa_seq.o:	radau_sa_seq.f90
	$(FC) -c $(FFLAGS) $(INCLUDEMUMPS) $(INCLUDEMPI) -DALLOW_NON_INIT radau_sa_seq.f90

radau_ad_seq.o:	radau_ad_seq.f90
	$(FC) -c $(FFLAGS) $(INCLUDEMUMPS) $(INCLUDEMPI) -DALLOW_NON_INIT radau_ad_seq.f90

radau_main.o:	radau_main.f90
	$(FC) -c $(FFLAGS) $(INCLUDEMUMPS) $(INCLUDEMPI) radau_main.f90
	
clean:
	$(RM) *.o *.mod

veryclean:	clean
	$(RM) main
