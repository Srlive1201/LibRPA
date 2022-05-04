#Makefile
include ../make.inc

# Project information
VERSION = 1.0.0
EXENAME = chi0_main.exe
# TODO: test the library
LIBNAME = librpa.a

LIBS = $(LAPACK_DIR) $(LAPACK_LIB)

object_main = lib_main.o
objects = \
		  aperiodic_chi0.o \
		  cal_periodic_chi0.o \
		  complexmatrix.o \
		  constants.o \
		  envs.o \
		  Gauss_Quadrature.o \
		  input.o \
		  matrix.o \
		  matrix3.o \
		  meanfield.o \
		  parallel_mpi.o \
		  profiler.o \
		  read_aims.o \
		  ri.o \
		  timefreq.o \
		  vector3_order.o

.PHONY: default clean lib exe

default: exe

exe: $(EXENAME)
	cp -a $(EXENAME) ../

lib: $(LIBNAME)
	cp -a $(LIBNAME) ../

clean:
	rm -rf *.o *.gch *.exe *.a

$(EXENAME): $(object_main) $(objects)
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LIBS)

$(LIBNAME): $(objects)
	$(AR) $(ARFLAGS) $(LIBNAME) $(objects) 
	$(AR) s $(LIBNAME)

.cpp.o:
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

envs.o: envs.cpp
	$(MPICXX) $(CXXFLAGS) -c $(PWD)/$< -o $@
