#Makefile
include ../make.inc

# Project information
VERSION = 1.0.0
EXENAME = chi0_main.exe
# TODO: test the library
LIBNAME = librpa.a

LIBS = $(LAPACK_LIBRARIES)

object_main = lib_main.o
objects = \
		  aperiodic_chi0.o \
		  cal_periodic_chi0.o \
		  chi0.o \
		  complexmatrix.o \
		  constants.o \
		  envs.o \
		  Gauss_Quadrature.o \
		  input.o \
		  mathtools.o \
		  matrix.o \
		  matrix3.o \
		  meanfield.o \
		  parallel_mpi.o \
		  pbc.o \
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
	rm -f *.o *.gch $(EXENAME) $(LIBNAME)
	rm -f ../$(EXENAME) ../$(LIBNAME)

$(EXENAME): $(object_main) $(objects)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LIBS)

$(LIBNAME): $(objects)
	$(AR) $(ARFLAGS) $@ $^
	$(AR) s $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

envs.o: envs.cpp
	$(CXX) $(CXXFLAGS) -c $(PWD)/$< -o $@
