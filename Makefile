#Makefile
CPP_MPI = mpiicpc
CPP=icpc

LAPACK_DIR = /opt/intel/mkl/lib/intel64
LAPACK_LIB = -L${LAPACK_DIR} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
OPTS = -std=c++11 -qopenmp
#OPTS_MPI= -cxx=${CPP}

global = cal_periodic_chi0.o aperiodic_chi0.o Gauss_Quadrature.o vector3_order.o matrix.o complexmatrix.o matrix3.o input.o parallel_mpi.o global_class.o
object = lib_main.o ${global}

chi0_main.exe: ${object}
	${CPP_MPI} -g -o chi0_main.exe ${object} ${LAPACK_LIB} ${OPTS}
.cpp.o:
	${CPP_MPI} ${OPTS} -c $< -o $@

.PHONY: clean
clean:
	rm -rf *.o *.gch *.exe
