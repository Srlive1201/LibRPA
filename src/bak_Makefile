#Makefile
CC = mpiicpc
LAPACK_DIR = /opt/intel/mkl/lib/intel64

global = cal_periodic_chi0.o aperiodic_chi0.o Gauss_Quadrature.o vector3_order.o matrix.o complexmatrix.o matrix3.o input.o
object = lib_main.o ${global}

chi0_main.exe: ${object}
	${CC} -o chi0_main.exe ${object} -L${LAPACK_DIR} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -std=c++11

.cpp.o:
	${CC} -c -std=c++11 $< -o $@

.PHONY: clean
clean:
	rm -rf *.o *.gch
