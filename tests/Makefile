include ../make.inc

# reset flags to debug flags
CXXFLAGS = -g -O0 -Wall -std=c++11 -qopenmp

.PHONY: default clean test serial_test para4_test

serial_objects = test_atoms.exe test_matrix.exe test_timefreq.exe test_pbc.exe \
		  test_blas_lapack.exe test_complexmatrix.exe test_meanfield.exe
para_objects = test_parallel_mpi.exe

test_targets = $(patsubst %.exe, run_%, $(serial_objects) $(para_objects))

test_timefreq_deps = timefreq,envs,complexmatrix,matrix,mathtools,constants
test_parallel_mpi_deps = parallel_mpi,matrix,complexmatrix,vector3_order
test_pbc_deps = pbc,vector3_order
test_meanfield_deps = meanfield,matrix,complexmatrix,constants

default: $(serial_objects) $(para_objects)

test: $(test_targets)

test_atoms.exe: test_atoms.o
	$(CXX) $(CXXFLAGS) -o $@ $^
run_test_atoms: test_atoms.exe
	./$<

test_pbc.exe: test_pbc.o $(shell echo ../src/{$(test_pbc_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^
run_test_pbc: test_pbc.exe
	./$<

test_timefreq.exe: test_timefreq.o $(shell echo ../src/{$(test_timefreq_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_timefreq: test_timefreq.exe
	./$<

test_matrix.exe: test_matrix.o ../src/matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_matrix: test_matrix.exe
	./$<

test_complexmatrix.exe: test_complexmatrix.o ../src/complexmatrix.o ../src/matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_complexmatrix: test_complexmatrix.exe
	./$<

test_parallel_mpi.exe: test_parallel_mpi.o $(shell echo ../src/{$(test_parallel_mpi_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_parallel_mpi: test_parallel_mpi.exe
	mpirun -np 4 ./$<

test_blas_lapack.exe: test_blas_lapack.o ../src/complexmatrix.o ../src/matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_blas_lapack: test_blas_lapack.exe
	./$<

test_meanfield.exe: test_meanfield.o $(shell echo ../src/{$(test_meanfield_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)
run_test_meanfield: test_meanfield.exe
	./$<

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../src -c $< -o $@

../src/envs.o: ../src/envs.cpp
	$(CXX) $(CXXFLAGS) -I../src -c $(PWD)/$< -o $@

clean:
	rm -f *.o *.exe
