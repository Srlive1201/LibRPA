include ../make.inc

.PHONY: default clean

objects = test_atoms.exe test_matrix.exe test_timefreq.exe test_parallel_mpi.exe test_pbc.exe

test_matrix_deps = timefreq,envs,matrix,complexmatrix,mathtools,constants
test_parallel_mpi_deps = parallel_mpi,matrix,complexmatrix,vector3_order
test_pbc_deps = pbc,vector3_order

default: $(objects)

test_atoms.exe: test_atoms.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_pbc.exe: test_pbc.o $(shell echo ../src/{$(test_pbc_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

test_timefreq.exe: test_timefreq.o $(shell echo ../src/{$(test_matrix_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)

test_matrix.exe: test_matrix.o ../src/matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)

test_parallel_mpi.exe: test_parallel_mpi.o $(shell echo ../src/{$(test_parallel_mpi_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../src -c $< -o $@

../src/envs.o: ../src/envs.cpp
	$(CXX) $(CXXFLAGS) -I../src -c $(PWD)/$< -o $@

clean:
	rm -f *.o *.exe
