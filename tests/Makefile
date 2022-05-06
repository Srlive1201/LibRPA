include ../make.inc

.PHONY: default clean

objects = test_atoms.exe test_matrix.exe test_timefreq.exe

test_matrix_deps = timefreq,envs,matrix,complexmatrix

default: $(objects)

test_atoms.exe: test_atoms.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_timefreq.exe: test_timefreq.o $(shell echo ../src/{$(test_matrix_deps)}.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)

test_matrix.exe: test_matrix.o ../src/matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LAPACK_LIBRARIES)

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../src -c $< -o $@

../src/envs.o: ../src/envs.cpp
	$(CXX) $(CXXFLAGS) -I../src -c $(PWD)/$< -o $@

clean:
	rm -f *.o *.exe
