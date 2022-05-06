include ../make.inc

.PHONY: default clean

objects = test_atoms.exe

default: $(objects)

test_atoms.exe: test_atoms.o
	$(CXX) $(CXXFLAGS) -o $@ $^

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../src -c $< -o $@

clean:
	rm -f *.o *.exe
