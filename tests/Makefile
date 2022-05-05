include ../make.inc

.PHONY: default clean

objects =

default: $(objects)

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../src -c $< -o $@

clean:
	rm -f *.o *.exe
