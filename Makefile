# define the compiler manually, because Cray
CXX=g++
ifneq (, $(shell which CC 2>/dev/null))
  CXX=CC
endif

CXXFLAGS=-O3 -march=native -std=c++2b -fopenmp

all: nvortex2d.bin nvortex3d.bin

%.bin: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f *.bin
