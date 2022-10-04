# define the compiler manually, because Cray
CXX=g++
ifneq (, $(shell which CC))
  CXX=CC
endif

all : nvortex2d.bin nvortex3d.bin

%.bin : %.cpp
	$(CXX) -march=native -O2 -o $@ $< -fopenmp

clean :
	rm -f nvortex2d.bin nvortex3d.bin
