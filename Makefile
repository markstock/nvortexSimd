# define the compiler manually, because Cray
CC=g++
#CC=CC

all : nvortex2d.bin nvortex3d.bin

%.bin : %.cpp
	$(CC) -march=native -O2 -o $@ $< -fopenmp

clean :
	rm -f nvortex2d.bin nvortex3d.bin
