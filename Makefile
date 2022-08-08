all : nvortex2d.bin

nvortex2d.bin : nvortex2d.cpp
	g++ -march=native -O2 -o nvortex2d.bin nvortex2d.cpp

clean :
	rm -f nvortex2d.bin
