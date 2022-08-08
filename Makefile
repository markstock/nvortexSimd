all : nvortex2d.bin nvortex3d.bin

%.bin : %.cpp
	g++ -march=native -O2 -o $@ $<

clean :
	rm -f nvortex2d.bin nvortex3d.bin
