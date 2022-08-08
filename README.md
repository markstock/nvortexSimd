# nvortexSimd
Simple N-body codes using C++ experimental simd support

## Compile and run
On a Linux workstation with GCC 10+, this should be as easy as

	make
	./nvortex2d.bin

## To do
We should be able to iterate over a `std::vector<float>` but cast it to a simd type just before the kernel.

## About
This will be a copy of [nvortexVc](https://github.com/Applied-Scientific-Research/nvortexVc) but using the C++-built-in experimental simd extensions.
Learn more about `experimental/simd` at [cppreference](https://en.cppreference.com/w/cpp/experimental/simd) and [stackoverflow](https://stackoverflow.com/questions/58584012/how-to-use-stdexperimentalsimd).

Note that GCC will auto-vectorize the 2D code, but not the 3D code (too many registers?). Unlike Vc, there is no vectorized square root yet, so the 3D code is much slower than the 2D code now.
