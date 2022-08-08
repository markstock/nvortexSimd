/*
 * nvortexSimd - test platform for SIMD-acceleration of an N-vortex solver
 *
 * Copyright (c) 2022 Mark J. Stock <markjstock@gmail.com>
*/

#include <iostream>
#include <random>
#include <chrono>
#include <experimental/simd>

using std::experimental::native_simd;

static void usage() {
  fprintf(stderr, "Usage: nvortex2d.bin [-n=<number>]\n");
  exit(1);
}

int main(int argc, char const *argv[]) {

  size_t ntarg = 10000;
  const bool timeit = true;

  if (argc > 1) {
    if (strncmp(argv[1], "-n=", 3) == 0) {
      int num = atof(argv[1] + 3);
      if (num < 1) usage();
      ntarg = num;
    }
  }

  std::cout << "experimental/simd test01\n";

  constexpr std::size_t VECREG_SIZE = native_simd<float>::size();
  std::cout << "  register size is " << VECREG_SIZE << " wide\n";

  std::cout << "  number of particles is " << ntarg << "\n";
  const size_t nvec = 1 + (ntarg-1)/VECREG_SIZE;
  std::vector<native_simd<float>> x(nvec),y(nvec),s(nvec),r(nvec);
  std::vector<native_simd<float>> u(nvec),v(nvec);
  std::cout << "  vector length is " << x.size() << " entries\n";

  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::mt19937 generator(12345);
  std::uniform_real_distribution<float> distribution(-1.f, 1.f);
  std::uniform_real_distribution<float> posdistro(0.1f, 1.f);

  // initialize to random
  for (size_t i=0; i<x.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    x[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<y.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    //y[i][j] = distribution(generator);
    y[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<s.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    //s[i][j] = distribution(generator);
    s[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<r.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    //r[i][j] = 0.1f * posdistro(generator);
    r[i][j] = (i*VECREG_SIZE+j < ntarg) ? 0.1f*posdistro(generator) : 1.f;
  }

  auto start = std::chrono::system_clock::now();

  // perform some math

  // loop over targets, one at a time
  for (size_t i=0; i<x.size(); ++i) {
  for (size_t ii=0; ii<std::min(VECREG_SIZE,ntarg-i*VECREG_SIZE); ++ii) {
    // same results as:
    //for (size_t ii=0; ii<VECREG_SIZE; ++ii) {

    native_simd<float> usum = 0.f;
    native_simd<float> vsum = 0.f;
    const native_simd<float> targrad = r[i][ii];
    const native_simd<float> tr2 = targrad*targrad;

    // loop over sources, vector at a time
    for (size_t j=0; j<x.size(); ++j) {

      const native_simd<float> dx = x[i][ii] - x[j];
      const native_simd<float> dy = y[i][ii] - y[j];
      const native_simd<float> dist = dx*dx + dy*dy + r[j]*r[j] + tr2;
      const native_simd<float> fac = s[j] / dist;

      usum -= dy * fac;
      vsum += dx * fac;
    }

    u[i][ii] = std::experimental::reduce(usum);
    v[i][ii] = std::experimental::reduce(vsum);
    if (not timeit and i*VECREG_SIZE+ii < 10) std::cout << "vel at " << x[i][ii] << " " << y[i][ii] << " is " << u[i][ii] << " " << v[i][ii] << "\n";
  }
  }

  const float flops = ntarg*(1+nvec*VECREG_SIZE*13);
  std::cout << "performed " << flops << " flops\n";

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;

  if (timeit) {
    std::cout << "work complete in " << (elapsed.count()/1.e+9) << " sec\n";
    std::cout << "performance is " << (flops / elapsed.count()) << " GF/s\n";
  }

  return 0;
}
