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

#define FLOAT float

static void usage() {
  fprintf(stderr, "Usage: nvortex3d.bin [-n=<number>]\n");
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

  std::cout << "experimental/simd nvortex3d\n";

  constexpr std::size_t VECREG_SIZE = native_simd<FLOAT>::size();
  std::cout << "  register size is " << VECREG_SIZE << " wide\n";

  std::cout << "  number of particles is " << ntarg << "\n";
  const size_t nvec = 1 + (ntarg-1)/VECREG_SIZE;
  std::vector<native_simd<FLOAT>> x(nvec),y(nvec),z(nvec);
  std::vector<native_simd<FLOAT>> sx(nvec),sy(nvec),sz(nvec),r(nvec);
  std::vector<native_simd<FLOAT>> u(nvec),v(nvec),w(nvec);
  std::cout << "  vector length is " << x.size() << " entries\n";

  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::mt19937 generator(12345);
  std::uniform_real_distribution<FLOAT> distribution(-1.f, 1.f);
  std::uniform_real_distribution<FLOAT> posdistro(0.1f, 1.f);

  const FLOAT thisstrmag = 1.0 / std::sqrt(ntarg);
  const FLOAT thisrad    = (2./3.) / std::pow(ntarg,1./3.);

  // initialize to random
  for (size_t i=0; i<x.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    x[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<y.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    y[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<z.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    z[i][j] = (i*VECREG_SIZE+j < ntarg) ? distribution(generator) : 0.f;
  }
  for (size_t i=0; i<sx.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    sx[i][j] = (i*VECREG_SIZE+j < ntarg) ? thisstrmag*distribution(generator) : 0.f;
  }
  for (size_t i=0; i<sy.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    sy[i][j] = (i*VECREG_SIZE+j < ntarg) ? thisstrmag*distribution(generator) : 0.f;
  }
  for (size_t i=0; i<sz.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    sz[i][j] = (i*VECREG_SIZE+j < ntarg) ? thisstrmag*distribution(generator) : 0.f;
  }
  for (size_t i=0; i<r.size(); ++i) for (size_t j=0; j<VECREG_SIZE; ++j) {
    r[i][j] = (i*VECREG_SIZE+j < ntarg) ? thisrad*posdistro(generator) : 1.f;
  }

  // do a smaller precursor calculation to warm things up
  std::vector<size_t> itargs = {{x.size()/2, x.size(), x.size()}};

  for (size_t thissize : itargs) {
  std::cout << "run with " << VECREG_SIZE*thissize << " parts" << std::endl;
  auto start = std::chrono::system_clock::now();

  // loop over targets, one at a time
  #pragma omp parallel for schedule(guided)
  for (size_t i=0; i<thissize; ++i) {
  for (size_t ii=0; ii<std::min(VECREG_SIZE,ntarg-i*VECREG_SIZE); ++ii) {
    // same results as:
    //for (size_t ii=0; ii<VECREG_SIZE; ++ii) {

    native_simd<FLOAT> usum = 0.f;
    native_simd<FLOAT> vsum = 0.f;
    native_simd<FLOAT> wsum = 0.f;
    const native_simd<FLOAT> targrad = r[i][ii];
    const native_simd<FLOAT> tr2 = targrad*targrad;

    // loop over sources, vector at a time
    for (size_t j=0; j<x.size(); ++j) {

      const native_simd<FLOAT> dx = x[i][ii] - x[j];
      const native_simd<FLOAT> dy = y[i][ii] - y[j];
      const native_simd<FLOAT> dz = z[i][ii] - z[j];
      const native_simd<FLOAT> dist = dx*dx + dy*dy + dz*dz + r[j]*r[j] + tr2;

      // correct kernel, with sqrt
      const native_simd<FLOAT> fac = 1.f / (dist*sqrt(dist));
      // must serialize these calcs
      //native_simd<FLOAT> fac;
      //for (size_t jj=0; jj<VECREG_SIZE; ++jj) {
      //  fac[jj] = 1.f / (dist[jj]*std::sqrt(dist[jj]));
      //}

      // wrong kernel, not a Green's function solution
      //const native_simd<FLOAT> fac = 1.f / (dist*dist);

      usum += fac * (dy*sz[j] - dz*sy[j]);
      vsum += fac * (dz*sx[j] - dx*sz[j]);
      wsum += fac * (dx*sy[j] - dy*sx[j]);
    }

    u[i][ii] = std::experimental::reduce(usum);
    v[i][ii] = std::experimental::reduce(vsum);
    w[i][ii] = std::experimental::reduce(wsum);
    if (not timeit and i*VECREG_SIZE+ii < 10) std::cout << "vel at " << x[i][ii] << " " << y[i][ii] << " " << z[i][ii] << " is " << u[i][ii] << " " << v[i][ii] << " " << w[i][ii] << "\n";
  }
  }

  const float flops = (thissize*VECREG_SIZE)*(1+nvec*VECREG_SIZE*29);
  std::cout << "  performed " << flops << " flops\n";

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;

  if (timeit) {
    std::cout << "  work complete in " << (elapsed.count()/1.e+9) << " sec\n";
    std::cout << "  performance is " << (flops / elapsed.count()) << " GF/s\n";
  }

  } // end loop over runs

  return 0;
}
