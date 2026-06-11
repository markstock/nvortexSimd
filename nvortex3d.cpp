/*
 * nvortexSimd - test platform for SIMD-acceleration of an N-vortex solver
 *
 * Copyright (c) 2022,6 Mark J. Stock <markjstock@gmail.com>
 */

#include <array>
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
  constexpr std::size_t JBLOCK = 64;   // source SIMD vectors per cache block; tune for L1 size
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

  // create temporary r squared array
  std::vector<native_simd<FLOAT>> r2(nvec);
  for (size_t i=0; i<r.size(); ++i) {
    r2[i] = r[i] * r[i];
  }

  // do a smaller precursor calculation to warm things up
  std::vector<size_t> itargs = {{x.size()/2, x.size(), x.size()}};

  for (size_t thissize : itargs) {
  std::cout << "run with " << VECREG_SIZE*thissize << " target parts" << std::endl;
  auto start = std::chrono::steady_clock::now();

  // loop over targets one i-block at a time; j blocked so each source vector is
  // loaded once and reused across all target lanes before moving to the next block
  #pragma omp parallel for schedule(static)
  for (size_t i=0; i<thissize; ++i) {
    const size_t iilim = std::min(VECREG_SIZE, ntarg-i*VECREG_SIZE);

    // Per-lane target positions and velocity accumulators
    std::array<native_simd<FLOAT>, VECREG_SIZE> usum, vsum, wsum;
    std::array<native_simd<FLOAT>, VECREG_SIZE> tx, ty, tz, tr2;
    for (size_t ii=0; ii<iilim; ++ii) {
      tx[ii]   = x[i][ii];
      ty[ii]   = y[i][ii];
      tz[ii]   = z[i][ii];
      tr2[ii] = r2[i][ii];
      usum[ii] = 0.f;
      vsum[ii] = 0.f;
      wsum[ii] = 0.f;
    }

    // Blocked source loop: source data for JBLOCK vectors stays in L1 cache
    // while all target lanes in this i-block consume it
    for (size_t jb=0; jb<x.size(); jb+=JBLOCK) {
      const size_t jend = std::min(jb+JBLOCK, x.size());
      for (size_t j=jb; j<jend; ++j) {

        // Load source data once; reused across all iilim target lanes below
        const native_simd<FLOAT> xj  = x[j];
        const native_simd<FLOAT> yj  = y[j];
        const native_simd<FLOAT> zj  = z[j];
        const native_simd<FLOAT> sxj = sx[j];
        const native_simd<FLOAT> syj = sy[j];
        const native_simd<FLOAT> szj = sz[j];
        const native_simd<FLOAT> r2j = r2[j];

        for (size_t ii=0; ii<iilim; ++ii) {
          const native_simd<FLOAT> dx   = tx[ii]  - xj;
          const native_simd<FLOAT> dy   = ty[ii]  - yj;
          const native_simd<FLOAT> dz   = tz[ii]  - zj;
          const native_simd<FLOAT> dist = dx*dx + dy*dy + dz*dz + r2j + tr2[ii];

          // correct kernel, with sqrt
          const native_simd<FLOAT> invDist = 1.f / sqrt(dist);
          const native_simd<FLOAT> fac     = invDist / dist;

          usum[ii] += fac * (dy*szj - dz*syj);
          vsum[ii] += fac * (dz*sxj - dx*szj);
          wsum[ii] += fac * (dx*syj - dy*sxj);
        }
      }
    }

    // Reduce across SIMD lanes and store
    for (size_t ii=0; ii<iilim; ++ii) {
      u[i][ii] = std::experimental::reduce(usum[ii]);
      v[i][ii] = std::experimental::reduce(vsum[ii]);
      w[i][ii] = std::experimental::reduce(wsum[ii]);
      if (not timeit and i*VECREG_SIZE+ii < 10) std::cout << "vel at " << x[i][ii] << " " << y[i][ii] << " " << z[i][ii] << " is " << u[i][ii] << " " << v[i][ii] << " " << w[i][ii] << "\n";
    }
  }

  auto end = std::chrono::steady_clock::now();
  const float sec = std::chrono::duration<double>(end - start).count();

  const float flops = (thissize*VECREG_SIZE)*(1+nvec*VECREG_SIZE*29);
  std::cout << "  performed " << flops << " flops\n";

  if (timeit) {
    std::cout << "  work complete in " << sec << " sec\n";
    std::cout << "  performance is " << (flops / sec) * 1e-9 << " GF/s\n";
  }

  } // end loop over runs

  return 0;
}
