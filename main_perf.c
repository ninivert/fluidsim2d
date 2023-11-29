// OpenMP: make USE_OMP=1 main_perf
// MPI: no

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "solver.h"
#include "utils.h"
#include "main_utils.h"
#include "spinny.h"

int main(int argc, char** argv) {
#ifdef _OPENMP
  #pragma omp single
  printf("[omp] max threads=%d\n", omp_get_max_threads());
#else
  printf("[no omp]\n");
#endif

  const uint nsteps = 1000;
  uint nx = 100, ny = 100;

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 100);
  sim_print(sim);

  long long t0 = current_timestamp();

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("[sim] iter %u/%u\n", istep, nsteps);
    spinny(sim);  // apply sources
    sim_step(sim);  // solve
  }

  printf("[sim] finished in %llu milliseconds\n", current_timestamp() - t0);
  write_results(sim);
  sim_free(sim);

  return 0;
}