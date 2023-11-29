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
  uint max_threads = 1;
  if (argc == 4) { max_threads = atoi(argv[3]); }
  omp_set_num_threads(max_threads);
  printf("[omp] using %d threads\n", omp_get_max_threads());
#else
  printf("[no omp]\n");
#endif

  const uint nsteps = 1000;
  uint nx = 100, ny = 100;
  if (argc >= 3) { nx = atoi(argv[1]); ny = atoi(argv[2]); }

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 100);
  sim_print(sim);

  long long t0 = current_timestamp();
  long long totaltime = 0;

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("[sim] iter %u/%u\n", istep, nsteps);
    long long tstep0 = current_timestamp();
    spinny(sim);  // apply sources
    sim_step(sim);  // solve
    totaltime += current_timestamp() - tstep0;
  }

  printf("[sim] finished in %llu milliseconds, avg local compute time per step = %f ms\n", current_timestamp() - t0, (double) totaltime / (double) nsteps);
  write_results(sim);
  sim_free(sim);

  return 0;
}