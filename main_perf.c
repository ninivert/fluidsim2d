#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "omp.h"

#include "solver.c"
#include "utils.h"

long long current_timestamp() {
  struct timeval te; 
  gettimeofday(&te, NULL); // get current time
  long long milliseconds = te.tv_sec*1000LL + te.tv_usec/1000; // calculate milliseconds
  // printf("milliseconds: %lld\n", milliseconds);
  return milliseconds;
}

int write_results(Sim* sim) {
  FILE *fp = fopen("sim.out", "w");
  if (!fp) {
    fprintf(stderr, "unable to open file");
    return 1;
  }
  fprintf(fp, "simulation with nx=%ui, ny=%ui, dt=%f, diff=%f, visc=%f\n", sim->nx, sim->ny, sim->dt, sim->diff, sim->visc);
  fprintf(fp, "rho\n");
  writearr(fp, sim->rho, sim->nx*sim->ny);
  fprintf(fp, "vx\n");
  writearr(fp, sim->vx, sim->nx*sim->ny);
  fprintf(fp, "vy\n");
  writearr(fp, sim->vy, sim->nx*sim->ny);
  fclose(fp);
}

int main(int argc, char** argv) {
  const nsteps = 1000;
  uint nx = 100, ny = 100;

  if (argc == 3) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  } else {
    printf("usage: ./main_perf nx ny\n");
    printf("using default values.\n");
  }

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 20);
  sim_print(sim);

  long long t0 = current_timestamp();

  printf("max threads:%d\n", omp_get_max_threads());

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("%u / %u\n", istep, nsteps);

    // spinny demo
    for (size_t y = 0; y < sim->ny; ++y) { for (size_t x = 0; x < sim->nx; ++x) {
      double xf = linterp(x, 0, nx-1, -1, 1), yf = linterp(y, 0, ny-1, -1, 1);
      double r = sqrt(xf*xf + yf*yf), theta = atan2(yf, xf);
      sim->vx_prev[IX(x, y)] = 0.01*r*sin(theta);
      sim->vy_prev[IX(x, y)] = -0.01*r*cos(theta);
    }}
    for (size_t i = 0; i < nx*ny; ++i) sim->rho_prev[i] = 0.0;
    sim->rho_prev[IX(nx/2+10,ny/2+10)] = 10.0;
    
    sim_step(sim);
  }

  printf("finished in %llu milliseconds\n", current_timestamp() - t0);

  write_results(sim);

  sim_free(sim);

  return 0;
}