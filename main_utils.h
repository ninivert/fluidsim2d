#pragma once

#include "solver.h"

#include <stdlib.h>
#include <stdio.h>

static void writearr(FILE* fp, double buf[], uint n) {
  for (uint i = 0; i < n; ++i) {
    // todo: a better way to dump the array than just the string
    fprintf(fp, "%.17f ", buf[i]);
  }
  fprintf(fp, "\n");
}

static long long current_timestamp() {
  struct timeval te; 
  gettimeofday(&te, NULL); // get current time
  long long milliseconds = te.tv_sec*1000LL + te.tv_usec/1000; // calculate milliseconds
  // printf("milliseconds: %lld\n", milliseconds);
  return milliseconds;
}

static int write_results(Sim* sim) {
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
  return 0;
}