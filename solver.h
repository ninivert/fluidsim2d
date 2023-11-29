
#pragma once

#include <string.h>

#include "utils.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define IX(i,j) ((i)+(nx)*(j))
#define SWAP(x0,x) { double * tmp=x0; x0=x; x=tmp; }

typedef struct {
  uint nx, ny;
  double diff, visc, dt;
  uint nrelax;  // TODO : use a threshold instead if this is -1 ?
  double *rho, *rho_prev, *vx, *vx_prev, *vx_, *vy, *vy_prev, *vy_;
} Sim;


Sim* sim_new(uint nx, uint ny, double dt, double diff, double visc, uint nrelax);
void sim_free(Sim* sim);
void sim_print(const Sim* sim);
void dens_step ( uint nx, uint ny, uint nrelax, double * x, double * x0, double * u, double * v, double diff, double dt );
void vel_step ( uint nx, uint ny, uint nrelax, double * u, double * v, double * u0, double * v0, double visc, double dt );
void sim_step(Sim* sim);