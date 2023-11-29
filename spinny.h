#pragma once

#include <math.h>

#include "solver.h"

static void spinny_dens(Sim* sim) {
  uint nx = sim->nx, ny = sim->ny;
  for (uint i = 0; i < nx*ny; ++i) sim->rho_prev[i] = 0.0;
  sim->rho_prev[IX(nx/2+10,ny/2+10)] = 10.0;
}

static void spinny_vel(Sim* sim) {
  uint nx = sim->nx, ny = sim->ny;
  for (size_t y = 0; y < ny; ++y) {
    for (size_t x = 0; x < nx; ++x) {
      double xf = linterp(x, 0, nx-1, -1, 1), yf = linterp(y, 0, ny-1, -1, 1);
      double r = sqrt(xf*xf + yf*yf), theta = atan2(yf, xf);
      sim->vx_prev[IX(x, y)] = 0.01*r*sin(theta);
      sim->vy_prev[IX(x, y)] = -0.01*r*cos(theta);
    }
  }
}

static void spinny(Sim* sim) {
  spinny_vel(sim);
  spinny_dens(sim);
}