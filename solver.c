// 2D incompressible fluid solver
// based on the work of : Jos Stam (jstam@aw.sgi.com), Creation Date : Jan 9 2003

#pragma once

#include "utils.h"
// #ifdef _OPENMP
// #include "omp.h"
// #endif

#define IX(i,j) ((i)+(nx)*(j))
#define SWAP(x0,x) { double * tmp=x0; x0=x; x=tmp; }

typedef struct {
  uint nx, ny;
  double diff, visc, dt;
  uint nrelax;  // TODO : use a threshold instead if this is -1 ?
  double *rho, *rho_prev, *vx, *vx_prev, *vy, *vy_prev;
} Sim;

Sim* sim_new(uint nx, uint ny, double dt, double diff, double visc, uint nrelax) {
  Sim* sim = malloc(sizeof(*sim));
  sim->nx = nx;
  sim->ny = ny;
  sim->dt = dt;
  sim->diff = diff;
  sim->visc = visc;
  sim->nrelax = nrelax;
  sim->vx = calloc(nx*ny, sizeof(*sim->vx));
  sim->vx_prev = calloc(nx*ny, sizeof(*sim->vx_prev));
  sim->vy = calloc(nx*ny, sizeof(*sim->vy));
  sim->vy_prev = calloc(nx*ny, sizeof(*sim->vy_prev));
  sim->rho = calloc(nx*ny, sizeof(*sim->rho));
  sim->rho_prev = calloc(nx*ny, sizeof(*sim->rho_prev));
  return sim;
}

void sim_free(Sim* sim) {
  free(sim->rho);
  free(sim->rho_prev);
  free(sim->vx);
  free(sim->vx_prev);
  free(sim->vy);
  free(sim->vy_prev);
  free(sim);
}

void sim_print(const Sim* sim) {
  printf("simulation with nx=%ui, ny=%ui, dt=%f, diff=%f, visc=%f\n", sim->nx, sim->ny, sim->dt, sim->diff, sim->visc);
}

void add_source (uint nx, uint ny, double *x, const double *s, double dt) {
  for (uint i=0 ; i < nx*ny ; ++i) x[i] += dt*s[i];
}

// TODO : use function pointers with callbacks instead of `b` to impose bcs
void set_bnd (uint nx, uint ny, int b, double *x) {
  double mult_lr = b == 1 ? -1.0 : 1.0;
  double mult_tb = b == 2 ? -1.0 : 1.0;

  // top and bottom
  for (uint i = 1; i < nx-1; ++i) {
    x[IX(i,0   )] = mult_tb * x[IX(i,   1)];
    x[IX(i,ny-1)] = mult_tb * x[IX(i,ny-2)];
  }

  // left and right
  for (uint j = 1; j < ny-1; ++j) {
    x[IX(0   ,j)] = mult_lr * x[IX(1,   j)];
    x[IX(nx-1,j)] = mult_lr * x[IX(nx-2,j)];
  }

  // set corners
  x[IX(0  , 0   )] = 0.5f*(x[IX(1,0      )]+x[IX(0   ,1   )]);
  x[IX(0  , ny-1)] = 0.5f*(x[IX(1,ny-1   )]+x[IX(0   ,ny-2)]);
  x[IX(nx-1,0   )] = 0.5f*(x[IX(nx-2,0   )]+x[IX(nx-1,1   )]);
  x[IX(nx-1,nx-1)] = 0.5f*(x[IX(nx-2,ny-1)]+x[IX(nx-1,ny-2)]);
}

void lin_solve (uint nx, uint ny, uint nrelax, int b, double * x, const double * x0, double a, double c ) {
  // int num_threads = omp_get_max_threads();
  for (uint k=0 ; k<nrelax ; ++k) {
    // #pragma omp parallel for schedule(static, (nx-2)/num_threads)
    for (uint i=1 ; i<nx-1 ; ++i) { for (uint j=1 ; j<ny-1 ; ++j) {
        x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
    }}
    set_bnd (nx, ny, b, x);
  }
}

void diffuse (uint nx, uint ny, uint nrelax, int b, double * x, const double * x0, double diff, double dt) {
  double a=dt*diff*nx*ny;
  lin_solve ( nx, ny, nrelax, b, x, x0, a, 1+4*a );
}

void advect (uint nx, uint ny, int b, double * d, const double * d0, const double * u, const double * v, double dt) {
  int i0, j0, i1, j1;
  double x, y, s0, t0, s1, t1;
  double dtx0 = dt*ny, dty0 = dt*ny;

  for (uint i=1 ; i<nx-1 ; ++i) { for (uint j=1 ; j<ny-1 ; ++j) {
    x = i-dtx0*u[IX(i,j)]; y = j-dty0*v[IX(i,j)];
    x = clamp(x, 0.5, nx-1.5); y = clamp(y, 0.5, ny-1.5);
    i0=(int)x; i1=i0+1;
    j0=(int)y; j1=j0+1;
    s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
    d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
  }}
  set_bnd (nx, ny, b, d);
}

void project (uint nx, uint ny, uint nrelax, double *u, double *v, double *p, double *div) {
  for (uint i=1 ; i<nx-1 ; ++i) { for (uint j=1 ; j<ny-1 ; ++j) {
    div[IX(i,j)] = -0.5*((u[IX(i+1,j)]-u[IX(i-1,j)])/nx+(v[IX(i,j+1)]-v[IX(i,j-1)])/ny);
    p[IX(i,j)] = 0.0;
  }}	
  set_bnd (nx, ny, 0, div); set_bnd (nx, ny, 0, p);
  lin_solve (nx, ny, nrelax, 0, p, div, 1.0, 4.0);

  for (uint i=1 ; i<nx-1 ; ++i) { for (uint j=1 ; j<ny-1 ; ++j) {
    u[IX(i,j)] -= 0.5*nx*(p[IX(i+1,j)]-p[IX(i-1,j)]);
    v[IX(i,j)] -= 0.5*ny*(p[IX(i,j+1)]-p[IX(i,j-1)]);
  }}
  set_bnd ( nx, ny, 1, u ); set_bnd ( nx, ny, 2, v );
}

void dens_step ( uint nx, uint ny, uint nrelax, double * x, double * x0, double * u, double * v, double diff, double dt ) {
  add_source ( nx, ny, x, x0, dt );
  SWAP ( x0, x ); diffuse ( nx, ny, nrelax, 0, x, x0, diff, dt );
  SWAP ( x0, x ); advect ( nx, ny, 0, x, x0, u, v, dt );
}

void vel_step ( uint nx, uint ny, uint nrelax, double * u, double * v, double * u0, double * v0, double visc, double dt ) {
  add_source ( nx, ny, u, u0, dt ); add_source ( nx, ny, v, v0, dt );
  SWAP ( u0, u ); diffuse ( nx, ny, nrelax, 1, u, u0, visc, dt );
  SWAP ( v0, v ); diffuse ( nx, ny, nrelax, 2, v, v0, visc, dt );
  project ( nx, ny, nrelax, u, v, u0, v0 );
  SWAP ( u0, u ); SWAP ( v0, v );
  advect ( nx, ny, 1, u, u0, u0, v0, dt ); advect ( nx, ny, 2, v, v0, u0, v0, dt );
  project ( nx, ny, nrelax, u, v, u0, v0 );
}


void sim_step(Sim* sim) {
  vel_step(sim->nx, sim->ny, sim->nrelax, sim->vx, sim->vy, sim->vx_prev, sim->vy_prev, sim->visc, sim->dt);
  dens_step(sim->nx, sim->ny, sim->nrelax, sim->rho, sim->rho_prev, sim->vx, sim->vy, sim->diff, sim->dt);
}