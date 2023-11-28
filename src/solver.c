#pragma once

#include "utils.h"

#define IX(i,j) ((i)+(N)*(j))
#define SWAP(x0,x) { double * tmp=x0; x0=x; x=tmp; }
#define FOR_EACH_CELL for ( i=1 ; i<N-1 ; i++ ) { for ( j=1 ; j<N-1 ; j++ ) {
#define END_FOR }}

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
  sim->rho = malloc(nx*ny*sizeof(*sim->rho));
  sim->rho_prev = malloc(nx*ny*sizeof(*sim->rho_prev));
  sim->vx = malloc(nx*ny*sizeof(*sim->vx));
  sim->vx_prev = malloc(nx*ny*sizeof(*sim->vx_prev));
  sim->vy = malloc(nx*ny*sizeof(*sim->vy));
  sim->vy_prev = malloc(nx*ny*sizeof(*sim->vy_prev));
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

// uint unravel(uint x, uint y, Arr2d* arr) {
//   return y*arr->nx + x;
// }
// #define A(arr, x, y) arr->data[unravel(x, y, arr)]

void add_source (uint N, double *x, const double *s, double dt) {
  int i, size=N*N; // todo nx*ny
  for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

void set_bnd ( int N, int b, double * x )
{
  int i;

  for ( i=1 ; i<N-1 ; i++ ) {
    x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
    x[IX(N-1,i)] = b==1 ? -x[IX(N-2,i)] : x[IX(N-2,i)];
    x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
    x[IX(i,N-1)] = b==2 ? -x[IX(i,N-2)] : x[IX(i,N-2)];
  }
  x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
  x[IX(0  ,N-1)] = 0.5f*(x[IX(1,N-1)]+x[IX(0  ,N-2)]);
  x[IX(N-1,0  )] = 0.5f*(x[IX(N-2,0  )]+x[IX(N-1,1)]);
  x[IX(N-1,N-1)] = 0.5f*(x[IX(N-2,N-1)]+x[IX(N-1,N-2)]);
}

void lin_solve ( int N, int b, double * x, double * x0, double a, double c )
{
  int i, j, k;

  for ( k=0 ; k<20 ; k++ ) {
    FOR_EACH_CELL
      x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
    END_FOR
    set_bnd ( N, b, x );
  }
}

void diffuse ( int N, int b, double * x, double * x0, double diff, double dt )
{
  double a=dt*diff*N*N;
  lin_solve ( N, b, x, x0, a, 1+4*a );
}

void advect ( int N, int b, double * d, double * d0, double * u, double * v, double dt )
{
  int i, j, i0, j0, i1, j1;
  double x, y, s0, t0, s1, t1, dt0;

  dt0 = dt*N;
  FOR_EACH_CELL
    x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
    if (x<0.5f) x=0.5f; if (x>N-1.5f) x=N-1.5f; i0=(int)x; i1=i0+1;
    if (y<0.5f) y=0.5f; if (y>N-1.5f) y=N-1.5f; j0=(int)y; j1=j0+1;
    s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
    d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
           s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
  END_FOR
  set_bnd ( N, b, d );
}

void project ( int N, double * u, double * v, double * p, double * div )
{
  int i, j;

  FOR_EACH_CELL
    div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
    p[IX(i,j)] = 0;
  END_FOR	
  set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

  lin_solve ( N, 0, p, div, 1, 4 );

  FOR_EACH_CELL
    u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
    v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
  END_FOR
  set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

void dens_step ( int N, double * x, double * x0, double * u, double * v, double diff, double dt )
{
  add_source ( N, x, x0, dt );
  SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt );
  SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt );
}

void vel_step ( int N, double * u, double * v, double * u0, double * v0, double visc, double dt )
{
  add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
  SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
  SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
  project ( N, u, v, u0, v0 );
  SWAP ( u0, u ); SWAP ( v0, v );
  advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
  project ( N, u, v, u0, v0 );
}


void sim_step(Sim* sim) {
  vel_step(sim->nx, sim->vx, sim->vy, sim->vx_prev, sim->vy_prev, sim->visc, sim->dt);
  dens_step(sim->nx, sim->rho, sim->rho_prev, sim->vx, sim->vy, sim->diff, sim->dt);
}