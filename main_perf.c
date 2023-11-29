#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "mpi.h"

#include "solver.c"
#include "utils.h"

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
  return 0;
}

enum Role { ROLE_DENS, ROLE_VEL };

int main(int argc, char** argv) {
  // init mpi
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int size, rank, rank_of_other;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  // openmp greeting
  #ifdef _OPENMP
    #pragma omp single
    printf("[omp] max threads:%d\n", omp_get_max_threads());
  #endif

  enum Role role = rank == 0 ? ROLE_DENS : ROLE_VEL;
  rank_of_other = rank == 0 ? 1 : 0;
  printf("[mpi %d/%d] role = %d\n", rank, size, role);

  // init simulation
  const uint nsteps = 1000;
  uint nx = 100, ny = 100;

  if (argc == 3) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  } else {
    if (rank == 0) {
      printf("usage: ./main_perf nx ny\n");
      printf("using default values.\n");
    }
  }

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 20);
  if (rank == 0) {
    sim_print(sim);
  }

  long long t0 = current_timestamp();

  MPI_Request request_v[2];

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("[mpi %d/%d] %u / %u\n", rank, size, istep, nsteps);

    // do a step
    if (role == ROLE_DENS) {
      for (size_t i = 0; i < nx*ny; ++i) sim->rho_prev[i] = 0.0;
      sim->rho_prev[IX(nx/2+10,ny/2+10)] = 10.0;
      dens_step(sim->nx, sim->ny, sim->nrelax, sim->rho, sim->rho_prev, sim->vx, sim->vy, sim->diff, sim->dt);
    } else {
      for (size_t y = 0; y < sim->ny; ++y) { for (size_t x = 0; x < sim->nx; ++x) {
        double xf = linterp(x, 0, nx-1, -1, 1), yf = linterp(y, 0, ny-1, -1, 1);
        double r = sqrt(xf*xf + yf*yf), theta = atan2(yf, xf);
        sim->vx_prev[IX(x, y)] = 0.01*r*sin(theta);
        sim->vy_prev[IX(x, y)] = -0.01*r*cos(theta);
      }}
      vel_step(sim->nx, sim->ny, sim->nrelax, sim->vx, sim->vy, sim->vx_prev, sim->vy_prev, sim->visc, sim->dt);
    }

    // vel process needs to update the vx,vy buffers of the density process
    // but vel is indep of density
    if (role == ROLE_DENS) {
      MPI_Irecv(sim->vx, sim->nx*sim->ny, MPI_DOUBLE, rank_of_other, 0, MPI_COMM_WORLD, &request_v[0]);
      MPI_Irecv(sim->vy, sim->nx*sim->ny, MPI_DOUBLE, rank_of_other, 0, MPI_COMM_WORLD, &request_v[1]);
    } else {
      MPI_Isend(sim->vx, sim->nx*sim->ny, MPI_DOUBLE, rank_of_other, 0, MPI_COMM_WORLD, &request_v[0]);
      MPI_Isend(sim->vy, sim->nx*sim->ny, MPI_DOUBLE, rank_of_other, 0, MPI_COMM_WORLD, &request_v[1]);
    }
    MPI_Waitall(2, request_v, MPI_STATUS_IGNORE);
  }

  printf("finished in %llu milliseconds\n", current_timestamp() - t0);

  // only the density thread has the full information
  if (role == ROLE_DENS) {
    write_results(sim);
  }

  sim_free(sim);

  MPI_Finalize();

  return 0;
}