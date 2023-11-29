#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "omp.h"
#include "mpi.h"

#include "solver.h"
#include "spinny.h"
#include "utils.h"
#include "main_utils.h"

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

  // assign roles
  enum Role role = rank == 0 ? ROLE_DENS : ROLE_VEL;
  rank_of_other = rank == 0 ? 1 : 0;
  printf("[mpi %d/%d] role = %d\n", rank, size, role);

  // openmp greeting
  #pragma omp single
  printf("[omp] num threads:%d\n", omp_get_num_threads());

  // init simulation
  const uint nsteps = 1000;
  uint nx = 100, ny = 100;

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 20);
  if (rank == 0) sim_print(sim);

  MPI_Request request_v[2];

  long long t0 = current_timestamp();

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("[sim %d/%d] iter %u/%u\n", rank, size, istep, nsteps);

    // apply sources and step
    if (role == ROLE_DENS) {
      spinny_dens(sim);
      dens_step(sim->nx, sim->ny, sim->nrelax, sim->rho, sim->rho_prev, sim->vx, sim->vy, sim->diff, sim->dt);
    } else {
      spinny_vel(sim);
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

  printf("[sim] finished in %llu milliseconds\n", current_timestamp() - t0);
  // only the density thread has the full information
  if (role == ROLE_DENS) write_results(sim);
  sim_free(sim);

  MPI_Finalize();

  return 0;
}