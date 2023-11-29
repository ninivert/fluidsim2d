#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "mpi.h"

#include "solver.h"
#include "spinny.h"
#include "utils.h"
#include "main_utils.h"

enum Role { ROLE_DENS, ROLE_VEL };

int main(int argc, char** argv) {
  // init mpi
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
  if (provided < MPI_THREAD_SINGLE) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // MPI_Init(NULL, NULL);
  int size, rank, rank_of_other;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // TODO: assert that MPI has been called with np=2

  // assign roles
  enum Role role = rank == 0 ? ROLE_DENS : ROLE_VEL;
  rank_of_other = rank == 0 ? 1 : 0;
  printf("[mpi %d/%d] hello from role = %d\n", rank, size, role);

  // openmp init
#ifdef _OPENMP
  // uint max_threads = omp_get_max_threads();
  // double alpha0 = 0.25;
  // if (argc == 4) { alpha0 = atof(argv[3]); }
  // alpha0 = clamp(alpha0, 0.0, 1.0);
  // const double alpha = role == ROLE_DENS ? alpha0 : 1-alpha0;
  // uint nthreads = clampi((int) (alpha*(double) max_threads), 1, max_threads-1);
  // printf("[mpi %d/%d] omp max threads=%d, running with %d (alpha=%f)\n", rank, size, max_threads, nthreads, alpha);
  uint nthreads_dens = 2, nthreads_vel = 6;
  if (argc == 5) { nthreads_dens = atoi(argv[3]); nthreads_vel = atoi(argv[4]); }
  uint nthreads = role == ROLE_DENS ? nthreads_dens : nthreads_vel;
  omp_set_num_threads(nthreads);
  printf("[mpi %d/%d] omp running with %d threads\n", rank, size, nthreads);
#else
  printf("[no omp]\n");
#endif

  // init simulation
  const uint nsteps = 1000;
  uint nx = 100, ny = 100;
  if (argc >= 3) { nx = atoi(argv[1]); ny = atoi(argv[2]); }

  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 100);
  if (rank == 0) sim_print(sim);

  MPI_Request request_v[2];

  long long t0 = current_timestamp();
  long long totaltime = 0;

  for (uint istep = 0; istep < nsteps; ++istep) {
    if (istep % 100 == 0) printf("[sim %d/%d] iter %u/%u\n", rank, size, istep, nsteps);

    long long tstep0 = current_timestamp();
    // apply sources and step
    if (role == ROLE_DENS) {
      spinny_dens(sim);
      dens_step(sim->nx, sim->ny, sim->nrelax, sim->rho, sim->rho_prev, sim->vx, sim->vy, sim->diff, sim->dt);
    } else {
      spinny_vel(sim);
      vel_step(sim->nx, sim->ny, sim->nrelax, sim->vx, sim->vy, sim->vx_prev, sim->vy_prev, sim->visc, sim->dt);
    }
    totaltime += current_timestamp() - tstep0;

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

  printf("[sim %d/%d] finished in %llu milliseconds, avg local compute time per step = %f ms\n", rank, size, current_timestamp() - t0, (double) totaltime / (double) nsteps);
  // only the density thread has the full information
  if (role == ROLE_DENS) write_results(sim);
  sim_free(sim);

  MPI_Finalize();

  return 0;
}