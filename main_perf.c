#include <stdio.h>
#include <stdlib.h>

#include "solver.c"
#include "utils.h"

int main(int argc, char** argv) {
  const nsteps = 1000;
  const uint nx = 200, ny = 100;
  Sim* sim = sim_new(nx, ny, 0.1, 0.00001, 0.0, 20);
  sim_print(sim);



  sim_free(sim);

  return 0;
}