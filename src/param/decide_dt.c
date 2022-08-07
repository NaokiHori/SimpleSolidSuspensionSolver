#include <math.h>
#include <float.h>
#include <mpi.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"


int param_decide_dt(param_t *param, const parallel_t *parallel, const fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double Re = param->Re;
  const double lx = param->lx;
  const double *dxc = param->dxc;
  const double dy = param->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  /* ! advective constraint, which will be used as dt ! 25 ! */
  double dt_adv;
  {
    // sufficiently small number
    const double small = 1.e-8;
    // grid-size/velocity
    dt_adv = 5.e-2; // max dt_adv
    // ux
    for(int j=1; j<=jsize; j++){
      for(int i=2; i<=itot; i++){
        // to avoid zero-division
        double vel = fmax(fabs(UX(i, j)), small);
        dt_adv = fmin(dt_adv, DXC(i)/vel);
      }
    }
    // uy
    for(int j=1; j<=jsize; j++){
      for(int i=1; i<=itot; i++){
        // to avoid zero-division
        double vel = fmax(fabs(UY(i, j)), small);
        dt_adv = fmin(dt_adv, dy/vel);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &dt_adv, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_adv *= param->safefactors[0];
  }
  /* ! diffusive constraints in each direction ! 13 ! */
  double dt_dif_x, dt_dif_y;
  {
    // find minimum grid size in x direction
    double dx;
    dx = lx;
    for(int i=2; i<=itot; i++){
      dx = fmin(dx, DXC(i));
    }
    dt_dif_x = 0.25*Re*pow(dx, 2.);
    dt_dif_y = 0.25*Re*pow(dy, 2.);
    dt_dif_x *= param->safefactors[1];
    dt_dif_y *= param->safefactors[1];
  }
  /* ! time step size is assigned ! 1 ! */
  param->dt = fmin(dt_adv, fmin(dt_dif_x, dt_dif_y));
  return 0;
}

