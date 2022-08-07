#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "linalg.h"


#define QX(I, J) (qx[((J)-1)*(itot)+((I)-1)])

int fluid_compute_potential(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int isize = parallel_get_size(itot, mpisize, mpirank);
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double dx = param->dx;
  const double dy = param->dy;
  double *psi = fluid->psi;
  double *qx = NULL;
  double *qy = NULL;
  fftw_plan fwrd = NULL;
  fftw_plan bwrd = NULL;
  /* allocate matrix and create fftw plans */
  {
    /* ! allocate buffers (different alignments, datatypes) ! 19 ! */
    // x-aligned (decomposed in y)
    qx = common_calloc(itot*jsize, sizeof(double));
    // y-aligned (decomposed in x)
    qy = common_calloc(isize*jtot, sizeof(double));
    /* ! create fftw plans ! 14 ! */
    fftw_iodim dims[1], hdims[1];
    fftw_r2r_kind kind[1];
    dims[0].n  = itot;
    dims[0].is = 1;
    dims[0].os = 1;
    hdims[0].n  = jsize;
    hdims[0].is = itot;
    hdims[0].os = itot;
    // forward transform
    kind[0] = FFTW_REDFT10;
    fwrd = fftw_plan_guru_r2r(1, dims, 1, hdims, qx, qx, kind, FFTW_ESTIMATE);
    // backward transform
    kind[0] = FFTW_REDFT01;
    bwrd = fftw_plan_guru_r2r(1, dims, 1, hdims, qx, qx, kind, FFTW_ESTIMATE);
  }
  /* ! compute right-hand-side ! 17 ! */
  const double gamma = param->rkcoefs[rkstep].gamma;
  const double dt = param->dt;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      double ux_xm = UX(i  , j  );
      double ux_xp = UX(i+1, j  );
      double uy_ym = UY(i  , j  );
      double uy_yp = UY(i  , j+1);
      QX(i, j) =
        1./(gamma*dt)*(
         +(ux_xp-ux_xm)/dx
         +(uy_yp-uy_ym)/dy
        );
    }
  }
  /* ! project to wave space ! 1 ! */
  fftw_execute(fwrd);
  /* ! transpose x-aligned matrix to y-aligned matrix ! 1 ! */
  parallel_transpose(itot, jtot, sizeof(double), MPI_DOUBLE, qx, qy);
  /* solve linear systems */
  {
    double *tdm_l = common_calloc(jtot, sizeof(double));
    double *tdm_c = common_calloc(jtot, sizeof(double));
    double *tdm_u = common_calloc(jtot, sizeof(double));
    double *tdm_q = common_calloc(jtot, sizeof(double));
    for(int i=0; i<isize; i++){
      /* ! compute eigenvalue of this i position ! 4 ! */
      int ioffset = parallel_get_offset(itot, mpisize, mpirank);
      double eigenvalue = -4./pow(dx, 2.)*pow(
          sin(M_PI*(i+ioffset)/(2.*itot)),
          2.
      );
      /* ! initialise tri-diagonal matrix ! 5 ! */
      for(int j=0; j<jtot; j++){
        tdm_l[j] = 1./dy/dy;
        tdm_u[j] = 1./dy/dy;
        tdm_c[j] = -tdm_l[j]-tdm_u[j]+eigenvalue;
      }
      /* ! solve linear system ! 7 ! */
      for(int j=0; j<jtot; j++){
        tdm_q[j] = qy[i*jtot+j];
      }
      my_dgtsv_p(jtot, tdm_l, tdm_c, tdm_u, tdm_q);
      for(int j=0; j<jtot; j++){
        qy[i*jtot+j] = tdm_q[j];
      }
    }
    common_free(tdm_l);
    common_free(tdm_c);
    common_free(tdm_u);
    common_free(tdm_q);
  }
  /* ! transpose y-aligned matrix to x-aligned matrix ! 1 ! */
  parallel_transpose(jtot, itot, sizeof(double), MPI_DOUBLE, qy, qx);
  /* ! project to physical space ! 1 ! */
  fftw_execute(bwrd);
  /* ! normalise and store result ! 5 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      PSI(i, j) = QX(i, j)/(2.*itot);
    }
  }
  // NOTE: psi and p are assumed to have the same array shape
  fluid_update_boundaries_p(param, parallel, psi);
  /* deallocate memories */
  common_free(qx);
  common_free(qy);
  fftw_destroy_plan(fwrd);
  fftw_destroy_plan(bwrd);
  return 0;
}

#undef QX

