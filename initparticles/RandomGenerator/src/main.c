#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "fileio.h"


static double gen_random(const double min, const double max){
  return (max-min)*rand()/RAND_MAX+min;
}

static double fmin3(const double val0, const double val1, const double val2){
  double retval = val0;
  if(val1 < retval){
    retval = val1;
  }
  if(val2 < retval){
    retval = val2;
  }
  return retval;
}

static int check_stats(double *mindist, double *vfrac, const double lx, const double ly, const int n_particles, const double *rs, const double *xs, const double *ys){
  *mindist = DBL_MAX;
  for(int n0 = 0; n0 < n_particles; n0++){
    double r0 = rs[n0];
    double x0 = xs[n0];
    double y0 = ys[n0];
    for(int n1 = 0; n1 < n0; n1++){
      double r1 = rs[n1];
      double x1 = xs[n1];
      double y1 = ys[n1];
      double xvec = fabs(x0-x1);
      double yvec = fmin3(
          fabs(y0-y1-ly),
          fabs(y0-y1   ),
          fabs(y0-y1+ly)
      );
      double dist =
        sqrt(
            +pow(xvec, 2.)
            +pow(yvec, 2.)
        )-(r0+r1);
      if(dist < *mindist){
        *mindist = dist;
      }
    }
  }
  *vfrac = 0.;
  for(int n = 0; n < n_particles; n++){
    double r = rs[n];
    *vfrac += M_PI*pow(r, 2.);
  }
  *vfrac /= lx*ly;
  return 0;
}

static int output(const int n_particles, const double *dens, const double *rs, const double *xs, const double *ys, const double *uxs, const double *uys, const double *vzs){
  fileio_w_0d_serial("..", "n_particles",   NPYIO_INT,    sizeof(int),    &n_particles);
  fileio_w_1d_serial("..", "particle_dens", NPYIO_DOUBLE, sizeof(double), n_particles, dens);
  fileio_w_1d_serial("..", "particle_rs",   NPYIO_DOUBLE, sizeof(double), n_particles, rs  );
  fileio_w_1d_serial("..", "particle_xs",   NPYIO_DOUBLE, sizeof(double), n_particles, xs  );
  fileio_w_1d_serial("..", "particle_ys",   NPYIO_DOUBLE, sizeof(double), n_particles, ys  );
  fileio_w_1d_serial("..", "particle_uxs",  NPYIO_DOUBLE, sizeof(double), n_particles, uxs );
  fileio_w_1d_serial("..", "particle_uys",  NPYIO_DOUBLE, sizeof(double), n_particles, uys );
  fileio_w_1d_serial("..", "particle_vzs",  NPYIO_DOUBLE, sizeof(double), n_particles, vzs );
  return 0;
}

int main(void){
  srand(0);
  const double lx = 1.;
  const double ly = 1.;
  const int n_particles = 16;
  double *dens = common_calloc(n_particles, sizeof(double));
  double *rs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  // density
  for(int n = 0; n < n_particles; n++){
    dens[n] = 1.;
  }
  // radius
  for(int n = 0; n < n_particles; n++){
    rs[n] = gen_random(0.10, 0.10);
  }
  // position, no overlaps
  for(int n0 = 0; n0 < n_particles; n0++){
regen:
    {
      double r0 = rs[n0];
      double x0 = gen_random(rs[n0], lx-rs[n0]);
      double y0 = gen_random(0.,     ly       );
      for(int n1 = 0; n1 < n0; n1++){
        double r1 = rs[n1];
        double x1 = xs[n1];
        double y1 = ys[n1];
        double xvec = fabs(x0-x1);
        double yvec = fmin3(
            fabs(y0-y1-ly),
            fabs(y0-y1   ),
            fabs(y0-y1+ly)
        );
        double dist =
          sqrt(
              +pow(xvec, 2.)
              +pow(yvec, 2.)
          )-(r0+r1);
        if(dist < 0.){
          goto regen;
        }
      }
      xs[n0] = x0;
      ys[n0] = y0;
      printf("particle %*d @ % .3f % .3f\n", 5, n0, x0, y0);
    }
  }
  {
    double mindist, vfrac;
    check_stats(&mindist, &vfrac, lx, ly, n_particles, rs, xs, ys);
    printf("minimum distance (should be positive): % .1e\n", mindist);
    printf("volume fraction                      : % .1e\n", vfrac);
  }
  // velocities, still
  memset(uxs, 0, sizeof(double)*n_particles);
  memset(uys, 0, sizeof(double)*n_particles);
  memset(vzs, 0, sizeof(double)*n_particles);
  // output
  output(n_particles, dens, rs, xs, ys, uxs, uys, vzs);
  // clean-up buffers
  common_free(dens);
  common_free(rs);
  common_free(xs);
  common_free(ys);
  common_free(uxs);
  common_free(uys);
  common_free(vzs);
  return 0;
}

