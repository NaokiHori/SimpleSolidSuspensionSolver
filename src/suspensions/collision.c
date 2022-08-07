#include <stdio.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "suspensions.h"


static int get_l(const int n_particles, const int n0){
  return (2*n_particles-n0-1)*(n0  )/2  ;
}

static int get_r(const int n_particles, const int n0){
  return (2*n_particles-n0-2)*(n0+1)/2-1;
}

static int get_particle_indices(const int n_particles, const int n, int *n0, int *n1){
  for(*n0 = 0; ; (*n0)++){
    if(get_l(n_particles, *n0) <= n && n <= get_r(n_particles, *n0)){
      break;
    }
  }
  *n1 = n-get_l(n_particles, *n0)+(*n0)+1;
  return 0;
}

static int get_my_range(const parallel_t *parallel, const int n_total, int *n_min, int *n_max){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  *n_min = parallel_get_offset(n_total, mpisize, mpirank);
  *n_max = parallel_get_size(n_total, mpisize, mpirank) + *n_min;
  return 0;
}

static double harmonic_average(const double v0, const double v1){
  double retval = 0.5*(1./v0+1./v1);
  return 1./retval;
}

static int compute_collision_force_p_p(const param_t *param, const int cnstep, particle_t *p0, particle_t *p1){
  const double p0den  = p0->den;
  const double p0r    = p0->r;
  const double p0mass = suspensions_compute_mass(p0den, p0r);
  const double p0x    = p0->x+p0->dx;
  const double p0y    = p0->y+p0->dy;
  const double p1den  = p1->den;
  const double p1r    = p1->r;
  const double p1x    = p1->x+p1->dx;
  const double p1y    = p1->y+p1->dy;
  const double p1mass = suspensions_compute_mass(p1den, p1r);
  // k: pre-factor (spring stiffness)
  const double mass = harmonic_average(p0mass, p1mass);
  const double reft = pow(harmonic_average(p0r, p1r), 1.5);
  const double k = mass*pow(M_PI, 2.)/pow(reft, 2.);
  // compute normal vector from particle 0 to 1 (N.B. periodicity in y)
  const double ly = param->ly;
  double nx, ny, norm;
  nx = p1x-p0x;
  norm = DBL_MAX;
  for(int periodic = -1; periodic <= 1; periodic++){
    double ny_ = p1y-p0y+ly*periodic;
    double norm_ = pow(nx, 2.)+pow(ny_, 2.);
    if(norm_ < norm){
      norm = norm_;
      ny = ny_;
    }
  }
  norm = sqrt(norm);
  nx /= norm;
  ny /= norm;
  double overlap_dist = p0r+p1r-norm;
  // impose force when overlapped
  if(overlap_dist > 0.){
    double cfx = k*overlap_dist*nx;
    double cfy = k*overlap_dist*ny;
    p0->cfx[cnstep] -= 1./p0mass*cfx;
    p0->cfy[cnstep] -= 1./p0mass*cfy;
    p1->cfx[cnstep] += 1./p1mass*cfx;
    p1->cfy[cnstep] += 1./p1mass*cfy;
  }
  return 0;
}

static int compute_collision_force_p_w(const double wallx, const int cnstep, particle_t *p){
  const double pden  = p->den;
  const double pr    = p->r;
  const double pmass = suspensions_compute_mass(pden, pr);
  const double px    = p->x+p->dx;
  // k: pre-factor (spring stiffness)
  const double mass = pmass;
  const double reft = pow(pr, 1.5);
  const double k = mass*pow(M_PI, 2.)/pow(reft, 2.);
  // compute normal vector from particle to the wall
  double nx = wallx-px;
  double norm = fabs(nx);
  nx /= norm;
  double overlap_dist = pr-norm;
  // impose force when overlapped
  if(overlap_dist > 0.){
    double cfx = k*overlap_dist*nx;
    p->cfx[cnstep] -= 1./pmass*cfx;
  }
  return 0;
}

int suspensions_compute_collision_force(const param_t *param, const parallel_t *parallel, const int cnstep, suspensions_t *suspensions){
  /*
   * NOTE: only the spring in the normal direction is considered for simplicity
   * Although this is sufficient to avoid over-penetrations between particles,
   *   obviously not collect from a physical perspective
   */
  const double lx = param->lx;
  const int n_particles = suspensions->n_particles;
  particle_t **particles = suspensions->particles;
  // reset collision force at this CN step
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->cfx[cnstep] = 0.;
    p->cfy[cnstep] = 0.;
  }
  // particle-particle collisions
  {
    const int n_total = n_particles*(n_particles-1)/2;
    int n_min, n_max;
    get_my_range(parallel, n_total, &n_min, &n_max);
    for(int n = n_min; n < n_max; n++){
      int n0, n1;
      get_particle_indices(n_particles, n, &n0, &n1);
      compute_collision_force_p_p(param, cnstep, particles[n0], particles[n1]);
    }
  }
  // particle-wall collisions
  {
    double wall_locations[2] = {0., lx};
    for(int wall_index = 0; wall_index < 2; wall_index++){
      double wall_location = wall_locations[wall_index];
      int n_min, n_max;
      get_my_range(parallel, n_particles, &n_min, &n_max);
      for(int n = n_min; n < n_max; n++){
        compute_collision_force_p_w(wall_location, cnstep, particles[n]);
      }
    }
  }
  // synchronise computed forcings
  {
    // prepare message buffer
    double *buffer = common_calloc(2*n_particles, sizeof(double));
    // pack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      buffer[2*n+0] = p->cfx[cnstep];
      buffer[2*n+1] = p->cfy[cnstep];
    }
    // sum up all
    MPI_Allreduce(MPI_IN_PLACE, buffer, 2*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // unpack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      p->cfx[cnstep] = buffer[2*n+0];
      p->cfy[cnstep] = buffer[2*n+1];
    }
    // clean-up buffer
    common_free(buffer);
  }
  return 0;
}

