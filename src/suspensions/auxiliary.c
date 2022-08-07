#include <math.h>
#include "common.h"
#include "suspensions.h"


// number of grid points outside circles
#define NEXTRA 3
#define POW2(x) ((x)*(x))

int suspensions_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max){
  *min = (int)(floor((grav_center-radius)/grid_size)-NEXTRA);
  *max = (int)(ceil ((grav_center+radius)/grid_size)+NEXTRA);
  *min = *min < lbound ? lbound : *min;
  *max = *max > ubound ? ubound : *max;
  return 0;
}

#define BETA 2.

static inline double compute_dist(const double grid_size, const double radius, const double px, const double py, const double x, const double y){
  double deltax = x-px;
  double deltay = y-py;
  double dist = radius-sqrt(POW2(deltax)+POW2(deltay));
  return dist/grid_size;
}

double suspensions_s_weight(const double grid_size, const double radius, const double px, const double py, const double x, const double y){
  double dist = compute_dist(grid_size, radius, px, py, x, y);
  return 0.5*BETA*(1.-POW2(tanh(BETA*dist)));
}

double suspensions_v_weight(const double grid_size, const double radius, const double px, const double py, const double x, const double y){
  double dist = compute_dist(grid_size, radius, px, py, x, y);
  return 0.5*(1.+tanh(BETA*dist));
}

#undef BETA

double suspensions_compute_volume(const double r){
  return M_PI*POW2(r);
}

double suspensions_compute_mass(const double den, const double r){
  double vol = suspensions_compute_volume(r);
  return den*vol;
}

double suspensions_compute_moment_of_inertia(const double den, const double r){
  double mass = suspensions_compute_mass(den, r);
  return 0.5*mass*POW2(r);
}

#undef NEXTRA
#undef POW2

