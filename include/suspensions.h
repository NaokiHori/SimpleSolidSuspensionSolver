#if !defined(SUSPENSIONS_H)
#define SUSPENSIONS_H

#include "structure.h"
#include "arrays/suspensions.h"

#define CNSTEPMAX 2

typedef struct particle_t_ {
  // fixed parameters
  // density (ratio), radius
  double den, r;
  // gravity center
  double x, y;
  double dx, dy;
  // translational and rotational velocities
  double ux, uy, vz;
  double dux, duy, dvz;
  // surface forces and torque at k+1/2 step
  // 1/m * \int_{S} a_i dS or 1/I * \int_{S} \epsilon_{ijk} \omega_j a_k dS
  // a_i^k = \alpha^k ( U_i^k - u_i^* ) / ( \gamma \Delta t )
  double fux, fuy, tvz;
  // internal inertia, e.g.,
  //   0: 1/C * \int_{V^{k  }} u_i^{k } dV^{k  }
  //   1: 1/C * \int_{V^{k+1}} u_i^{**} dV^{k+1}
  // similar to the rotational component
  // C is a pre-factor, mass or momentum of inertia
  double iux[CNSTEPMAX], iuy[CNSTEPMAX], ivz[CNSTEPMAX];
  // collision force based on k (0) and k+1 (1) step particle positions and velocities
  double cfx[CNSTEPMAX], cfy[CNSTEPMAX];
} particle_t;

struct suspensions_t_ {
  int n_particles;
  particle_t **particles;
  // responses of surface forces and torque on the momentum fields
  double *dux, *duy;
};

/* constructor and destructor */
extern suspensions_t *suspensions_init(const param_t *param, const parallel_t *parallel);
extern int suspensions_finalise(suspensions_t *suspensions);

/* called by main update routine */
extern int suspensions_reset_particle_increments(suspensions_t *suspensions);
extern int suspensions_compute_inertia(const param_t *param, const parallel_t *parallel, const int cnstep, const fluid_t *fluid, suspensions_t *suspensions);
extern int suspensions_compute_collision_force(const param_t *param, const parallel_t *parallel, const int cnstep, suspensions_t *suspensions);
extern int suspensions_exchange_momentum(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, suspensions_t *suspensions);
extern int suspensions_increment_particles(const param_t *param, const int rkstep, suspensions_t *suspensions);
extern int suspensions_update_momentum_fleid(const param_t *param, const parallel_t *parallel, fluid_t *fluid, const suspensions_t *suspensions);
extern int suspensions_update_particles(const param_t *param, suspensions_t *suspensions);

/* other supportive functions */
extern int suspensions_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max);
extern double suspensions_compute_volume(const double r);
extern double suspensions_compute_mass(const double den, const double r);
extern double suspensions_compute_momentum_of_inertia(const double den, const double r);

extern double suspensions_s_weight(const double grid_size, const double radius, const double px, const double py, const double x, const double y);
extern double suspensions_v_weight(const double grid_size, const double radius, const double px, const double py, const double x, const double y);

#endif // SUSPENSIONS_H
