#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"


static int increment_particle_velocities(const param_t *param, const int rkstep, suspensions_t *suspensions){
  const double Fr = param->Fr;
  const double dt = param->dt;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    const double den = p->den;
    p->dux =
      +0.5*gamma*dt*(p->cfx[0]+p->cfx[1]) // collision
      +          dt* p->fux               // boundary force
      +(p->iux[1]-p->iux[0])              // internal inertia
      +1.0*gamma*dt*(1.-den)/pow(Fr, 2.)  // buoyancy
      ;
    p->duy =
      +0.5*gamma*dt*(p->cfy[0]+p->cfy[1]) // collision
      +          dt* p->fuy               // boundary force
      +(p->iuy[1]-p->iuy[0])              // internal inertia
      ;
    p->dvz =
      +          dt* p->tvz               // boundary torque
      +(p->ivz[1]-p->ivz[0])
      ;
  }
  return 0;
}

static int increment_particle_positions(const param_t *param, const int rkstep, suspensions_t *suspensions){
  const double dt = param->dt;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    double ux0 = p->ux;
    double uy0 = p->uy;
    double ux1 = p->ux+p->dux;
    double uy1 = p->uy+p->duy;
    p->dx = 0.5*gamma*dt*(ux0+ux1);
    p->dy = 0.5*gamma*dt*(uy0+uy1);
  }
  return 0;
}

int suspensions_increment_particles(const param_t *param, const int rkstep, suspensions_t *suspensions){
  increment_particle_velocities(param, rkstep, suspensions);
  increment_particle_positions(param, rkstep, suspensions);
  return 0;
}

