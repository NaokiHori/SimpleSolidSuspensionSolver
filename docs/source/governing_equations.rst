
.. include:: /references.txt

###################
Governing equations
###################

In addition to the Navier-Stokes equations

.. math::
   \der{u_i}{x_i} &= 0, \\
   \der{u_i}{t}
   +
   u_j \der{u_i}{x_j}
   &=
   -
   \der{p}{x_i}
   +
   \frac{1}{Re} \der{}{x_j} \der{u_i}{x_j}
   +
   f_i,

which govern the behaviour of the liquid, we consider the Newton-Euler equations to describe the motions of rigid particles.
In particular, the translational motion (velocity of the gravity center) is given by

.. math::

   \rho_p V_p \der{U_i}{t}
   =
   -
   \int_{\partial V_p} a_i d S_p^k
   +
   \der{}{t} \int_{V_p} u_i d V_p
   +
   F_i.

Taking the outer product with :math:`r_i` yields the equation for the rotational motion.

:math:`a_i` comes from the fluid-particle interactions (derivation can be found in e.g., |BREUGEM2012|).
There are several ways to formulate :math:`a_i`, and one of the most successful ways is the direct forcing on Lagrange points defined on each particle surface (|UHLMANN2005|).
Here, on the other hand, we exchange the momentum between particles and fluid directly on the Eulerian field (e.g., |KAJISHIMA2002|):

.. math::

   a_i
   \equiv
   w
   \frac{
      U_i + \epsilon_{ijk} \Omega_j r_k
      -
      u_i
   }{\Delta t},

so that the no-slip / no-penetration conditions are imposed on the particle surface, where :math:`w` is an appropriate weight.
In this project, we use

.. math::

   f^{\prime} \left( d \right)
   =
   \frac{1}{2} \beta \left\{ 1 - \tanh^2 \left( \beta d \right) \right\}

as an approximation of Dirac delta, where :math:`d` is the signed distance from the particle surface normalised by the reference grid size :math:`\Delta \equiv \Delta x \equiv \Delta y`.
One of the indefinite integral of :math:`f^{\prime} \left( d \right)` is

.. math::

   f \left( d \right)
   =
   \frac{1}{2} \left\{ 1 + \tanh \left( \beta d \right) \right\},

which is used in the THINC method (one type of the volume-of-fluid methods) as an phase indicator.

:math:`\beta` controls the sharpness of the surface, the larger :math:`\beta` is, the sharper the surface is, and :math:`\beta \rightarrow \infty` in theory.
In order to enforce the exact surface velocity at the surface location, :math:`\beta` is fixed to :math:`2` so that :math:`f^{\prime} \left( d = 0 \right) = 1`.

*********
Collision
*********

The lubrication force between particles (or particle and the wall) tends to be underestimated especially when the resolution between the objects is insufficient, and as a result particles can penetrate (too much) to each other.
In order to avoid this, this library employs a collision model between objects:

.. math::

   m \der{U_i}{t} + k X_i = 0,

in which a spring whose spring constant is :math:`k` is assumed between objects and a corresponding force pushes away the objects to each other (note that :math:`X_i` is the penetration in the normal direction :math:`i`).
This force is only activated when a penetration is observed, which is the so-called soft-sphere collision model (|TSUJI1993|).

A general solution of the above equation with appropriate boundary conditions

.. math::

   X_i \left( t = 0 \right) = X_i \left( t = T \right) = 0

leads

.. math::

   X_i = A \sin \sqrt{\frac{k}{m}} t,

where

.. math::

   k = m \frac{\pi^2}{T^2}

and thus the repellent force as a function of the penetration depth in the normal direction :math:`x_i` is given by

.. math::

   f_i = k x_i.

Since the particles are rigid, :math:`T = 0` in theory, which is impractical within the current numerical scope (so-called event-driven collision model deals with :math:`T = 0`, which is not easy to handle as well as the fluid motion).

As a remedy, in this project, I consider a timescale

.. math::

   T = \sqrt{\frac{\rho r^3}{\sigma}} \propto r^\frac{3}{2},

where :math:`r` is a particle-based lengthscale (e.g., harmonic average of the colliding particle radii).

Note that this is analogous to the capillary timescale for deformable objects and (again) :math:`\sigma \rightarrow \infty` in theory.
Although giving larger :math:`T` generally stabilise the integration process especially when the suspension volume fraction is large (e.g., :math:`\gt 40 \%`), the results might be unphysical because of large penetration depths.

The above model does not include a dumper, indicating that the collision is perfectly elastic (restitution coefficient is :math:`1`) and no dissipation is involved.
Obviously this is not true in most cases and also tangential elasticity, dumping, and sliding motions should be included (see |COSTA2015| for extensive analyses).

In this project, however, I neglect all these aspects just for simplicity and include only springs in the normal direction to obtain something *plausible*.

