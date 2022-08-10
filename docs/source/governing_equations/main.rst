
.. include:: /references.txt

###################
Governing equations
###################

****************
Momentum balance
****************

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

   \rho_p V_p \oder{U_i}{t}
   =
   -
   \int_{\partial V_p} a_i d S_p^k
   +
   \oder{}{t} \int_{V_p} u_i d V_p
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

   f^{\prime} \left( x \right)
   =
   \frac{1}{2} \beta \left\{ 1 - \tanh^2 \left( \beta x \right) \right\}

as an approximation of Dirac delta, where :math:`x` is the signed distance from the particle surface normalised by the reference grid size :math:`\Delta \equiv \Delta x \equiv \Delta y`.
One of the indefinite integral of :math:`f^{\prime} \left( x \right)` is

.. math::

   f \left( x \right)
   =
   \frac{1}{2} \left\{ 1 + \tanh \left( \beta x \right) \right\},

which is used in the THINC method (one type of the volume-of-fluid methods, see e.g., |XIAO2005|) as an phase indicator.

.. image:: data/indicator.pdf
   :width: 600

:math:`\beta` controls the sharpness of the surface, the larger :math:`\beta` is, the sharper the surface is, and :math:`\beta \rightarrow \infty` in theory.
For the time being, :math:`\beta` is fixed to :math:`2` so that :math:`f^{\prime} \left( d = 0 \right) = 1`, whose effects should be further investigated.

*********
Collision
*********

The lubrication force coming from the hydrodynamic interaction between particles (or a particle and a wall) tends to be underestimated especially when the resolution between the objects is insufficient, and as a result particles can penetrate (too much) to each other.
In order to avoid this, this library employs a collision model between objects, which is based on a simple spring model:

.. math::

   m \oder{U_i}{t} + k X_i = 0,

in which a spring whose spring constant is :math:`k` is assumed between objects and a corresponding force whose direction is parallel to the normal and magnitude is proportional to the penetration depth pushes away the objects to each other.
This force is only activated when a penetration is observed and try to resolve the penetration gently (in :math:`T`), which is the so-called soft-sphere collision model (|TSUJI1993|).

A general solution of the above equation with boundary conditions

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

Since the particles are rigid, :math:`T = 0` in theory, which is impractical within the current numerical scope (so-called event-driven collision model deals with :math:`T = 0`, which is not easy to integrate in time simultaneously as the fluid behaviour).

As a remedy, in this project, I consider a timescale

.. math::

   T = \sqrt{\frac{\rho r^3}{\sigma}} \propto We^{\frac{1}{2}} r^{\frac{3}{2}},

which is analogous to the capillary timescale for deformable objects whose :math:`\sigma` is the surface tension coefficient (and :math:`We` is a pseudo Weber number).
Note that :math:`r` is a particle-based length scale (e.g., harmonic average of the colliding particle radii).

In theory, :math:`\sigma \rightarrow \infty` (or equivalently :math:`T \rightarrow 0`) since the objects are rigid.
Giving larger :math:`T` stabilises the integration process because the repellent force becomes smaller, which is especially useful when the suspension volume fraction is relatively large (e.g., :math:`\gt 40 \%`).
However, the outcome might be unphysical because of the large penetration depths, and thus a proper timescale should be determined by the user eventually.

The above spring model does not include a dumper, indicating that the collision is perfectly elastic (restitution coefficient is :math:`1`) and no dissipation is involved.
Obviously this is not true in most cases and also elastic, dumping, and sliding motions in the tangential direction should be included (see |COSTA2015| for extensive analyses).

In this project, however, I neglect all these aspects just for simplicity, and include only springs in the normal directions to obtain *plausible* results.

