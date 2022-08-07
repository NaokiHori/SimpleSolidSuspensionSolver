
.. include:: /references.txt

##############
Implementation
##############

The evolution of momentum from a Runge-Kutta step :math:`k` to :math:`k+1` is given as

.. math::

   \frac{
      u_i^{k+1}
      -
      u_i^{k  }
   }{\Delta t}
   =
   \left( rhs \right)_i
   +
   f_i^{k+\frac{1}{2}},

where :math:`f_i` is the response from the immersed object whose (temporal) formulation will be given later, and

.. math::

   \left( rhs \right)_i
   \equiv
   -
   \gamma^k \dder{p}{x_i}^{k}
   +
   \alpha^k \left( adv + dif \right)_i^{k  }
   +
   \beta^k  \left( adv + dif \right)_i^{k-1}.

We request that, on the boundary, the fluid velocity at :math:`k+1` step :math:`u_i^k` is equal to the velocity of the boundary :math:`u_i^{bnd}`, giving

.. math::

   f_i^{k+\frac{1}{2}}
   \equiv
   \frac{
      u_i^{bnd}
      -
      u_i^{k  }
   }{\Delta t}
   -
   \left( rhs \right)_i.

We can simplify the right-hand-side by employing the prediction velocity :math:`u_i^*`, which is

.. math::

   u_i^{*  }
   =
   u_i^{k  }
   +
   \Delta t \left( rhs \right)_i

and thus the forcing leads

.. math::

   f_i^{k+\frac{1}{2}}
   \equiv
   \frac{
      u_i^{bnd}
      -
      u_i^{*  }
   }{\Delta t}.

See |UHLMANN2005| for details.

The update procedure is as follows.

#. Predict

   .. math::

      \frac{
         u_i^{*  }
         -
         u_i^{k  }
      }{\Delta t}
      =
      \left( rhs \right)_i,

   which is taken care of by :c-lang:`fluid_update_velocity`.

#. Compute forcing

   .. math::

      f_i^{k+\frac{1}{2}}
      \equiv
      \frac{
         u_i^{bnd}
         -
         u_i^{*  }
      }{\Delta t},

   which is taken care of by :c-lang:`suspensions_exchange_momentum`.

#. Update velocity

   .. math::

      \frac{
         u_i^{**}
         -
         u_i^{* }
      }{\Delta t}
      =
      f_i^{k+\frac{1}{2}},

   which is taken care of by :c-lang:`suspensions_update_momentum_fleid`.

#. Enforce continuity

   The last velocity :math:`u_i^{**}` does not satisfy the continuity condition in general and thus we need to enforce it by solving a Poisson equation, which is identical to the single-phase counterpart.

#. Update particles

   Once the fluid variables are updated, particle velocities and positions are updated.In order to stabilise, we update them with a semi-implicit method (Crank-Nicolson scheme) iteratively.

   In particular, we solve

   .. math::

      \delta U_i
      \equiv
      U_i^{k+1}
      -
      U_i^{k  }
      =
      -
      \frac{\gamma \Delta t}{m_p} \int a_i dS
      +
      \frac{1}{m_p} \int u_i dV
      +
      \frac{\gamma \Delta t}{2} \left( F_i^{k+1} + F_i^{k  } \right)

   and

   .. math::

      \delta X_i
      \equiv
      X_i^{k+1}
      -
      X_i^{k  }
      =
      \frac{\gamma \Delta t}{2} \left(
         U_i^{k+1}
         +
         U_i^{k  }
      \right)

   iteratively until the particle locations converge (notice that :math:`F_i^{k+1}` is a function of :math:`X_i^{k+1}`).

   The computations of the right-hand-sides are taken care of by :c-lang:`suspensions_increment_particles`, while the update process (left-hand-sides) are done by :c-lang:`suspensions_update_particles`, respectively.

