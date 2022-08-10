
.. include:: /references.txt

##############
Implementation
##############

********************
Temporal integration
********************

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

***************
Parallelisation
***************

Since this library is built on `a MPI-parallelised DNS solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_, Eulerian domain is decomposed in :math:`y` direction.
Regarding particles, update procedures such as momentum exchange are conducted on the same Eulerian domain and thus the parallelisation is straightforward.
The velocity and position of the gravity centers are stored in a Lagrangian way and shared among all processes.
Since they are :math:`\mathcal{O} \left( N_p \right)` operations, the computational cost is negligibly small compared to the fluid solver :math:`\mathcal{O} \left( N_x \times N_y \right)`: proportional to the number of grid points.

Collisions (judgements and computations of repellent forces), on the other hand, require relatively heavy computational loads since we need :math:`\mathcal{O} \left( N_p^2 \right)` operations.
Although it is still much smaller than :math:`N_x \times N_y`, it can be fairly heavy to be completed by a main process when the volume fraction or the domain size is large.
Thus this library parallelises the procedure simply as the same way as the decomposing the Eulerian domain.
Namely, there are following collision pairs to be considered:

.. csv-table::

   :file: collision-table.csv

Here the left-most column :math:`n_0` and the top-most row :math:`n_1` show indices of the particles.
The other numbers show the indices of pairs, which are separated to each processor and considered independently.

In practice, we need to translate the index of collision pair :c-lang:`n` to the corresponding particle indices :math:`n_0` and :math:`n_1`, which are handled by :c-lang:`get_particle_indices`.

