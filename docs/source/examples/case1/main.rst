
.. _example_shear_migration:

.. include:: /references.txt

#########################
Migration in a shear flow
#########################

Comparison with the reference data by |FENG1994| is described, who simulated the lateral migration of a neutrally-buoyant circular object in a two-dimensional wall-bounded shear flow.

*************
Configuration
*************

We consider a square domain (:math:`l_x = l_y = 1`), in which a circular object (whose radius is :math:`0.125`) is positioned at :math:`x = 0.25, y = 0`.
The particle is stationary at :math:`t = 0`, while the liquid (viscosity is adjusted so that the bulk Reynolds number leads :math:`40`) and the wall have a profile of the laminar Couette flow.
Spatial resolutions are varied from :math:`16` to :math:`64` grids per domain size (:math:`4` to :math:`16` grids per diameter) to check the spatial convergence.

In practice, the configuration is specified as follows (:math:`16` grids per domain):

.. literalinclude:: config/exec16.sh
   :language: sh

.. note::
   ``src/fluid/boundary_conditions.c`` and ``src/fluid/init.c`` are modified so that the specified initial conditions (Couette flow profile) are given and the moving walls are treated, respectively.

*******
Results
*******

.. note::
   Conducted by GitHub Actions.

.. literalinclude:: data/16/ci.txt
   :language: text

.. literalinclude:: data/32/ci.txt
   :language: text

.. literalinclude:: data/64/ci.txt
   :language: text

Lateral (wall-normal, :math:`x` direction) migration of a particle as a function of the simulation time :math:`t` is shown below.

.. image:: data/result.pdf
   :width: 800

