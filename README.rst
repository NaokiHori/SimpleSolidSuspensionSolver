##############################
Simple Solid Suspension Solver
##############################

For circular (2D) or spherical (3D) objects, please refer to `Simple IBM Solver <https://github.com/NaokiHori/SimpleIBMSolver>`_.

|License|_ |CI|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleSolidSuspensionSolver
.. _License: https://opensource.org/licenses/MIT

.. |CI| image:: https://github.com/NaokiHori/SimpleSolidSuspensionSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleSolidSuspensionSolver/actions/workflows/ci.yml

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleSolidSuspensionSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleSolidSuspensionSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleSolidSuspensionSolver/blob/main/docs/source/snapshot.png
   :width: 800
   :target: https://youtu.be/iuO5CxvAlio

********
Overview
********

This library numerically solves the motion of rigid bodies (ellipses) governed by the Newton-Euler equations suspended in viscous liquid by means of the finite-difference and the immersed boundary methods.

Please refer to `Simple Navier-Stokes Solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_ and `the documentation of this library <https://naokihori.github.io/SimpleSolidSuspensionSolver/index.html>`_ for details.

This library is still under development, and especially the documentation is far from being sufficient.

********
Features
********

* MPI-parallelised
* Eulerian-based (no Lagrangian points) IBM
* Stable collision model between ellipses

**********
Dependency
**********

`Docker <https://www.docker.com>`_

***********
Quick start
***********

#. Create working directory

   .. code-block:: console

      $ mkdir /path/to/your/working/directory
      $ cd    /path/to/your/working/directory

#. Fetch source

   .. code-block:: console

      $ git clone https://github.com/NaokiHori/SimpleSolidSuspensionSolver
      $ cd SimpleSolidSuspensionSolver

#. Build

   .. code-block:: console

      $ docker build -t simplenavierstokessolver:latest .
      $ docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner simplenavierstokessolver:latest
      $ make output
      $ make all

#. Run

   .. code-block:: console

      $ mpirun -n 2 ./a.out

********
Examples
********

Several examples can be found in the documentation.

#. `Migration of a circular object in a shear flow <https://naokihori.github.io/SimpleSolidSuspensionSolver/examples/case1/main.html>`_

#. `Segr??-Silberberg effect <https://naokihori.github.io/SimpleSolidSuspensionSolver/examples/case2/main.html>`_

#. `Rotation of an ellipse in a shear flow <https://naokihori.github.io/SimpleSolidSuspensionSolver/examples/case3/main.html>`_

#. `Suspension in a plane Poiseuille flow <https://naokihori.github.io/SimpleSolidSuspensionSolver/examples/case4/main.html>`_

