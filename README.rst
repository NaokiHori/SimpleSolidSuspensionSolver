########################
Simple Suspension Solver
########################

|License|_ |Documentation|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleSuspensionSolver
.. _License: https://opensource.org/licenses/MIT

.. |Documentation| image:: https://github.com/NaokiHori/SimpleSuspensionSolver/actions/workflows/documentation.yml/badge.svg
.. _Documentation: https://github.com/NaokiHori/SimpleSuspensionSolver/actions/workflows/documentation.yml

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleSuspensionSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleSuspensionSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleSuspensionSolver/blob/main/docs/source/snapshot.png
   :width: 800
   :target: https://youtu.be/EyaXi0o0GZ0

********
Overview
********

This library numerically solves the motion of rigid bodies governed by the Newton-Euler equations suspended in viscous liquid by means of the finite-difference and the immersed boundary methods.

Please refer to `Simple Navier-Stokes Solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_ and `the documentation of this library <https://naokihori.github.io/SimpleSuspensionSolver/index.html>`_ for details.

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

      $ git clone https://github.com/NaokiHori/SimpleSuspensionSolver
      $ cd SimpleSuspensionSolver

#. Build

   .. code-block:: console

      $ docker build -t simplenavierstokessolver:latest .
      $ docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner simplenavierstokessolver:latest
      $ make output
      $ make all

#. Run

   .. code-block:: console

      $ mpirun -n 2 ./a.out

