############
Introduction
############

This library is built on `Simple Navier-Stokes Solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_, but I made the following changes.

   * Rigid body motions

      In addition to the working liquid, finite-size particles are also considered and their translational and rotational motions are integrated in time.

   * No temperature coupling

      Temperature field is out of my focus in this project.

   * Uniform grid spacings in all directions

      In addition to the walls, boundary layers are formed on particle surfaces.
      Although their thicknesses are totally different (especially for high Reynolds number flows), using stretched grid has less advantage compared to the original single-phase Rayleigh-Bénard flows.
      Also stretched grid can break the angular momentum balances of the particle motions.
      Thus we adopt uniform grid spacings not only in :math:`y` but in :math:`x` directions.

   * Efficiency

      Since grid points are uniformly distributed in :math:`x` direction, we can use the Discrete Cosine Transform (in the wall-normal direction) instead of the Discrete Fourier Transform (in the homogeneous direction) to solve Poisson equations.
      By doing so, we can halve the number of parallel matrix transform, which is the most computationally demanding procedure.

They are explained in this document, while other parts are omitted since they are identical to the original library.

.. note::

   This library is developed just for fun in my summer holiday.
   Although it is based on `a publicaton <https://www.sciencedirect.com/science/article/pii/S0045793021003716>`_, several major changes are made and thus this library is unofficial.

