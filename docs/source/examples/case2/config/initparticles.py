## initialise particles (density, radius, positions, velocities)

import sys
import numpy as np


def initialise(n_particles):
    ## arrays to be returned
    # densities
    dens = [1.]
    # major axis
    as_ = [0.125]
    # minor axis
    bs_ = [0.125]
    # x locations
    xs = [0.4]
    # y locations
    ys = [0.]
    # z angles
    azs = [0.]
    # x velocities
    uxs = [0.]
    # y velocities
    uys = [0.]
    # z angular velocities
    vzs = [0.]
    return dens, as_, bs_, xs, ys, azs, uxs, uys, vzs

def output(root, n_particles, dens, as_, bs_, xs, ys, azs, uxs, uys, vzs):
    n_particles = np.array(n_particles, dtype="int32")
    dens = np.array(dens, dtype="float64")
    as_  = np.array(as_,  dtype="float64")
    bs_  = np.array(bs_,  dtype="float64")
    xs   = np.array(xs,   dtype="float64")
    ys   = np.array(ys,   dtype="float64")
    azs  = np.array(azs,  dtype="float64")
    uxs  = np.array(uxs,  dtype="float64")
    uys  = np.array(uys,  dtype="float64")
    vzs  = np.array(vzs,  dtype="float64")
    np.save(f"{root}/n_particles", n_particles)
    np.save(f"{root}/particle_dens", dens)
    np.save(f"{root}/particle_as",    as_)
    np.save(f"{root}/particle_bs",    bs_)
    np.save(f"{root}/particle_xs",     xs)
    np.save(f"{root}/particle_ys",     ys)
    np.save(f"{root}/particle_azs",   azs)
    np.save(f"{root}/particle_uxs",   uxs)
    np.save(f"{root}/particle_uys",   uys)
    np.save(f"{root}/particle_vzs",   vzs)

if __name__ == "__main__":
    n_particles = 1
    dens, as_, bs_, xs, ys, azs, uxs, uys, vzs = initialise(n_particles)
    root = sys.argv[1]
    output(root, n_particles, dens, as_, bs_, xs, ys, azs, uxs, uys, vzs)

