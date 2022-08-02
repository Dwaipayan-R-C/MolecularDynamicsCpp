#include <gtest/gtest.h>
#include "../header/gupta.h"
#include "../header/neighbors.h"


TEST(GuptaTest, Forces) {
    constexpr int nx = 2, ny = 2, nz = 2;
    constexpr double lattice_constant = 1.5;
    constexpr double cutoff = 5.0;
    constexpr double delta = 0.0001;  // difference used for numerical (finite difference) computation of forces

    double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
    
    NeighborList neighbor_list(cutoff);
    Positions_t positions(3, nx * ny * nz);
    Atoms atoms(positions);
    atoms.masses = 1;
    // we create a cubic lattice with random displacements
    atoms.positions.setRandom();  // random numbers between -1 and 1
    atoms.positions *= 0.1;
    for (int x{0}, i{0}; x < nx; ++x) {
        for (int y{0}; y < ny; ++y) {
            for (int z{0}; z < nz; ++z, ++i) {
                atoms.positions(0, i) += x * lattice_constant;
                atoms.positions(1, i) += y * lattice_constant;
                atoms.positions(2, i) += z * lattice_constant;
            }
        }
    }

    neighbor_list.update(atoms);
    atoms.forces.setZero();
    double e0{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences approximation
    for (int i{0}; i < atoms.nb_atoms(); ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            neighbor_list.update(atoms);
            double eplus{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            neighbor_list.update(atoms);
            double eminus{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            // check whether finite-difference and analytic forces agree
            EXPECT_NEAR(fd_force, forces0(j, i), 1e-5);
        }
    }
}