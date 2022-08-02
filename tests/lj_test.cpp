#include <gtest/gtest.h>

#include "../header/atom_structure.h"
#include "../header/lj_direct_summation.h"
#include "../header/lj_kinetic_and_potential.h"

TEST(LJDirectSummationTest, Forces) {
    const size_t nb_atoms = 10;
    constexpr double epsilon = 0.7; // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.00001; // difference used for numerical (finite difference) computation of forces
    string names;
    double potential = 0;
    
    // Atoms atoms(nb_atoms);
    Positions_t positions(3, nb_atoms);
    Atoms atoms = {positions};  
    atoms.masses = 1;  
    atoms.positions.setRandom(); // random numbers between -1 and 1
    // atoms.forces.setZero();
    // compute and store energy of the indisturbed configuration
    double e0 = Potential(atoms,10,epsilon,sigma);
    lj_direct_summation_force(atoms, 10,epsilon,sigma);
    Forces_t forces0{atoms.forces};
    
    // loop over all atoms and compute forces from a finite differences
    // approximation
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{Potential(atoms,10,epsilon,sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{Potential(atoms,10,epsilon,sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};
            
            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);

            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}



