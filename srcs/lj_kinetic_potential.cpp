//
// Created by Dwaipayan on 20-05-2022.
//

#include "../header/lj_kinetic_and_potential.h"
#include "../header/neighbors.h"

// Function to calculate kinetic energy of the system
//  Function to calculate the kinetic energy of the system
double Kinetic(Atoms &atoms, double rc, double eps, double sigma) {

    double v2, kin;
    kin = 0.;

    for (int i = 0; i < atoms.nb_atoms(); i++) {
        v2 = 0.;

        for (int j = 0; j < 3; j++) {
            v2 += atoms.velocities(j, i) * atoms.velocities(j, i);
        }
        kin += atoms.masses(0) * v2 / 2.;
    }

    return kin;
}


// Function to calculate the potential energy of the system
double Potential(Atoms &atoms, double rc, double eps, double sigma) {
    double quot, r2, rnorm, term1, term2, Pot;
    int i, j, k;

    Pot = 0.;
    NeighborList neighbor_list(rc);
    neighbor_list.update(atoms);
    for (auto [i, j]: neighbor_list) {
        r2 = 0.;
        for (k = 0; k < 3; k++) {
            r2 += (atoms.positions(k, i) - atoms.positions(k, j)) *
                  (atoms.positions(k, i) - atoms.positions(k, j));
        }
        rnorm = sqrt(r2);
        quot = sigma / rnorm;
        term1 = pow(quot, 12.);
        term2 = pow(quot, 6.);
        Pot += 4 * eps * (term1 - term2);
    }
    Pot = Pot / 2;

    return Pot;
}

// Created by Dwaipayan on 20-05-2022.






