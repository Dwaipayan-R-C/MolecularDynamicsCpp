
#include "Eigen/Core"
#include "header/verlet.h"
#include "header/berendsen_thermostat.h"


using namespace Eigen;

void initialize(Atoms &atoms) {
    int n, p, i, j, k;
    double pos;

    // Number of atoms in each direction
    n = int(ceil(pow(nb_atoms, 1.0/3)));

    //  spacing between atoms along a given direction
    pos = 20 / n;

    //  index for number of particles assigned positions
    p = 0;
    //  initialize positions
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++) {
                if (p<nb_atoms) {

                    atoms.positions(0,p) = (i + 0.5)*pos;
                    atoms.positions(1,p) = (j + 0.5)*pos;
                    atoms.positions(2,p) = (k + 0.5)*pos;
                }
                p++;
            }
        }
    }
}
tuple<Array3Xd, Array3Xd> Verlet_one(double dt, Atoms &atoms) {
    int i;
    //  Update positions and velocity with current velocity and acceleration

    for (i = 0; i < nb_atoms; i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
        atoms.positions.col(i) += atoms.velocities.col(i) * dt + 0.5*atoms.forces.col(i)*dt*dt;
    }
//    std::cout << atoms.velocities.sum()<< std::endl;
    return make_tuple(atoms.velocities, atoms.positions);
}


Array3Xd Verlet_two(double dt, Atoms &atoms) {
    int i;
    //  Update velocity with updated acceleration
    for (i = 0; i < nb_atoms; i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
    }
//    std::cout << atoms.velocities.sum()<< std::endl;
//    atoms.velocities = berendsen_thermostat(atoms,10, timescale, 10*timescale);
    return atoms.velocities;
}

