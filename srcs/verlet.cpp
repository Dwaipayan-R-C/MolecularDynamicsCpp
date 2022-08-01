
#include "Eigen/Core"
#include "../header/verlet.h"
#include "../header/berendsen_thermostat.h"


using namespace Eigen;


void Verlet_one(double dt, Atoms &atoms) {
    int i;
    //  Update positions and velocity with current velocity and acceleration

    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
        atoms.positions.col(i) += atoms.velocities.col(i)*dt;
    }

}


void Verlet_two(double dt, Atoms &atoms) {
    int i;
    //  Update velocity with updated acceleration
    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
    }


}

