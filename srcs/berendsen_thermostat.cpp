//
// Created by Dwaipayan on 26-05-2022.
//
#include <iostream>
#include "../header/berendsen_thermostat.h"

Eigen::Array3Xd berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                                     double relaxation_time, double kb) {

    double vSqdSum, lambda;
    int i, j;
    vSqdSum = 0.;
    double theoretical_kinetic = 1.5*atoms.nb_atoms()*1*Tin;
    for (i = 0; i < nb_atoms; i++) {
        for (j=0; j<3; j++) {
            vSqdSum += atoms.velocities(j,i)*atoms.velocities(j,i);
        }
    }
    lambda = sqrt(1+((3*kb*temperature/(mass*vSqdSum))-1)*timestep/relaxation_time);

    for (i = 0; i < nb_atoms; i++) {
        atoms.velocities.col(i) *= lambda;
    }
    return atoms.velocities;
}
