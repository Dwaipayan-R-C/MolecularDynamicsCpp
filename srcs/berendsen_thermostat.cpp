//
// Created by Dwaipayan on 26-05-2022.
//
#include <iostream>
#include "../header/berendsen_thermostat.h"

Eigen::Array3Xd berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                                     double relaxation_time, double kb, double Tin) {

    double vSqdSum, lambda;
    int i, j;
    vSqdSum = 0.;    
    for (i = 0; i < atoms.nb_atoms(); i++) {
        for (j=0; j<3; j++) {
            vSqdSum += atoms.velocities(j,i)*atoms.velocities(j,i);
        }
    }
    lambda = sqrt(1+((3*kb*temperature*atoms.nb_atoms()/(atoms.masses(0)*vSqdSum))-1)*timestep/relaxation_time);

    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) *= lambda;
    }
    return atoms.velocities;
}
