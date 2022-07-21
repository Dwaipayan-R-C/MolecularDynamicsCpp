//
// Created by Dwaipayan on 26-05-2022.
//
#include "Eigen/Core"
#include "variable_def.h"
#include "atom_structure.h"

#ifndef YAMD_BERENDSEN_THERMOSTAT_H
#define YAMD_BERENDSEN_THERMOSTAT_H

Eigen::Array3Xd berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                                     double relaxation_time, double kb);

#endif //YAMD_BERENDSEN_THERMOSTAT_H
