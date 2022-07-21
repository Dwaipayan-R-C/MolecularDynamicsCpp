//
// Created by Dwaipayan on 06-05-2022.
//

#include "Eigen/Core"
#include "variable_def.h"
#include "atom_structure.h"


#ifndef YAMD_VERLET_H
#define YAMD_VERLET_H


/*
 * Verlet_one is the prediction step with updated position and velocity by dt/2
 * Argument - timestep, Atoms
 */
void initialize(Atoms &atoms);
tuple<Eigen::Array3Xd,Eigen::Array3Xd>Verlet_one(double dt, Atoms &atoms);


/*
 * Verlet_two is the propagation step with updated velocity by dt
 * Argument - timestep, Atoms
 */
Eigen::Array3Xd Verlet_two(double dt, Atoms &atoms);



#endif //YAMD_VERLET_H
