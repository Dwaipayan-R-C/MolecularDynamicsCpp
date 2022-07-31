//
// Created by Dwaipayan on 20-05-2022.
//
#include "Eigen/Core"
#include "variable_def.h"
#include "atom_structure.h"
#ifndef YAMD_LJ_KINETIC_AND_POTENTIAL_H
#define YAMD_LJ_KINETIC_AND_POTENTIAL_H

/*
 * Function to calculate the Kinetic energy of the atoms
 * Argument - Atoms
 */
double Kinetic(Atoms &atoms, double rc, double eps, double sigma) ;

/*
 * Function to calculate the Lennard Jones Potential energy of the atoms
 * Argument - Atoms
 */
double Potential(Atoms &atoms, double rc, double eps, double sigma);


#endif //YAMD_LJ_KINETIC_AND_POTENTIAL_H
//#include "Eigen/Core"
//#include "variable_def.h"
//#include "atom_structure.h"
//#ifndef YAMD_LJ_KINETIC_AND_POTENTIAL_H
//#define YAMD_LJ_KINETIC_AND_POTENTIAL_H
//
///*
// * Function to calculate the Kinetic energy of the atoms
// * Argument - Atoms
// */
//double Kinetic(Atoms &atoms) ;
//
///*
// * Function to calculate the Lennard Jones Potential energy of the atoms
// * Argument - Atoms
// */
//double Potential(Atoms &atoms);
//
//
//#endif //YAMD_LJ_KINETIC_AND_POTENTIAL_H