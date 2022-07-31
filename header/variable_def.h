#ifndef YAMD_VARIABLE_DEF_H
#define YAMD_VARIABLE_DEF_H
#include "Eigen/Core"
#include "atom_structure.h"
#include <iostream>

using namespace Eigen;
// XYZ File constants
const std::string filename = "/trajectory";
const std::string file_extension = ".xyz";

// Simulation Variable Declaration
// const double mass = 1; //  g/mol
const double mass = 196.96657 * 103.6; //  g/mol
//gold = 196.96657
const int nb_steps = 1000;
const int save_every = 100;
const int nb_atoms = 923;
// const int nb_atoms = 923;
//const int nb_atoms = 3871;
const double Tin = 20;
// const double kb = 1; // eV/K
const double kb = 8.617333262 * pow(10,-5); // eV/K
const double rc = 1;
const double delQ = 30; // for 150 in EKin vs Timestep
const int tau = 1000;
const int measurement_gap = tau/2;

// Lennard Jones variables
const double eps = 1;
const double sigma = 1;
const double timescale = 1;  // 1 = 1 fs // for 2 it ran in EKin vs Timestep

//const double timescale = .001 * pow((mass * pow(sigma,2)/eps),(1/2));
#endif //YAMD_VARIABLE_DEF_H
