//
// Created by Dwaipayan on 20-05-2022.
//

#include "Eigen/Core"
#include "variable_def.h"
#include "atom_structure.h"
#include <iostream>

using namespace Eigen;
#ifndef YAMD_LJ_DIRECT_SUMMATION_H
#define YAMD_LJ_DIRECT_SUMMATION_H

Array3Xd lj_direct_summation_force(Atoms &atoms);
#endif //YAMD_LJ_DIRECT_SUMMATION_H
//#include "Eigen/Core"
//#include "variable_def.h"
//#include <iostream>
//
//using namespace Eigen;
//#ifndef YAMD_LJ_DIRECT_SUMMATION_H
//#define YAMD_LJ_DIRECT_SUMMATION_H
//
//Array3Xd lj_direct_summation_force(Array3Xd &pos);
//#endif //YAMD_LJ_DIRECT_SUMMATION_H
