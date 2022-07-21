//
// Created by Dwaipayan on 15-05-2022.
//
#include "Eigen/Core"

#ifndef YAMD_ATOM_STRUCTURE_H
#define YAMD_ATOM_STRUCTURE_H
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
//using Names_t = std::array<10>;
//std::string Names_t[10];
using namespace std;

struct Atoms {
//    Names_t atom_name ;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(Positions_t &p)
            : positions{p},
              velocities{3, p.cols()},
              forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();        }
};
#endif //YAMD_ATOM_STRUCTURE_H
