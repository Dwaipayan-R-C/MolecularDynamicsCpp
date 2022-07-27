//
// Created by Dwaipayan on 15-05-2022.
//
#include "Eigen/Core"

#ifndef YAMD_ATOM_STRUCTURE_H
#define YAMD_ATOM_STRUCTURE_H
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
//using Names_t = std::array<10>;
//std::string Names_t[10];
using namespace std;

struct Atoms {
//    Names_t atom_name ;
    Masses_t masses;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(Positions_t &p)
            : positions{p},
              velocities{3, p.cols()},
              forces{3, p.cols()},
              masses{p.cols()}{
        velocities.setZero();
        forces.setZero();
    }
    Atoms(Positions_t &p, Velocities_t &v)
            : positions{p},
              velocities{v},
              forces{3, p.cols()},
              masses{p.cols()}{
        velocities.setZero();
        forces.setZero();
    }
    void resize(const int size){
        positions.conservativeResize(3, size);
        velocities.conservativeResize(3, size);
        forces.conservativeResize(3, size);
        masses.conservativeResize(size);
    }
    Atoms(const Positions_t &p, double mass) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setConstant(mass);
    }

    Eigen::Index nb_atoms() const {
        return positions.cols();        }
};
#endif //YAMD_ATOM_STRUCTURE_H
