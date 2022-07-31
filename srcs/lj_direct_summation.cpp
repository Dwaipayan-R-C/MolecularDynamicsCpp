//
// Created by Dwaipayan on 20-05-2022.
//
#include <iostream>
#include "../header/lj_direct_summation.h"
#include "../header/neighbors.h"


void lj_direct_summation_force(Atoms &atoms, double rc, double eps, double sigma){
    int i, j, k;
    double f, rSqd;
    double rij[3];
    Eigen::Vector3d distance_vector;
    Array3Xd accln(3,atoms.nb_atoms());
    accln = 0;
    Array3Xd pos = atoms.positions;
    NeighborList neighbor_list(rc);
    neighbor_list.update(atoms);
    for (auto[i, j]: neighbor_list) {
        if (i < j) {
            rSqd = 0;
            distance_vector = (pos.col(i) - pos.col(j));
            rSqd = pow(distance_vector.norm(),2);
            f = 24*eps*(2*pow(sigma,12)*pow(rSqd, -7) - pow(sigma,6)*pow(rSqd, -4));
            accln.col(i) += (pos.col(i) - pos.col(j))*f/atoms.masses(0);
            accln.col(j) -= (pos.col(i) - pos.col(j))*f/atoms.masses(0);
        }
    }
    atoms.forces = accln;
}


