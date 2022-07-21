//
// Created by Dwaipayan on 20-05-2022.
//
#include <iostream>
#include "header/lj_direct_summation.h"
#include "header/neighbors.h"

Array3Xd lj_direct_summation_force(Atoms &atoms){
    int i, j, k;
    double f, rSqd;
    double rij[3];
    Eigen::Vector3d distance_vector;
    Array3Xd accln(3,nb_atoms);
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
            accln.col(i) += (pos.col(i) - pos.col(j))*f/mass;
            accln.col(j) -= (pos.col(i) - pos.col(j))*f/mass;
        }
    }
    return accln;
}
//
// Created by Dwaipayan on 20-05-2022.
//

//#include "header/lj_direct_summation.h"
//
//Array3Xd lj_direct_summation_force(Array3Xd &pos){
//    int i, j, k;
//    double f, rSqd;
//    double rij[3];
//    Eigen::Vector3d distance_vector;
//    Array3Xd accln(3,nb_atoms);
//    accln = 0;
//    for (i = 0; i < nb_atoms-1; i++) {   // loop over all distinct pairs i,j
//        for (j = i+1; j < nb_atoms; j++) {
//            // initialize r^2 to zero
//            rSqd = 0;
//            distance_vector = (pos.col(i) - pos.col(j));
//            rSqd = pow(distance_vector.norm(),2);
//            f = 24*eps*(2*pow(sigma,12)*pow(rSqd, -7) - pow(sigma,6)*pow(rSqd, -4));
////            f = 48*eps*(pow(sigma/distance_vector.norm(),6)*(2*pow(sigma/distance_vector.norm(),6)-0.5))/rSqd;
//            accln.col(i) += (pos.col(i) - pos.col(j))*f*mass;
//            accln.col(j) -= (pos.col(i) - pos.col(j))*f*mass;
//
//        }
//    }
//    return accln;
//}