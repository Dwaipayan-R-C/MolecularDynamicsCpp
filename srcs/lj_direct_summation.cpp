/*
* Copyright 2021 Dwaipayan Roy Chowdhury
*
* ### MIT license
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
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


