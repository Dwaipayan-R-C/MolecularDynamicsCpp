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
#include "../header/lj_kinetic_and_potential.h"
#include "../header/neighbors.h"

// Function to calculate kinetic energy of the system
//  Function to calculate the kinetic energy of the system
double Kinetic(Atoms &atoms, double rc, double eps, double sigma) {

    double v2, kin;
    atoms.kin_energy.setZero();
    kin = 0.;

    for (int i = 0; i < atoms.nb_atoms(); i++) {
        v2 = 0.;
        for (int j = 0; j < 3; j++) {
            v2 += atoms.velocities(j, i) * atoms.velocities(j, i);
            atoms.kin_energy(i) += 0.5 * atoms.masses(i)*atoms.velocities(j, i) * atoms.velocities(j, i);
        }        
        kin += atoms.masses(0) * v2 / 2.;        
    }
    
    return atoms.kin_energy.sum();
}


// Function to calculate the potential energy of the system
double Potential(Atoms &atoms, double rc, double eps, double sigma) {
    double quot, r2, rnorm, term1, term2, Pot;
    int i, j, k;

    Pot = 0.;
    NeighborList neighbor_list(rc);
    neighbor_list.update(atoms);
    for (auto [i, j]: neighbor_list) {
        r2 = 0.;
        for (k = 0; k < 3; k++) {
            r2 += (atoms.positions(k, i) - atoms.positions(k, j)) *
                  (atoms.positions(k, i) - atoms.positions(k, j));
        }
        rnorm = sqrt(r2);
        quot = sigma / rnorm;
        term1 = pow(quot, 12.);
        term2 = pow(quot, 6.);
        Pot += 4 * eps * (term1 - term2);
    }
    Pot = Pot / 2;

    return Pot;
}

// Created by Dwaipayan on 20-05-2022.






