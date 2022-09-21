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
#include "../header/berendsen_thermostat.h"

Eigen::Array3Xd berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                                     double relaxation_time, double kb, double Tin) {

    double vSqdSum, lambda;
    int i, j;
    vSqdSum = 0.;    
    for (i = 0; i < atoms.nb_atoms(); i++) {
        for (j=0; j<3; j++) {
            vSqdSum += atoms.velocities(j,i)*atoms.velocities(j,i);
        }
    }
    lambda = sqrt(1+((3*kb*temperature*atoms.nb_atoms()/(atoms.masses(0)*vSqdSum))-1)*timestep/relaxation_time);

    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) *= lambda;
    }
    return atoms.velocities;
}
