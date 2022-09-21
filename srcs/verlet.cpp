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
#include "Eigen/Core"
#include "../header/verlet.h"
#include "../header/berendsen_thermostat.h"


using namespace Eigen;


void Verlet_one(double dt, Atoms &atoms) {
    int i;
    //  Update positions and velocity with current velocity and acceleration

    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
        atoms.positions.col(i) += atoms.velocities.col(i)*dt;
    }

}


void Verlet_two(double dt, Atoms &atoms) {
    int i;
    //  Update velocity with updated acceleration
    for (i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * dt;
    }


}

