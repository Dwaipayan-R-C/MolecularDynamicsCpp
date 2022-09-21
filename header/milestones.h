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

#include <iomanip>
#include <iostream>
#include "verlet.h"
#include <functional>
#include "lj_direct_summation.h"
#include "lj_kinetic_and_potential.h"
#include "atom_structure.h"
#include "xyz.h"
#include "berendsen_thermostat.h"
#include "gupta.h"
#include "neighbors.h"
#include <fstream>
#include <cstdlib>
#include <chrono>
#include "mpi.h"
#include "domain.h"

void milestone4(int steps = 2000,double mass=1,double sigma=1, double eps=1, int save_gap=1, double rc =5);
void milestone5(int steps = 2000,double mass=1,double sigma=1, double eps=1, int save_gap=100,  double Tin = 20, double target_temp = 0.6, double boltzmann_kb = 1, double rc = 5);
void milestone6(int steps = 1000,double mass=1,double sigma=1, double eps=1, double boltzmann_kb = 1, double rc = 5);
void milestone7(int steps = 10000, double mass = 196.96657 * 103.6, double delQ = 150, double boltzmann_kb =8.617333262 * pow(10,-5),double  timestep = 2,double rc = 5, int tau = 2000, int save_every = 1000,double sigma=1, double eps=1);

void milestone7_heatCap(int steps = 6000, double mass = 196.96657 * 103.6,  double boltzmann_kb =8.617333262 * pow(10,-5),double  timestep = 2,double rc = 5, int tau = 1000, int save_every = 1000,double sigma=1, double eps=1);
void milestone8(int argc, char *argv[], int steps = 10000 , double mass = 196.96657 * 103.6, double  timestep = 1,double rc =  5, int save_every = 1000, double sigma=1, double eps=1);
void milestone9(int argc, char *argv[] , double mass = 196.96657 * 103.6, double  timestep = 2,double rc =  5, int save_every = 100, bool scale_length = true, int scale_every = 10000, double scale_rate = 0.002 * 141.942, int target_temp = 0);
void lj_potential_distance(int steps = 13000,double mass=1,double sigma=1, double eps=1, int save_gap=10, double rc =5);
void energy_drift(int steps=4000, double mass=196.96657 * 103.6, double boltzmann_kb =8.617333262 * pow(10,-5),double  timestep = 10, double rc=5, double sigma=1, double eps=1);