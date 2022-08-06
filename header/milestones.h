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

void milestone4(int steps = 2000,double mass=1,double sigma=1, double eps=1, int save_gap=100, double rc =5);
void milestone5(int steps = 2000,double mass=1,double sigma=1, double eps=1, int save_gap=100,  double Tin = 20, double target_temp = 0.6, double boltzmann_kb = 1, double rc = 5);
void milestone6(int steps = 1000,double mass=1,double sigma=1, double eps=1, double boltzmann_kb = 1, double rc = 5);
void milestone7(int steps = 100000, double mass = 196.96657 * 103.6, double delQ = 0.5, double boltzmann_kb =8.617333262 * pow(10,-5),double  timestep = 2,double rc = 5, int tau = 500, int save_every = 500,double sigma=1, double eps=1);
void milestone8(int argc, char *argv[], int steps = 50000 , double mass = 196.96657 * 103.6, double  timestep = 1,double rc =  5, int save_every = 1000, double sigma=1, double eps=1);