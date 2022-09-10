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
void milestone9(int argc, char *argv[], int steps = 800000 , double mass = 196.96657 * 103.6, double  timestep = 2,double rc =  5, int save_every = 100, bool scale_length = true, int scale_every = 20000, double scale_rate = 0.005 * 144.25, int target_temp = 0);
void lj_potential_distance(int steps = 13000,double mass=1,double sigma=1, double eps=1, int save_gap=10, double rc =5);
void energy_drift(int steps=4000, double mass=196.96657 * 103.6, double boltzmann_kb =8.617333262 * pow(10,-5),double  timestep = 10, double rc=5, double sigma=1, double eps=1);