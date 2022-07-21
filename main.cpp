#include <iomanip>
#include <iostream>
#include "header/verlet.h"
#include <functional>
#include "header/lj_direct_summation.h"
#include "header/lj_kinetic_and_potential.h"
#include "header/atom_structure.h"
#include "header/xyz.h"
#include "header/berendsen_thermostat.h"
#include "header/gupta.h"
#include "header/neighbors.h"


using namespace Eigen;

//check starts
//constexpr int nx = 2, ny = 2, nz = 2;
//constexpr double lattice_constant = 1.5;
//constexpr double cutoff = 5.0;
//constexpr double delta = 0.0001;  // difference used for numerical (finite difference) computation of forces
//
//double A = 0.2061, xi = 1.790, p = 10.229, q = 4.036, re = 4.079/sqrt(2);
//Positions_t positions(3,nx * ny * nz);
//NeighborList neighbor_list(cutoff);
//
//Atoms atoms(positions);
// check ends
// Global Variable Declaration
string names;
double potential = 0;
double kinetic_energy = 0;
double total_energy = 0;
int count_relax = 0;
int relax_value = 0;
Positions_t positions(3,nb_atoms);
Atoms atoms = {
        positions
};
double calculateAverage(double numbers[], int count )
{
    double sum = 0;         //sum is used to add all the values in the array
    double average;

    for (int i= 0 ; i < count; i++)
        sum += numbers[i];

    average = sum /count;
    return average;
}

int main() {

    // reads the initial atomic position and velocity from xyz file
    // assigns to the corresponding
    clock_t start, end;
    start = clock();
    auto [positions, velocities]{read_xyz_with_velocities("../../xyz/cluster_923.xyz")};

    atoms.positions=positions;
    atoms.velocities = 0;

    // gets the initial force from Lennard Jones potential
//    atoms.forces = lj_direct_summation_force(atoms);
    int number = 0;
    double temperature = 0;
    double totalEnergy = 0;
    int k = 0;
    double alpha = 0;
    for (int i=0; i<nb_steps;i++){
        // Verlet prediction 1
        tuple<Array3Xd, Array3Xd> verlet_one_return_data = Verlet_one(timescale, atoms);
        atoms.velocities = get<0>(verlet_one_return_data);
        atoms.positions = get<1>(verlet_one_return_data);

        //Increase the timestep
        //implement the formula of constant latent heat
        // Onserve the energy diagrams

            if(i%tau == 0 && i != 0) {
                // avg_temperature[k] = temperature/20;
//            std::cout << temperature/tau << std::endl;
                std::cout << "[" << totalEnergy / (relax_value) << " , " << temperature / (relax_value) << "]," << std::endl;
                relax_value = 0;
                count_relax = 0;
                temperature = 0;
                totalEnergy = 0;
                // Determine alpha
                // del Q = Ekinetic' - Ekinetic
                // del Q = Ekinetic (sq of alpha - 1)
                alpha = sqrt((delQ / kinetic_energy) + 1);
                atoms.velocities = atoms.velocities * alpha;
                k = k + 1;
            }


        // Update force for Lennard jones
        // atoms.forces = lj_direct_summation_force(atoms);
        // GUPTA
        NeighborList neighbor_list(rc);
        neighbor_list.update(atoms);
        potential = gupta(atoms,neighbor_list,rc);

        // Verlet 2 propagation with updated force
        atoms.velocities = Verlet_two(timescale, atoms);
        // atoms.velocities = berendsen_thermostat(atoms,1500000, timescale, 10*timescale, kb);

        // Calculate Kinetic and Lennard Jones Potential Energy
        // potential = Potential(atoms);
        kinetic_energy = Kinetic(atoms);
        total_energy =  potential + kinetic_energy;
//        std::cout << total_energy<< std::endl;
        // save xyz file
        if(i%save_every == 0){
            write_xyz(filename+ to_string(number) + file_extension, atoms);
            number = number+1;
//            std::cout << temperature/tau << std::endl;
//            std::cout << totalEnergy/tau << std::endl;
        }
//        std::cout << total_energy << std::endl;
        if(count_relax>=measurement_gap){
            temperature += kinetic_energy * 2/(3*kb* nb_atoms);
            totalEnergy += total_energy;
            relax_value+=1;
//            count_relax = 0;
        }
        count_relax +=1;

//        totalEnergy += kinetic_energy;


//        std::cout <<"[ "<< total_energy<<" , "<<i*timescale <<"],"<< std::endl;
//        std::cout << kinetic_energy << std::endl;
    }

    end = clock();
    printf ("Time taken: %f secs\n",((float) end - start)/CLOCKS_PER_SEC);
    return 0;
}
//
//
//#include "header/atom_structure.h"
//#include "header/neighbors.h"
//#include "header/xyz.h"
//
//#include <gtest/gtest.h>
//
//
////TEST(GuptaTest, Forces) {
////    // we create a cubic lattice with random displacements
////    atoms.positions.setRandom();  // random numbers between -1 and 1
////    atoms.positions *= 0.1;
////    for (int x{0}, i{0}; x < nx; ++x) {
////        for (int y{0}; y < ny; ++y) {
////            for (int z{0}; z < nz; ++z, ++i) {
////                atoms.positions(0, i) += x * lattice_constant;
////                atoms.positions(1, i) += y * lattice_constant;
////                atoms.positions(2, i) += z * lattice_constant;
////            }
////        }
////    }
////
////    neighbor_list.update(atoms);
////    atoms.forces.setZero();
////    double e0{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
////    Forces_t forces0{atoms.forces};
////
////    // loop over all atoms and compute forces from a finite differences approximation
////    for (int i{0}; i < atoms.nb_atoms(); ++i) {
////        // loop over all Cartesian directions
////        for (int j{0}; j < 3; ++j) {
////            // move atom to the right
////            atoms.positions(j, i) += delta;
////            neighbor_list.update(atoms);
////            double eplus{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
////            // move atom to the left
////            atoms.positions(j, i) -= 2 * delta;
////            neighbor_list.update(atoms);
////            double eminus{gupta(atoms, neighbor_list, cutoff, A, xi, p, q, re)};
////            // move atom back to original position
////            atoms.positions(j, i) += delta;
////
////            // finite-differences forces
////            double fd_force{-(eplus - eminus) / (2 * delta)};
////
////            // check whether finite-difference and analytic forces agree
////            EXPECT_NEAR(fd_force, forces0(j, i), 1e-5);
////        }
////    }
//////    Names_t names{{"H", "H", "H", "H"}};
////    Positions_t positions(3, 4);
////    positions << 0, 1, 0, 0,
////            0, 0, 1, -1,
////            0, 0, 0, 0;
////
////    Atoms atoms(positions);
////    NeighborList neighbor_list(1.5);
////    auto &[seed, neighbors]{neighbor_list.update(atoms)};
////
////    // All atoms except 3 and 4 are neighbors of each other
////    EXPECT_EQ(neighbor_list.nb_neighbors(), 10);
////
////    EXPECT_EQ(neighbor_list.nb_neighbors(0), 3);
////    EXPECT_EQ(neighbor_list.nb_neighbors(1), 3);
////    EXPECT_EQ(neighbor_list.nb_neighbors(2), 2);
////    EXPECT_EQ(neighbor_list.nb_neighbors(3), 2);
////
////    EXPECT_TRUE((neighbors(Eigen::seq(seed(0), seed(1) - 1)) == Eigen::Array3i{3, 1, 2}).all());
////    EXPECT_TRUE((neighbors(Eigen::seq(seed(1), seed(2) - 1)) == Eigen::Array3i{3, 0, 2}).all());
////    EXPECT_TRUE((neighbors(Eigen::seq(seed(2), seed(3) - 1)) == Eigen::Array2i{0, 1}).all());
////    EXPECT_TRUE((neighbors(Eigen::seq(seed(3), seed(4) - 1)) == Eigen::Array2i{0, 1}).all());
////}
//
//// Modify the functions in the test
//// Check the force comp - maybe it's less.
//// Use Eigen.position.col subtraction and then .norm() to calculate the difference
//
//
//
//
///*
//* Copyright 2021 Lars Pastewka
//*
//* ### MIT license
//*
//* Permission is hereby granted, free of charge, to any person obtaining a copy
//* of this software and associated documentation files (the "Software"), to deal
//* in the Software without restriction, including without limitation the rights
//* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//* copies of the Software, and to permit persons to whom the Software is
//* furnished to do so, subject to the following conditions:
//*
//* The above copyright notice and this permission notice shall be included in
//* all copies or substantial portions of the Software.
//*
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//* SOFTWARE.
//*/
//
////#include <gtest/gtest.h>
////#include "header/atom_structure.h"
////#include "header/lj_direct_summation.h"
////#include "header/lj_kinetic_and_potential.h"
////
////TEST(LJDirectSummationTest, Forces) {
////    const size_t nb_atoms = 10;
////    constexpr double epsilon = 1;  // choose different to 1 to pick up missing factors
////    constexpr double sigma = .7;
////    constexpr double delta = 0.0001;  // difference used for numerical (finite difference) computation of forces
////    string names;
////    double potential;
////    double kinetic_energy;
////    double total_energy;
////    Positions_t positions(3, nb_atoms);
////    Atoms atoms = {
////            positions
////    };
//////    atoms.positions = 0;
//////    atoms.positions(0,0) = 0;
//////    atoms.positions(0,1) = 0.5;
//////    std::cout << atoms.positions <<std::endl;
////    atoms.positions.setRandom();  // random numbers between -1 and 1
////
////// compute and store energy of the indisturbed configuration
////    double e0 = Potential(atoms);
////    Array3Xd forces0 = lj_direct_summation_force(atoms.positions);
////    std::cout << forces0 <<std::endl;
////// loop over all atoms and compute forces from a finite differences approximation
////    for (int i{0}; i < nb_atoms; ++i) {
////// loop over all Cartesian directions
////        for (int j{0}; j < 3; ++j) {
////// move atom to the right
////            atoms.positions(j, i) += delta;
////            double eplus{Potential(atoms)};
////// move atom to the left
////            atoms.positions(j, i) -= 2 * delta;
////            double eminus{Potential(atoms)};
////// move atom back to original position
////            atoms.positions(j, i) += delta;
////
////// finite-differences forces
////            double fd_force{-(eplus - eminus) / (2 * delta)};
////            std::cout << fd_force <<std::endl;
//////            std::cout << abs(fd_force - forces0(j, i))  <<std::endl;
//////            std::cout << forces0(j, i) <<std::endl;
//////            std::cout << "\n" <<std::endl;
////// check whether finite-difference and analytic forces agree
////            if (abs(forces0(j, i)) > 1e-10) {
////                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / abs(forces0(j, i)), 0, 1e-5);
////            } else {
////                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
////            }
////        }
////    }
////}
