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

#include "../header/milestones.h"
#include <list>
using namespace std::chrono;
using namespace Eigen;
using std::ofstream;

/*
* This function is written for milestone 4. Here we check energy conservation.
*/
void milestone4(int steps, double mass, double sigma, double eps, int save_gap,
                double rc) {
    double potential = 0;
    double kinetic_energy = 0;
    double total_energy = 0;
    int number = 0;

    std::ofstream outdata("../data/milestone4.dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{read_xyz_with_velocities("../xyz/lj54.xyz")};
    Atoms atoms{positions, velocities};
    atoms.masses = 1;
    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));
    std::cout << "[ Kinetic , Potential , Total ]" << std::endl;

    // loop starts for the simulation with time
    for (int i = 0; i < steps; i++) {
        // Calculates verlet 1
        Verlet_one(timestep, atoms);
        // Calculates Lennard jones potential
        lj_direct_summation_force(atoms, rc, eps, sigma);
        // Calculates Verlet 2
        Verlet_two(timestep, atoms);
        // Calculates the potential
        potential = Potential(atoms, rc, eps, sigma);
        // Calculates the kinetic energy
        kinetic_energy = Kinetic(atoms, rc, eps, sigma);
        double total_energy = kinetic_energy + potential;
        outdata << "[ " << kinetic_energy << " , " << potential << " , "
                << total_energy << " ]," << std::endl;

        /*
            Outputs the data file and xyz file to the output folder. 
        */
        if (i % save_gap == 0) {
            std::cout << "[ " << kinetic_energy << " , " << potential << " , "
                      << total_energy << " ]," << std::endl;
            write_xyz("../xyz_output/milestone4/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }
    }

    outdata.close();
    std::cout << "Ran Velocity verlet successfully " << std::endl;
    std::cout << "Output has been written to data/milestone4.dat in the format "
                 "- Kinetic Energy, Potential Energy, Total Energy"
              << std::endl;
}

/*
* This function is written for milestone 5. Here we check berendsen thermostat.
*/
void milestone5(int steps, double mass, double sigma, double eps, int save_gap,
                double Tin, double target_temp, double boltzmann_kb,
                double rc) {
    double potential = 0;
    double kinetic_energy = 0;
    double total_energy = 0;
    double temperature = 0;
    double totalEnergy = 0;

    std::ofstream outdata("../data/milestone5.dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{read_xyz_with_velocities("../xyz/lj54.xyz")};
    Atoms atoms{positions, velocities};
    atoms.masses = mass;
    int number = 0;
    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));
    std::cout << "Current temperature , target temperature" << std::endl;
    for (int i = 0; i < steps; i++) {
        Verlet_one(timestep, atoms);
        lj_direct_summation_force(atoms, rc, eps, sigma);
        Verlet_two(timestep, atoms);
        if (i % 10 == 0) {
            atoms.velocities = berendsen_thermostat(
                atoms, target_temp, timestep, 10 * timestep, boltzmann_kb, Tin);
        }
        potential = Potential(atoms, rc, eps, sigma);
        kinetic_energy = Kinetic(atoms, rc, eps, sigma);

        double current_temp = kinetic_energy * 2 / (3 * boltzmann_kb);
        double total_energy = kinetic_energy + potential;
        outdata << "[ " << current_temp << " ," << target_temp << "],"
                << std::endl;

        // save xyz file
        if (i % save_gap == 0) {
            std::cout << "[ " << current_temp << " ," << target_temp << "],"
                      << std::endl;
            write_xyz("../xyz_output/milestone5/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }
    }

    outdata.close();
    std::cout << "Ran Berdendsen Thermostat successfully and output has been "
                 "written to data/milestone5.dat in the format - Total Energy "
                 ", Current temperature, Target temperature"
              << std::endl;
}

/*
* This function is written for milestone 6. Here we see  time variation with
* cluster size and cut-off radius.
*/
void milestone6(int steps, double mass, double sigma, double eps,
                double boltzmann_kb, double rc) {
    int number = 0;
    int count_relax = 0;
    int relax_value = 0;
    double potential = 0;
    double kinetic_energy = 0;
    double total_energy = 0;
    double temperature = 0;
    double totalEnergy = 0;

    std::ofstream outdata("../data/milestone6.dat");

    std::list<int> rc_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions,
          velocities]{read_xyz_with_velocities("../xyz/cluster_923.xyz")};
    Atoms atoms{positions, velocities};
    atoms.masses = mass;

    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));

    for (int x : rc_list) {
        auto start = high_resolution_clock::now();
        NeighborList neighbor_list(x);
        for (int i = 0; i < steps; i++) {
            Verlet_one(timestep, atoms);
            neighbor_list.update(atoms);
            lj_direct_summation_force(atoms, rc, eps, sigma);
            Verlet_two(timestep, atoms);
            potential = Potential(atoms, rc, eps, sigma);
            kinetic_energy = Kinetic(atoms, rc, eps, sigma);
            double total_energy = kinetic_energy + potential;
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        // in milisecond
        outdata << "[ " << x << " ," << duration.count() / 1000 << " ],"
                << std::endl;

        std::cout << "Running" << std::endl;
    }

    outdata.close();
    std::cout << "Ran cutoff radius successfully and output has been written "
                 "to data/milestone6.dat in the format - radius, time"
              << std::endl;
}

/*
* This function is written for milestone 7. Here we perform energy analysis and several
* simulations for kinetic energy and potential energy with time, total energy with temperature.
* We also extract the heat capacity, melting point and latent heat from energy vs temperature
* plot for different clusters and verify the property type i.e. extrinsic/intrinsic
*/
void milestone7(int steps, double mass, double delQ, double boltzmann_kb,
                double timestep, double rc, int tau, int save_every,
                double sigma, double eps) {
    std::string cluster_name = "923";
    int measurement_gap = tau / 2;
    int number = 0;
    int count_relax = 0;
    int relax_value = 0;
    double potential = 0;
    double kinetic_energy = 0;
    double total_energy = 0;
    double temperature = 0;
    double totalEnergy = 0;
    int k = 0;
    double alpha = 0;

    std::ofstream outdata("../data/milestone_7/data_" + cluster_name + ".dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{
        read_xyz_with_velocities("../xyz/custom_" + cluster_name + ".xyz")};
    Atoms atoms{positions};
    atoms.velocities = 0;
    atoms.masses = mass;

    auto start = high_resolution_clock::now();

    NeighborList neighbor_list(rc);
    for (int i = 0; i < steps; i++) {
        Verlet_one(timestep, atoms);
        neighbor_list.update(atoms);
        potential = gupta(atoms, neighbor_list, rc);
        // potential = atoms.energies.sum();
        Verlet_two(timestep, atoms);
        kinetic_energy = Kinetic(atoms, rc, eps, sigma);
        double total_energy = kinetic_energy + potential;

        // save xyz file
        if (i % save_every == 0) {
            write_xyz("../xyz_output/milestone7/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }

        if (count_relax >= measurement_gap) {
            temperature +=
                kinetic_energy * 2 / (3 * boltzmann_kb * atoms.nb_atoms());
            totalEnergy += total_energy;
            relax_value += 1;
        }
        // in milisecond
        outdata << "[ " << kinetic_energy << " ," << potential << " ,"
                << timestep * i << " ]," << std::endl;
        if (i % tau == 0 && i != 0) {

            std::cout << "[" << totalEnergy / (relax_value) << " , "
                      << temperature / (relax_value) << "]," << std::endl;
            // // in milisecond
            // outdata << "[ " << totalEnergy / (relax_value) << " ,"
            //         << temperature / (relax_value) << " ]," << std::endl;
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
        count_relax += 1;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    outdata.close();

    std::cout << "Embedded Atom potential ran successfully " << std::endl;
    std::cout << "Output has been written to data/milestone7.dat in the format "
                 "-  Average Total Energy (eV) vs Average temperature"
              << std::endl;
}

/*
* This function is written for milestone 8. Here we perform the MPI simulation and check 
* if our energy is conserved after implementing the parallel computing with 8 cores and 
* executing the periodic boundary conditions.
*/
void milestone8(int argc, char *argv[], int steps, double mass, double timestep,
                double rc, int save_every, double sigma, double eps) {
    MPI_Init(&argc, &argv);
    int number = 0;
    double potential = 0;
    double kinetic_energy = 0;
    int k = 0;
    double alpha = 0;

    Domain domain(MPI_COMM_WORLD, {30, 30, 30},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});

    std::ofstream outdata("../data/milestone8.dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    auto [positions,
          velocities]{read_xyz_with_velocities("../xyz/cluster_923.xyz")};
    Atoms atoms{positions};
    atoms.velocities = 0;
    atoms.masses = mass;
    atoms.energies = 0;
    atoms.kin_energy = 0;

    auto start = high_resolution_clock::now();
    NeighborList neighbor_list(rc);
    domain.enable(atoms); // divides into subdomains
    for (int i = 0; i < steps; i++) {
        Verlet_one(timestep, atoms);
        // Neighbour list looks for ghost nodes which is why we exchange
        // atoms and updates ghost nodes.
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * rc);
        neighbor_list.update(atoms);
        potential = gupta(atoms, neighbor_list, rc);
        Verlet_two(timestep, atoms);
        kinetic_energy = Kinetic(atoms, rc, eps, sigma);
        double global_potential = 0;
        double global_kinetic = 0;
        double local_kinetic = 0;
        double local_potential = 0;
        for (int k = 0; k < domain.nb_local(); k++) {
            local_potential += atoms.energies(k);
            local_kinetic += atoms.kin_energy(k);
        }

        MPI_Reduce(&local_potential, &global_potential, 1, MPI_DOUBLE, MPI_SUM,
                   0, MPI_COMM_WORLD);

        MPI_Reduce(&local_kinetic, &global_kinetic, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        domain.disable(atoms);
        if (domain.rank() == 0) {
            double total_energy = global_potential + global_kinetic;
            outdata << "[" << global_potential << " , " << global_kinetic
                    << " , " << total_energy << "]," << std::endl;
        }
        if (i % save_every == 0 && domain.rank() == 0) {
            double total_energy = global_potential + global_kinetic;

            std::cout << "[" << global_potential << " , " << global_kinetic
                      << " , " << total_energy << "]," << std::endl;
            write_xyz("../xyz_output/milestone8/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }
        domain.enable(atoms);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    outdata.close();
    MPI_Finalize();
    std::cout << "MPI ran successfully " << std::endl;
    std::cout << "Output has been written to data/milestone8.dat in the format "
                 "-  Average Total Energy (eV) vs Average temperature"
              << std::endl;
}

/*
* This function is written for milestone 9. Here we take two nanowires and study the 
* deformation at varying strain rates. We also present a stress strain curve to show 
* the elastic range and if this matches with the expected results
*/
void milestone9(int argc, char *argv[], double mass, double timestep,
                double rc, int save_every, bool scale_length, int scale_every,
                double scale_rate, int target_temp) {

    MPI_Init(&argc, &argv);
    int number = 0;
    int k = 0;
    Eigen::Array3d scale;
    double domain_size = 141.942;
    double z_scale = domain_size;
    double strain = 0;

    Domain domain(MPI_COMM_WORLD, {50, 50, 141.942},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});
    // 50, 50, 144.25 142.086 141.942
    std::ofstream outdata("../data/milestone9_100_002_small.dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{
        read_xyz_with_velocities("../xyz/whiskers/whisker_small_002.xyz")};
    Atoms atoms{positions};
    atoms.velocities = 0;
    atoms.masses = mass;
    atoms.energies = 0;
    atoms.kin_energy = 0;
    outdata << "step"
            << ", "
            << "strain"
            << ", "
            << "right"
            << ", "
            << "left"
            << ", \n";
    auto start = high_resolution_clock::now();
    NeighborList neighbor_list(rc);
    domain.enable(atoms); // divides into subdomains
    int i=0;
    while (z_scale<=154.433) {
        if (i >= 20000 && scale_length == true && i % scale_every == 0) {
            z_scale += scale_rate;
            scale << 50, 50, z_scale;
            domain.scale(atoms, scale);
        }

        Verlet_one(timestep, atoms);
        // Neighbour list looks for ghost nodes which is why we exchange
        // atoms and updates ghost nodes.
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * rc);
        neighbor_list.update(atoms);
        gupta(atoms, neighbor_list, rc);
        Verlet_two(timestep, atoms);
        double kinetic_energy = Kinetic(atoms, rc, 1, 1);
        double global_kinetic = 0;
        double local_kinetic = 0;

        berendsen_thermostat(atoms, target_temp, timestep, 1000, 8.617333262e-5,
                             0);
        for (int k = 0; k < domain.nb_local(); k++) {
            local_kinetic += atoms.kin_energy(k);
        }
        MPI_Reduce(&local_kinetic, &global_kinetic, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        double gl_force = 0;
        double gr_force = 0;
        double l_force = 0;
        double r_force = 0;
        if (domain.rank() == 0) {
            for (int k = domain.nb_local(); k < atoms.positions.cols(); ++k) {
                if (atoms.positions(2, k) < 0) {
                    l_force += atoms.forces(2, k) * mass;
                }
            }
        }
        if (domain.rank() == domain.size() - 1) {
            for (int l = domain.nb_local(); l < atoms.positions.cols(); ++l) {
                if (atoms.positions(2, l) > z_scale) {
                    r_force += atoms.forces(2, l) * mass;
                }
            }
        }

        gl_force = MPI::allreduce(l_force, MPI_SUM, MPI_COMM_WORLD);
        gr_force = MPI::allreduce(r_force, MPI_SUM, MPI_COMM_WORLD);

        domain.disable(atoms);
        if (domain.rank() == 0 && i < 10000 && i % 100 == 0) {
            std::cout << global_kinetic * 2 /
                             (3 * 8.617333262e-5 * atoms.nb_atoms())
                      << std::endl;
        }
        if (i % save_every == 0 && domain.rank() == 0) {
            strain = (z_scale - domain_size) / domain_size;
            outdata << i << ", " << strain << ", " << gr_force << ", "
                    << gl_force << ", \n";
            std::cout << i << ", " << z_scale << ", " << strain << ", "
                      << gr_force << ", " << gl_force << std::endl;
            write_xyz("../xyz_output/milestone9_100_002_small/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }
        i++;
        domain.enable(atoms);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    outdata.close();
    MPI_Finalize();
    std::cout << "Milestone 9 ran successfully " << std::endl;
}

void lj_potential_distance(int steps, double mass, double sigma, double eps,
                           int save_gap, double rc) {
    double potential = 0;
    double kinetic_energy = 0;
    double total_energy = 0;
    int number = 0;

    std::ofstream outdata("../data/lj_distance.dat");

    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    int nbAtoms = 2;
    Positions_t positions(3, nbAtoms);
    Atoms atoms = {positions};
    atoms.masses = 1;
    atoms.positions(0, 0) = 0;
    atoms.positions(0, 1) = 1;
    // atoms.forces = 0;
    atoms.masses = 1;
    double timestep = .0001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));
    atoms.velocities = 0;
    for (int i = 0; i < steps; i++) {
        Verlet_one(timestep, atoms);
        lj_direct_summation_force(atoms, rc, eps, sigma);
        Verlet_two(timestep, atoms);
        potential = Potential(atoms, rc, eps, sigma);
        kinetic_energy = Kinetic(atoms, rc, eps, sigma);
        double total_energy = kinetic_energy + potential;

        std::cout << atoms.positions(0, 0) << " , " << atoms.positions(0, 1)
                  << std::endl;
        // std::cout  << potential  << std::endl;
        outdata << "[ " << atoms.positions(0, 1) + atoms.positions(0, 1)
                << " , " << potential << " , " << total_energy << " ],"
                << std::endl;
        // atoms.positions(0,1) += 0.5;
        // save xyz file
        if (i % save_gap == 0) {

            write_xyz("../xyz_output/lj_distance/" + filename +
                          to_string(number) + file_extension,
                      atoms);
            number = number + 1;
        }
    }

    outdata.close();
}
