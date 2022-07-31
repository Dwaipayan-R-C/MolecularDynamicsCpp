#include "../header/milestones.h"
#include <list>
using namespace std::chrono;
using namespace Eigen;
using std::ofstream;

string names;
double potential = 0;
double kinetic_energy = 0;
double total_energy = 0;
int count_relax = 0;
int relax_value = 0;

void milestone4()
{

    int steps = 10000;
    double mass = 1;
    std::ofstream outdata("../data/milestone4.dat");

    if (!outdata)
    { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{read_xyz_with_velocities("../../xyz/lj54.xyz")};
    Atoms atoms{
        positions, velocities};
    atoms.masses = 1;
    int number = 0;
    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));
    int save_gap = 500;
    outdata << "[" << std::endl;
    for (int i = 0; i < steps; i++)
    {
        Verlet_one(timestep, atoms);
        lj_direct_summation_force(atoms);
        Verlet_two(timestep, atoms);
        potential = Potential(atoms);
        kinetic_energy = Kinetic(atoms);
        double total_energy = kinetic_energy + potential;
        outdata << "[ " << kinetic_energy
                << " , "
                << potential
                << " , "
                << total_energy
                << " ],"
                << std::endl;

        // save xyz file
        if (i % save_gap == 0)
        {
            std::cout << i << " steps completed" << std::endl;
            write_xyz("../xyz_output/milestone4/" + filename + to_string(number) + file_extension, atoms);
            number = number + 1;
        }
    }
    outdata << "]" << std::endl;
    outdata.close();
    std::cout << "Ran Velocity verlet successfully and output has been written to data/milestone4.dat in the format - Kinetic Energy, Potential Energy, Total Energy" << std::endl;
}

void milestone5()
{

    double mass = 1;
    int steps = 5000;
    double target_temp = 5000;
    double boltzmann_kb = 1;
    std::ofstream outdata("../data/milestone5.dat");

    if (!outdata)
    { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{read_xyz_with_velocities("../../xyz/lj54.xyz")};
    Atoms atoms{
        positions, velocities};
    atoms.masses = mass;
    int number = 0;
    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));
    int save_gap = 200;

    for (int i = 0; i < steps; i++)
    {
        Verlet_one(timestep, atoms);
        lj_direct_summation_force(atoms);
        Verlet_two(timestep, atoms);
        atoms.velocities = berendsen_thermostat(atoms, target_temp, timestep, 10 * timestep, boltzmann_kb);
        potential = Potential(atoms);
        kinetic_energy = Kinetic(atoms);

        double current_temp = (2 * kinetic_energy) / 3;
        double total_energy = kinetic_energy + potential;
        outdata
            << current_temp
            << " ,"
            << target_temp
            << std::endl;

        // save xyz file
        if (i % save_gap == 0)
        {
            std::cout << i << " steps completed" << std::endl;
            write_xyz("../xyz_output/milestone5/" + filename + to_string(number) + file_extension, atoms);
            number = number + 1;
        }
    }

    outdata.close();
    std::cout << "Ran Berdendsen Thermostat successfully and output has been written to data/milestone5.dat in the format - Total Energy , Current temperature, Target temperature" << std::endl;
}

void milestone6()
{
    double mass = 1;
    int steps = 1000;
    // double target_temp = 5000;
    double boltzmann_kb = 1;
    std::ofstream outdata("../data/milestone6.dat");

    std::list<int> rc_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    if (!outdata)
    { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions, velocities]{read_xyz_with_velocities("../../xyz/cluster_923.xyz")};
    Atoms atoms{
        positions, velocities};
    atoms.masses = mass;
    int number = 0;
    double timestep = .001 * pow((mass * pow(sigma, 2) / eps), (1 / 2));    

    for (int x : rc_list)
    {
        
        auto start = high_resolution_clock::now();

        NeighborList neighbor_list(x);
        for (int i = 0; i < steps; i++)
        {
            Verlet_one(timestep, atoms);
            neighbor_list.update(atoms);
            lj_direct_summation_force(atoms);
            Verlet_two(timestep, atoms);
            potential = Potential(atoms);
            kinetic_energy = Kinetic(atoms);
            double total_energy = kinetic_energy + potential;            
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        // in milisecond
        outdata << "[ "
            << x
            << " ,"
            << duration.count()/1000
            << " ],"
            << std::endl;
        
        std::cout
            << "Running"
            << std::endl;
    }

    outdata.close();
    std::cout << "Ran cutoff radius successfully and output has been written to data/milestone6.dat in the format - radius, time" << std::endl;
}

// int main(int argc, char *argv[])
// {
//     //     MPI_Init(&argc, &argv);

//     //    Domain domain(MPI_COMM_WORLD, {30, 30, 30}, {1, 1, 2}, {0, 0, 1});
//     //     Domain domain(MPI_COMM_WORLD, {40, 40, 40}, {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});
//     // reads the initial atomic position and velocity from xyz file
//     // assigns to the corresponding

//     clock_t start, end;
//     start = clock();
//     auto [positions, velocities]{read_xyz_with_velocities("../../xyz/li_54.xyz")};
//     Atoms atoms{
//         positions, velocities};
//     atoms.masses = mass;
//     // gets the initial force from Lennard Jones potential
//     // atoms.forces = lj_direct_summation_force(atoms);
//     int number = 0;
//     double temperature = 0;
//     double totalEnergy = 0;
//     int k = 0;
//     double alpha = 0;
//     double distnace;
//     NeighborList neighbor_list(rc);

//     for (int i = 0; i < nb_steps; i++)
//     {
//         int rank;
//         // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         // domain.enable(atoms);

//         // Verlet prediction 1
//         Verlet_one(timescale, atoms);

//         // Update force for Lennard jones
//         // atoms.forces = lj_direct_summation_force(atoms);
//         // GUPTA
//         //        std::cout<<atoms.nb_atoms()<<std::endl;
//         // domain.exchange_atoms(atoms);

//         //        std::cout<<atoms.nb_atoms()<<std::endl;
//         // domain.update_ghosts(atoms, 2 * rc);

//         //        std::cout<<atoms.nb_atoms()<<std::endl;
//         neighbor_list.update(atoms);

//         //        std::cout<<"neighbor_list.update worked "<<rank<<std::endl;
//         //        std::cout<<atoms.nb_atoms()<<std::endl;
//         potential = gupta(atoms, neighbor_list, rc);

//         // Verlet 2 propagation with updated force
//         Verlet_two(timescale, atoms);
//         //        std::cout<<"Verlet_two worked "<<rank<<std::endl;
//         // atoms.velocities = berendsen_thermostat(atoms,1500000, timescale, 10*timescale, kb);

//         // Calculate Kinetic and Lennard Jones Potential Energy
//         // potential = Potential(atoms);
//         kinetic_energy = Kinetic(atoms);
//         //        std::cout<<potential<<rank<<std::endl;
//         total_energy = potential + kinetic_energy;

//         // domain.disable(atoms);

//         // save xyz file
//         if (i % save_every == 0)
//         {
//             write_xyz(filename + to_string(number) + file_extension, atoms);
//             number = number + 1;
//         }
//         //        std::cout << total_energy << std::endl;
//         if (count_relax >= measurement_gap)
//         {
//             temperature += kinetic_energy * 2 / (3 * kb * nb_atoms);
//             totalEnergy += total_energy;
//             relax_value += 1;
//             //            count_relax = 0;
//         }

//         if (i % tau == 0 && i != 0)
//         {
//             std::cout << "[" << totalEnergy / (relax_value) << " , " << temperature / (relax_value) << "]," << std::endl;
//             relax_value = 0;
//             count_relax = 0;
//             temperature = 0;
//             totalEnergy = 0;
//             // Determine alpha
//             // del Q = Ekinetic' - Ekinetic
//             // del Q = Ekinetic (sq of alpha - 1)
//             alpha = sqrt((delQ / kinetic_energy) + 1);
//             atoms.velocities = atoms.velocities * alpha;
//             k = k + 1;
//         }
//         count_relax += 1;
//     }

//     end = clock();
//     printf("Time taken: %f secs\n", ((float)end - start) / CLOCKS_PER_SEC);
//     //     MPI_Finalize();
//     return 0;
// }
