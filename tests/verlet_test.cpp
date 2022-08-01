#include <gtest/gtest.h>
#include "../header/verlet.h"




TEST(VerletTest, SingleAtomVerletConstantForce)
{
    //single Atom test
    int nbAtoms = 1;
    Positions_t  positions(3,nbAtoms);
    Velocities_t  velocities(3, nbAtoms);
    Forces_t forces(3,nbAtoms);
    
    Atoms atoms{
        positions, velocities};
    atoms.positions.setZero();
    atoms.velocities.setZero();
    atoms.forces.setZero();
    atoms.forces.row(0) = 1.0;
    atoms.forces.row(0) = 1.0;
    double timestep = 1.0;
    
    //run fkt
    for (int i = 0; i < 10; i++)
    {
        Verlet_one(timestep, atoms);
        Verlet_two(timestep, atoms); 
    }      
    //check
    ASSERT_NEAR(atoms.positions(0),50,1e-6);
    ASSERT_NEAR(atoms.velocities(0),10,1e-6);
}


TEST(VerletTest, MultipleAtomVerletConstantForce)
{
    //single Atom test
    int nbAtoms = 10;
    Positions_t  positions(3,nbAtoms);
    positions.setZero();
    Atoms atoms(positions);
    atoms.forces.row(0) = 1.0;
    atoms.forces.row(1) = 1.0;
    atoms.forces.row(2) = 1.0;
    double timestep = 1.0;
    //run fkt
    for (int i = 0; i < 10; i++)
    {
        Verlet_one(timestep, atoms);
        Verlet_two(timestep, atoms); 
    } 
    //check
    for(int i = 0; i < nbAtoms; ++i)
    {
        //check positions
        ASSERT_NEAR(atoms.positions(0,i),50,1e-6);
        ASSERT_NEAR(atoms.positions(1,i),50,1e-6);
        ASSERT_NEAR(atoms.positions(2,i),50,1e-6);
        //check velocities
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
    }
}