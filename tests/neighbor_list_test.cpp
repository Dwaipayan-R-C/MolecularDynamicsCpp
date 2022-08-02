#include "../header/atom_structure.h"
#include "../header/neighbors.h"
#include "../header/xyz.h"

#include <gtest/gtest.h>

TEST(NeighborsTest, NeighborListTest) {    
    Positions_t positions(3, 4);
    positions << 0, 1, 0, 0,
            0, 0, 1, -1,
            0, 0, 0, 0;

    Atoms atoms(positions);
    NeighborList neighbor_list(1.5);
    auto &[seed, neighbors]{neighbor_list.update(atoms)};

    // All atoms except 3 and 4 are neighbors of each other
    EXPECT_EQ(neighbor_list.nb_neighbors(), 10);

    EXPECT_EQ(neighbor_list.nb_neighbors(0), 3);
    EXPECT_EQ(neighbor_list.nb_neighbors(1), 3);
    EXPECT_EQ(neighbor_list.nb_neighbors(2), 2);
    EXPECT_EQ(neighbor_list.nb_neighbors(3), 2);

    EXPECT_TRUE((neighbors(Eigen::seq(seed(0), seed(1) - 1)) == Eigen::Array3i{3, 1, 2}).all());
    EXPECT_TRUE((neighbors(Eigen::seq(seed(1), seed(2) - 1)) == Eigen::Array3i{3, 0, 2}).all());
    EXPECT_TRUE((neighbors(Eigen::seq(seed(2), seed(3) - 1)) == Eigen::Array2i{0, 1}).all());
    EXPECT_TRUE((neighbors(Eigen::seq(seed(3), seed(4) - 1)) == Eigen::Array2i{0, 1}).all());
}