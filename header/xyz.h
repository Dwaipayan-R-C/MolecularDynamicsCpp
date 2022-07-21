
#ifndef YAMD_XYZ_H
#define YAMD_XYZ_H
#include <iostream>
#include <fstream>

#include "atom_structure.h"


/*
 * Type Names_t, if not defined use
 * using Names_t = std::vector<std::string>;
 */
//using Names_t = std::vector<std::string>;
/*
 * Read positions from an XYZ file. XYZ files are structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z
 *         where Name is some name for the atom and X Y Z the position
 */
std::tuple<Positions_t> read_xyz(const std::string &filename);

/*
 * Read positions and velocities from an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
std::tuple<Positions_t, Velocities_t> read_xyz_with_velocities(const std::string &filename);

/*
 * Write positions and velocities to an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
void write_xyz(fstream &file, Atoms& atoms);

/*
 * Write positions and velocities to an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
void write_xyz(const std::string &filename, Atoms& atoms);

#endif //YAMD_XYZ_H