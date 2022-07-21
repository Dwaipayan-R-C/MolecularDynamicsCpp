//
// Created by Dwaipayan on 14-05-2022.
//
#include "initialize.h"
#ifndef YAMD_FORCE_H
#define YAMD_FORCE_H

Array3Xd calculateAcceleration(Array3Xd pos){
    int i, j, k;
    double f, rSqd;
    double rij[3]; // position of i relative to j

    Array3Xd accln(3,nb_atoms);

    for (i = 0; i < nb_atoms; i++) {  // set all accelerations to zero
        for (k = 0; k < 3; k++) {
            accln(k,i) = 0;
        }
    }

    for (i = 0; i < nb_atoms-1; i++) {   // loop over all distinct pairs i,j
        for (j = i+1; j < nb_atoms; j++) {
            // initialize r^2 to zero
            rSqd = 0;
            for (k = 0; k < 3; k++) {
                //  component-by-componenent position of i relative to j

                rij[k] = (pos(k,i) - pos(k,j));
                //  sum of squares of the components
                rSqd += rij[k] * rij[k];

            }

            //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            //  Since the r has become squared to 14 and 8 in denominator, we multiply one rij with f
            //  This is done to reduce the computation time
//            f = 24 * eps* (2 * (pow(sigma,-7))*pow(rSqd, -7) - (pow(sigma,-7))*pow(rSqd, -4));
            f = (2 * pow(rSqd, -7) - pow(rSqd, -4));
//            std::cout << f <<std::endl;
            for (k = 0; k < 3; k++) {
                //  from F = ma, where m = 1 in natural units!
                accln(k,i) += rij[k] * f;
//                std::cout << "Acceleration, Pair Distance" <<std::endl;
//                std::cout << accln(k,i) <<std::endl;
//                std::cout << rij[k] <<std::endl;
//
//                std::cout << "\n" <<std::endl;
                accln(k,j) -= rij[k] * f;
            }
        }
    }
    return accln;
}
#endif //YAMD_FORCE_H
