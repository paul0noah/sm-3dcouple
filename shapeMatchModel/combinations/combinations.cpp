//
//  combinations.cpp
//  helper
//
//  Created by Paul Rötzer on 28.04.21.
//

#include "combinations.hpp"
#include "helper/utils.hpp"


/*function computeCombinations(...)
 Consider the edge matrix of shape X (both orientations bc triangle mesh) and the one of shape Y (single edge orientation bc 3D contour)
*/
void Combinations::computeCombinations() {
    const int numVerticesX = EX.maxCoeff() + 1;
    productspace = Eigen::MatrixXi(EY.rows() * ( 2 * EX.rows() + numVerticesX), 4);
    productspace.setZero();
    piEY = Eigen::MatrixXi(EY.rows() * ( 2 * EX.rows() + numVerticesX), 2);
    piEY.setZero();

    const long nEy = EY.rows();
    long numadded = 0;
    for (int i = 0; i < nEy; i++) {
        for (int j = 0; j < EX.rows(); j++) {

            // intralayer
            productspace.row(numadded) << EY(i, 0), EY(i, 0), EX(j, 0), EX(j, 1);
            piEY(numadded, 0) = i;
            piEY(numadded, 1) = i;
            numadded++;

            // interlayer
            productspace.row(numadded) << EY(i, 0), EY(i, 1), EX(j, 0), EX(j, 1);
            piEY(numadded, 0) = i;
            piEY(numadded, 1) = (i+1) % nEy;
            numadded++;

        }

        for (int j = 0; j < numVerticesX; j++) {
            // interlayer
            productspace.row(numadded) << EY(i, 0), EY(i, 1), j, j;
            piEY(numadded, 0) = i;
            piEY(numadded, 1) = (i+1) % nEy;
            numadded++;
        }

    }
    combosComputed = true;
}


Combinations::Combinations(Eigen::MatrixXi& EX, Eigen::MatrixXi& EY) : EX(EX), EY(EY) {
    combosComputed = false;
}



Eigen::MatrixXi Combinations::getProductSpace() {
    if (!combosComputed) {
        computeCombinations();
    }
    return productspace;
}

Eigen::MatrixXi Combinations::getPiEy() {
    if (!combosComputed) {
        computeCombinations();
    }
    return piEY;
}
