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
    int currentstartidx = 0;
    for (int i = 0; i < nEy; i++) {
        if (EY(i, 0) == -1) {
            productspace.row(numadded) << -1, -1, -1, -1;
            piEY(numadded, 0) = i;
            piEY(numadded, 1) = i;
            numContours++;
            currentstartidx = i+1;
            numadded++;
            continue;
        }
        for (int j = 0; j < EX.rows(); j++) {

            // intralayer
            productspace.row(numadded) << EY(i, 0), EY(i, 0), EX(j, 0), EX(j, 1);
            piEY(numadded, 0) = i;
            piEY(numadded, 1) = i;
            numadded++;

            // interlayer
            productspace.row(numadded) << EY(i, 0), EY(i, 1), EX(j, 0), EX(j, 1);
            piEY(numadded, 0) = i;
            if (i+1 < nEy) {
                if (EY(i+1, 0) == -1) {
                    piEY(numadded, 1) = currentstartidx;
                }
                else {
                    piEY(numadded, 1) = i+1;
                }
            }
            else {
                piEY(numadded, 1) = currentstartidx;
            }
            numadded++;

        }

        for (int j = 0; j < numVerticesX; j++) {
            // interlayer
            productspace.row(numadded) << EY(i, 0), EY(i, 1), j, j;
            piEY(numadded, 0) = i;
            if (i+1 < nEy) {
                if (EY(i+1, 0) == -1) {
                    piEY(numadded, 1) = currentstartidx;
                }
                else {
                    piEY(numadded, 1) = i+1;
                }
            }
            else {
                piEY(numadded, 1) = currentstartidx;
            }
            numadded++;
        }

    }
    if (DEBUG_COMBINATIONS) std::cout << "[COMBOS] Detected " << numContours << " closed contours" << std::endl;
    productspace.conservativeResize(numadded, 4);
    combosComputed = true;
}


Combinations::Combinations(Eigen::MatrixXi& EX, Eigen::MatrixXi& EY) : EX(EX), EY(EY) {
    combosComputed = false;
    numContours = 1;
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

int Combinations::getNumContours() const {
    return numContours;
}
