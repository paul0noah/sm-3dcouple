//
//  ShapeMatchModel.cpp
//  dual-decompositions
//
//  Created by Paul Rötzer on 23.05.21.
//

#include "shapeMatchModel.hpp"
#include <fstream>
#include <iostream>
#include <cstdio>
#include "helper/utils.hpp"
#include <chrono>
#include <filesystem>
#include <algorithm>



void ShapeMatchModelDijkstra::generate() {

    if (verbose) std::cout << "[ShapeMM] Generating Shape Match Model..." << std::endl;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    if (verbose) std::cout << "[ShapeMM]   > Product Space" << std::endl;
    Combinations combos(EX, EY);
    productspace = combos.getProductSpace();
    Eigen::MatrixXi piEy = combos.getPiEy();
    
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    if (verbose) std::cout << "[ShapeMM]   Done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "  [ms])" << std::endl;
    if (verbose) std::cout << "[ShapeMM]   > Constraints" << std::endl;
    Constraints constr(EX, EY, productspace, piEy, couplingConstraints);
    const auto constrVectors = constr.getConstraints();
    AI  = std::get<0>(constrVectors);
    AJ  = std::get<1>(constrVectors);
    AV  = std::get<2>(constrVectors);
    RHS = std::get<3>(constrVectors);
    numCouplingConstraints = constr.getNumCouplingConstr();
    
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    if (verbose) std::cout << "[ShapeMM]   Done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "  [ms])" << std::endl;
    if (verbose) std::cout << "[ShapeMM]   > Energies" << std::endl;
    DeformationEnergy defEnergy(VX, VY, productspace, FeatDiffMatrix, lineIntegral);
    energy = defEnergy.getDeformationEnergy();
    
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    if (verbose) std::cout << "[ShapeMM]   Done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << "  [ms])" << std::endl;
    if (verbose) std::cout << "[ShapeMM] Done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t1).count() << "  [ms])" << std::endl;

    modelGenerated = true;
}



ShapeMatchModelDijkstra::ShapeMatchModelDijkstra(Eigen::MatrixXd& iVX, Eigen::MatrixXi& iEX, Eigen::MatrixXd& iVY, Eigen::MatrixXi& iEY, Eigen::MatrixXd& iFeatDiffMatrix, bool iCouplingConstraints, bool iLineIntegral) {
    VX = iVX;
    EX = iEX;
    VY = iVY;
    EY = iEY;
    FeatDiffMatrix = iFeatDiffMatrix;
    modelGenerated = false;
    verbose = true;
    numCouplingConstraints = 0;
    couplingConstraints = iCouplingConstraints;
    lineIntegral = iLineIntegral;
}


ShapeMatchModelDijkstra::~ShapeMatchModelDijkstra() {
    
}


Eigen::MatrixXd ShapeMatchModelDijkstra::getCostVector() {
    if (!modelGenerated) {
        generate();
    }
    return energy;
}

std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> ShapeMatchModelDijkstra::getAVectors() {
    if (!modelGenerated) {
        generate();
    }
    return std::make_tuple(AI, AJ, AV);
}

Eigen::MatrixXi ShapeMatchModelDijkstra::getRHS() {
    if (!modelGenerated) {
        generate();
    }
    return RHS;
}

Eigen::MatrixXi ShapeMatchModelDijkstra::getProductSpace() {
    if (!modelGenerated) {
        generate();
    }
    return productspace;
}

int ShapeMatchModelDijkstra::getNumCouplingConstraints() {
    return numCouplingConstraints;
}
