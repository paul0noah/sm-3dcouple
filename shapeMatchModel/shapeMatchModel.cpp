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
    Constraints constr(EX, EY, productspace, piEy, couplingConstraints, otherSelfIntersections);
    const auto constrVectors = constr.getConstraints();
    AI  = std::get<0>(constrVectors);
    AJ  = std::get<1>(constrVectors);
    AV  = std::get<2>(constrVectors);
    RHS = std::get<3>(constrVectors);
    AIleq  = std::get<4>(constrVectors);
    AJleq  = std::get<5>(constrVectors);
    AVleq  = std::get<6>(constrVectors);
    RHSleq = std::get<7>(constrVectors);
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



ShapeMatchModelDijkstra::ShapeMatchModelDijkstra(Eigen::MatrixXd& iVX, Eigen::MatrixXi& iEX, Eigen::MatrixXd& iVY, Eigen::MatrixXi& iEY, Eigen::MatrixXd& iFeatDiffMatrix, bool iCouplingConstraints, bool iOtherSelfIntersections, bool iLineIntegral) {
    VX = iVX;
    EX = iEX;
    VY = iVY;
    EY = iEY;
    FeatDiffMatrix = iFeatDiffMatrix;
    modelGenerated = false;
    verbose = true;
    numCouplingConstraints = 0;
    couplingConstraints = iCouplingConstraints;
    otherSelfIntersections = iOtherSelfIntersections;
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

std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> ShapeMatchModelDijkstra::getAleqVectors() {
    if (!modelGenerated) {
        generate();
    }
    return std::make_tuple(AIleq, AJleq, AVleq);
}

Eigen::MatrixXi ShapeMatchModelDijkstra::getRHSleq() {
    if (!modelGenerated) {
        generate();
    }
    return RHSleq;
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

Eigen::MatrixXi ShapeMatchModelDijkstra::getSortedMatching(const Eigen::MatrixXi& indicatorVector) {
    const long maxNumEdgesOnLevel = productspace.rows() / EY.rows();


    const long numCycleConstr = EY.rows();
    const long numInOutConstr = EY.rows() * VX.rows();
    std::vector<Eigen::Triplet<int>> in_out_entries;
    in_out_entries.reserve(AI.rows());
    for (long i = 0; i < AI.rows(); i++) {
        if (AI(i) >= numCycleConstr && AI(i) < numCycleConstr + numInOutConstr) {
            in_out_entries.push_back(Eigen::Triplet<int>(AI(i) - (int) numCycleConstr, AJ(i), AV(i)));
        }
    }

    Eigen::SparseMatrix<int, Eigen::RowMajor> in_out_rm(numInOutConstr, productspace.rows());
    in_out_rm.setFromTriplets(in_out_entries.begin(), in_out_entries.end());
    Eigen::SparseMatrix<int, Eigen::ColMajor> in_out_cm(numInOutConstr, productspace.rows());
    in_out_cm.setFromTriplets(in_out_entries.begin(), in_out_entries.end());


    long firstNonZeroIdx = -1;
    long numMatches = 0;
    for (long i = 0; i < indicatorVector.rows(); i++) {
        if (indicatorVector(i) == 1) {
            numMatches++;
            if (firstNonZeroIdx == -1)
                firstNonZeroIdx = i;
        }
    }
    Eigen::MatrixXi matchingSorted(numMatches, 4);
    matchingSorted = -matchingSorted.setOnes();
    Eigen::MatrixXi nodeUsed(numInOutConstr, 1); nodeUsed.setZero();
    //Eigen::MatrixXi sortedIndices(numMatches, 1); sortedIndices = -sortedIndices.setOnes();

    long idx = firstNonZeroIdx;
    for (long i = 0; i < numMatches; i++) {
        //sortedIndices(i, 1) = (int) idx;
        //std::cout << idx << ": " <<productspace.row(idx) << std::endl;
        matchingSorted.row(i) = productspace.row(idx);
        long row = -1;
        int currentVal = 0;
        bool newNodeFound = false;
        for (typename Eigen::SparseMatrix<int, Eigen::ColMajor>::InnerIterator it(in_out_cm, idx); it; ++it) {
            if (nodeUsed(it.row(), 0) == 0 && it.value() == -1) {
                //std::cout << "  " << it.value() << std::endl;
                row = it.row();
                nodeUsed(row, 0) = 1;
                currentVal = it.value();
                newNodeFound = true;
                break;
            }
        }
        if (!newNodeFound) {
            if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "[ShapeMM] Did not find new node, aborting" << std::endl;
            break;
        }

        for (typename Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(in_out_rm, row); it; ++it) {
            if (it.value() == -currentVal && indicatorVector(it.col(), 0) > 0) {
                idx = it.col();
                break;
            }
        }
    }

    return matchingSorted;
}
