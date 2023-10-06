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

Eigen::MatrixXi ShapeMatchModelDijkstra::getSortedMatching(const Eigen::MatrixXi& indicatorVector) {
    const long maxNumEdgesOnLevel = productspace.rows() / EY.rows();
    Eigen::MatrixXi matchingSorted(indicatorVector.rows(), 4);
    matchingSorted = -matchingSorted.setOnes();


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
        if (indicatorVector(i) > 0) {
            if (firstNonZeroIdx == -1)
                firstNonZeroIdx = i;
            numMatches++;
        }
    }
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

/*
    const long maxNumEdgesOnLevel = productspace.rows() / EY.rows();
    Eigen::MatrixXi matchingSorted(indicatorVector.rows(), 4);
    matchingSorted = -matchingSorted.setOnes();

    long numMatches = 0;
    for (long i = 0; i < indicatorVector.rows(); i++) {
        if (indicatorVector(i, 0) > 0) {
            matchingSorted.row(numMatches) = productspace.row(i);
            numMatches++;
            //std::cout << productspace.row(i) << std::endl;
        }
    }
    matchingSorted.conservativeResize(numMatches, 4);

    // by construction of product space we dont have to sort Y indices
    // we only have to resort rows according to X indices
    Eigen::MatrixXi tempRow(1, 4);
    for (long i = 1; i < numMatches; i++) {
        if (matchingSorted(i-1, 3) != matchingSorted(i, 2)) {
            if (DEBUG_SHAPE_MATCH_MODEL) std::cout << matchingSorted.row(i) << std::endl;
            if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "next: " << matchingSorted.row(i+1) << std::endl;
            tempRow.row(0) = matchingSorted.row(i);


            // find the actual next element
            long jLoopCounter = 0;
            for (long j = i+1; j < numMatches; j++) {
                if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "     checking: " << matchingSorted.row(j) << std::endl;
                if (matchingSorted(i-1, 3) == matchingSorted(j, 2) && matchingSorted(i-1, 1) == matchingSorted(j, 0)) {
                    matchingSorted.row(i) = matchingSorted.row(j);
                    matchingSorted.row(j) = tempRow.row(0);
                    break;
                }

                // make sure we dont spin forever
                if (jLoopCounter >= maxNumEdgesOnLevel) {
                    if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "[ShapeMM] Error searching for correct next edge, spinning to long :( " << std::endl;
                    break;
                }
                jLoopCounter++;
            }
        }
    }


    // resolve conflicts which might occur in the end of the matching
    std::vector<long> conflictIndices; conflictIndices.reserve(1000);
    for (long i = 1; i < numMatches; i++) {
        if (matchingSorted(i-1, 3) != matchingSorted(i, 2)) {
            conflictIndices.push_back(i);
        }
    }
    if (conflictIndices.size() == 0) {
        if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "[ShapeMM] No conflicts detected" << std::endl;
    }
    else {
        return matchingSorted;
        Eigen::MatrixXi resolvedSortedMatching = matchingSorted;
        long closingIdx = -1;
        for (auto& conflictIndex : conflictIndices) {
            if (conflictIndex > closingIdx && closingIdx != -1) {
                continue;
            }
            else {
                closingIdx = -1;
            }
            const int conflictstart2d = matchingSorted(conflictIndex, 0);
            const int conflictstart3d = matchingSorted(conflictIndex, 2);
            // search for closing loop
            for (long i = conflictIndex + 1; i < numMatches; i++) {
                if (matchingSorted(i, 1) == conflictstart2d && matchingSorted(i, 3) == conflictstart3d) {
                    closingIdx = i;
                }
            }
            if (closingIdx == -1) {
                if (DEBUG_SHAPE_MATCH_MODEL) std::cout << "[ShapeMM] Couldnt resolve conflict. aborting" << std::endl;
                break;
            }
            else {
                for (long i = 0; i < numMatches; i++) {
                    if (resolvedSortedMatching(i, 1) == conflictstart2d && resolvedSortedMatching(i, 3) == conflictstart3d) {
                        Eigen::MatrixXi tempRows = matchingSorted.block(i+1, 0, numMatches - conflictIndex, 4);
                        resolvedSortedMatching.block(i, 0, conflictIndex - closingIdx, 4) = matchingSorted.block(conflictIndex, 0, conflictIndex - closingIdx, 4);
                        resolvedSortedMatching.block(i + conflictIndex - closingIdx, 0, numMatches - conflictIndex, 4) = tempRows;
                        break;
                    }
                }
                matchingSorted = resolvedSortedMatching;
            }
        }

    }


    return matchingSorted;*/
}
