#include "constraints.hpp"
#include <vector>

Constraints::Constraints(Eigen::MatrixXi& EX, Eigen::MatrixXi& EY, Eigen::MatrixXi& productspace, Eigen::MatrixXi& piEy, bool coupling) :
    EX(EX), EY(EY), productspace(productspace), coupling(coupling), piEy(piEy) {

    nVX = EX.maxCoeff()+1;
    numCouplingConstraints = 0;
}


std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> Constraints::getConstraints() {
    const long nVY = EY.rows();
    Eigen::MatrixXi EY2Idx(EY.maxCoeff(), 1);
    EY2Idx = -EY2Idx.setOnes();
    for (int i = 0; i < nVY; i++) {
        EY2Idx(EY(i, 0), 0) = i;
    }
    const long nEX = EX.rows();
    long nnzConstraints = nVY * (nEX + nVX) + 2 * nVY * (2 * nEX + 2 * nVX);
    long numRowsConstraints = nVY + nVY * nVX;
    Eigen::ArrayXi VYcount;
    std::vector<std::vector<int>> couplingLayers;
    int numLayerCouples = 0;
    if (coupling) {
        const long numVY = EY.maxCoeff() + 1;
        for (int i = 0; i < numVY; i++) {
            couplingLayers.push_back(std::vector<int>());
        }
        // find double or triple vertices
        VYcount = Eigen::ArrayXi(numVY, 1);
        VYcount.setZero();
        for (int i = 0; i < EY.rows(); i++) {
            VYcount(EY(i, 0), 0) = VYcount(EY(i, 0), 0) + 1;
            couplingLayers[EY(i, 0)].push_back(i);
        }

        for (int i = 0; i < numVY; i++) {
            if (VYcount(i, 0) > 1) {
                numLayerCouples += VYcount(i, 0) - 1;
            }
        }

        if (numLayerCouples == 0) {
            if (DEBUG_CONSTRAINTS) std::cout << "[CONSTR] disabling coupling bc no intersections found" << std::endl;
            // no intersections found => no coupling constraints
            coupling = false;
        }
        else {
            numCouplingConstraints = numLayerCouples;
            if (DEBUG_CONSTRAINTS) std::cout << "[CONSTR] found "<< numLayerCouples << " coupling constraints" << std::endl;
            nnzConstraints += numLayerCouples * 2 * (nEX * 2 + nVX);
            numRowsConstraints += numLayerCouples * nVX;
        }
    }

    Eigen::MatrixXi I(nnzConstraints, 1); I.setZero();
    Eigen::MatrixXi J(nnzConstraints, 1); J.setZero();
    Eigen::MatrixXi V(nnzConstraints, 1); V.setZero();
    long idx = 0;

    Eigen::MatrixXi RHS(numRowsConstraints, 1);
    RHS.setZero();

    long offset = 0;
    // projection on inter layer edges for each layer has to be == 1 (for cycle)
    for (long i = 0; i < productspace.rows(); i++) {
        const bool interLayerEdge = productspace(i, 0) != productspace(i, 1);
        if (interLayerEdge) {
            const int rowIdx = piEy(i);//;EY2Idx(productspace(i, 0), 0);
            if (DEBUG_CONSTRAINTS && EY2Idx(productspace(i, 0), 0) < 0) {
                std::cout << "[CONSTR] EY2Idx(productspace(i, 0), 0) < 0 first loop" << std::endl;
            }

            I(idx, 0) = rowIdx;
            J(idx, 0) = (int) i;
            V(idx, 0) = 1;
            idx++;
            if (rowIdx > offset) {
                offset = rowIdx;
            }
        }

        if (DEBUG_CONSTRAINTS && idx >= nnzConstraints) {
            std::cout << "[CONSTR] idx beyond expected nonzeros in first loop" << std::endl;
            break;
        }
    }
    for (long i = 0; i < offset + 1; i++) {
        RHS(i, 0) = 1;
    }
    // the rest of RHS is zero i guess


    const long offset2 = offset + 1;
    // in and out of each node has to be 0 constraint
    Eigen::MatrixXi debugCounter(numRowsConstraints, 1);
    for (long i = 0; i < productspace.rows(); i++) {
        //const long inNodeIdx = EY2Idx(productspace(i, 1), 0) * nVX + productspace(i, 3); // we go into this node
        const long inNodeIdx = piEy(i, 0) * nVX + productspace(i, 2); // we go into this node
        if (DEBUG_CONSTRAINTS && EY2Idx(productspace(i, 1), 0) < 0) {
            std::cout << "[CONSTR] EY2Idx(productspace(i, 1), 0) < 0 second loop inNodeIdx" << std::endl;
        }
        //const long outNodeIdx = EY2Idx(productspace(i, 0), 0) * nVX + productspace(i, 2); // we leave this node
        const long outNodeIdx = piEy(i, 1) * nVX + productspace(i, 3); // we leave this node
        if (DEBUG_CONSTRAINTS && EY2Idx(productspace(i, 0), 0) < 0) {
            std::cout << "[CONSTR] EY2Idx(productspace(i, 0), 0) < 0 second loop outNodeIdx" << std::endl;
        }

        if (DEBUG_CONSTRAINTS) {
            if (inNodeIdx >= (nVY + 1) * nVX) {
                std::cout << "[CONSTR] inNodeIdx  >= nVY * nVX " << inNodeIdx << " productspace(i, 1), productspace(i, 3) "
                          << productspace(i, 1) << ", " << productspace(i, 3)<<  std::endl;
                continue;
            }
            if (outNodeIdx >= (nVY + 1) * nVX) {
                std::cout << "[CONSTR] outNodeIdx  >= nVY * nVX " << outNodeIdx << " productspace(i, 0), productspace(i, 2) "
                          << productspace(i, 0) << ", " << productspace(i, 2)<<  std::endl;
                continue;
            }
        }

        // in means +1
        I(idx, 0) = (int) (offset2 + inNodeIdx);
        J(idx, 0) = (int) i;
        V(idx, 0) = 1;
        idx++;

        // out means -1
        I(idx, 0) = (int) (offset2 + outNodeIdx);
        J(idx, 0) = (int) i;
        V(idx, 0) = -1;
        idx++;


        if (offset2 + inNodeIdx > offset) {
            offset = offset2 + inNodeIdx;;
        }

        if (offset2 + outNodeIdx > offset) {
            offset = offset2 + outNodeIdx;;
        }

        if (DEBUG_CONSTRAINTS && idx >= nnzConstraints) {
            std::cout << "[CONSTR] idx beyond expected nonzeros in second loop i="<<i << std::endl;
            break;
        }
    }

    if (DEBUG_CONSTRAINTS && offset + 1 - offset2 != nVY * nVX) {
        std::cout << "[CONSTR] offset - offset2 != nVY * nVX, offset - offset2 = "<< offset - offset2 << ", " << nVY * nVX << std::endl;
    }

    const long offset3 = offset + 1;
    // coupling of layers if the resemble same Y vertex and coupling is enabled
    if (coupling) {
        const long numEdgesOnLayer = (2 * nEX + nVX);
        int nthCoupling = 0;
        for (long i = 0; i < nVY; i++) {
            if (VYcount(i, 0) > 1) {
                int baseCoupleLayer = -1;
                for (auto& layer : couplingLayers.at(i)) {
                    if (baseCoupleLayer == -1) {
                        baseCoupleLayer = layer;
                        continue;
                    }


                    // any outgoing edge of vertex in X in layer i must be matched by any outgoing edge in baseCoupleLayer
                    const long edgeOffsetBase = baseCoupleLayer * numEdgesOnLayer;
                    const long edgeOffsetCouple = layer * numEdgesOnLayer;
                    for (long j = 0; j < numEdgesOnLayer; j++) {
                        const long rowIdx = nthCoupling * nVX + productspace(j, 2); // actually only on 0-th layer but should be fine

                        // base layer +1
                        I(idx, 0) = (int) (offset3 + rowIdx);
                        J(idx, 0) = (int) (edgeOffsetBase + j);
                        V(idx, 0) = 1;
                        idx++;

                        // coupled layer -1
                        I(idx, 0) = (int) (offset3 + rowIdx);
                        J(idx, 0) = (int) (edgeOffsetCouple + j);
                        V(idx, 0) = -1;
                        idx++;

                    }

                    nthCoupling++;
                }

                if (DEBUG_CONSTRAINTS && idx >= nnzConstraints) {
                    std::cout << "[CONSTR] idx beyond expected nonzeros in third loop" << std::endl;
                    break;
                }
            }
        }

    }

    I.conservativeResize(idx, 1);
    J.conservativeResize(idx, 1);
    V.conservativeResize(idx, 1);
    return std::make_tuple(I, J, V, RHS);
}



int Constraints::getNumCouplingConstr() {
    return numCouplingConstraints;
}