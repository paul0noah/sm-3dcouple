//
//  ShapeMatchModel.hpp
//  dual-decompositions
//
//  Created by Paul Rötzer on 23.05.21.
//

#ifndef ShapeMatchModel_hpp
#define ShapeMatchModel_hpp

#include "helper/shape.hpp"
#include "shapeMatchModel/energyComputation/deformationEnergy.hpp"
#include "shapeMatchModel/constraintsComputation/constraints.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ctime>
#include <ILP_input.h>

#define DEBUG_SHAPE_MATCH_MODEL true


class ShapeMatchModelDijkstra {
private:
    bool verbose;
    Eigen::MatrixXd VX;
    Eigen::MatrixXi EX;
    Eigen::MatrixXd VY;
    Eigen::MatrixXi EY;
    Eigen::MatrixXd FeatDiffMatrix;

    Eigen::MatrixXi AI;
    Eigen::MatrixXi AJ;
    Eigen::MatrixXi AV;
    Eigen::MatrixXi RHS;

    Eigen::MatrixXi AIleq;
    Eigen::MatrixXi AJleq;
    Eigen::MatrixXi AVleq;
    Eigen::MatrixXi RHSleq;

    Eigen::MatrixXi productspace;
    Eigen::MatrixXi piEy;
    Eigen::MatrixXd energy;
    bool modelGenerated;
    int numCouplingConstraints;
    bool couplingConstraints;
    bool lineIntegral;
    bool otherSelfIntersections;
    
public:
    ShapeMatchModelDijkstra(Eigen::MatrixXd& iVX, Eigen::MatrixXi& iFX, Eigen::MatrixXd& iVY, Eigen::MatrixXi& iEY, Eigen::MatrixXd& iFeatDiffMatrix, bool iCouplingConstraints, bool iOtherSelfIntersections, bool iLineIntegral);
    ~ShapeMatchModelDijkstra();
    void generate();

    Eigen::MatrixXd getCostVector();
    std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> getAVectors();
    std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> getAleqVectors();
    Eigen::MatrixXi getRHS();
    Eigen::MatrixXi getRHSleq();
    Eigen::MatrixXi getProductSpace();
    int getNumCouplingConstraints();
    Eigen::MatrixXi getSortedMatching(const Eigen::MatrixXi& indicatorVector);
    LPMP::ILP_input getIlpObj();

};

#endif /* ShapeMatchModel_hpp */
