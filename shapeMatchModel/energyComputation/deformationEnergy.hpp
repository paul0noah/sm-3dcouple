//
//  deformationEnergy.hpp
//  dual-decompositions
//
//  Created by Paul Rötzer on 31.03.21.
//

#ifndef deformationEnergy_hpp
#define deformationEnergy_hpp
#include <Eigen/Dense>
#include "helper/shape.hpp"

#define DEBUG_DEFORMATION_ENERGY false

class DeformationEnergy {
private:
    Eigen::MatrixXd& VX;
    Eigen::MatrixXd& VY;

    Eigen::MatrixXi& productspace;

    Eigen::MatrixXd& featDiffMatrix;
    bool lineIntegral;
    
    bool computed;
    Eigen::MatrixXd defEnergy;
    void computeEnergy();
    
public:
    DeformationEnergy(Eigen::MatrixXd& VX, Eigen::MatrixXd& VY, Eigen::MatrixXi& productspace, Eigen::MatrixXd& featDiffMatrix, bool iLineItegral);
    Eigen::MatrixXd getDeformationEnergy();
};

#endif /* deformationEnergy_hpp */

