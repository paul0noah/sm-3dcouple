//
//  combinations.hpp
//  helper
//
//  Created by Paul Rötzer on 28.04.21.
//

#ifndef combinations_hpp
#define combinations_hpp

#include <Eigen/Dense>
#include "helper/shape.hpp"
#define DEBUG_COMBINATIONS true

class Combinations {
private:
    bool combosComputed;
    Eigen::MatrixXi& EX;
    Eigen::MatrixXi& EY;
    Eigen::MatrixXi productspace;
    Eigen::MatrixXi piEY;
    int numContours;

public:
    void init();
    Combinations(Eigen::MatrixXi& EX, Eigen::MatrixXi& EY);
    
    void computeCombinations();
    Eigen::MatrixXi getProductSpace();
    Eigen::MatrixXi getPiEy();
    int getNumContours() const;
};

#endif /* combinations_hpp */
