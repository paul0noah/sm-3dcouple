//
//  ShapeMatchModel.cpp
//  dual-decompositions
//
//  Created by Paul Rötzer on 23.05.21.
//

#include "shapeMatchModel.hpp"

LPMP::ILP_input ShapeMatchModelDijkstra::getIlpObj() {
    LPMP::ILP_input ilp = LPMP::ILP_input();
    Eigen::MatrixXd objective = getCostVector();

    // build constraint matrices
    const std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> Avecs = getAVectors();
    const Eigen::MatrixXi AI = std::get<0>(Avecs);
    const Eigen::MatrixXi AJ = std::get<0>(Avecs);
    const Eigen::MatrixXi AV = std::get<0>(Avecs);

    const std::tuple<Eigen::MatrixXi, Eigen::MatrixXi, Eigen::MatrixXi> Avecsleq = getAleqVectors();
    const Eigen::MatrixXi AleqI = std::get<0>(Avecsleq);
    const Eigen::MatrixXi AleqJ = std::get<0>(Avecsleq);
    const Eigen::MatrixXi AleqV = std::get<0>(Avecsleq);

    const Eigen::MatrixXi rhs = getRHS();
    const Eigen::MatrixXi rhsleq = getRHSleq();

    std::vector<Eigen::Triplet<int>> entriesA;
    entriesA.reserve(AI.rows());
    for (long i = 0; i < AI.rows(); i++) {
        entriesA.push_back(Eigen::Triplet<int>(AI(i, 0), AJ(i, 0), AV(i, 0)));
    }
    Eigen::SparseMatrix<int, Eigen::RowMajor> A(AI.maxCoeff()+1, objective.rows()); A.setFromTriplets(entriesA.begin(), entriesA.end());

    std::vector<Eigen::Triplet<int>> entriesAleq;
    entriesAleq.reserve(AleqI.rows());
    for (long i = 0; i < AleqI.rows(); i++) {
        entriesAleq.push_back(Eigen::Triplet<int>(AleqI(i, 0), AleqJ(i, 0), AleqV(i, 0)));
    }
    Eigen::SparseMatrix<int, Eigen::RowMajor> Aleq(AleqI.maxCoeff()+1, objective.rows()); Aleq.setFromTriplets(entriesAleq.begin(), entriesAleq.end());


    // Add variables to ilp
    for (int i = 0; i < objective.rows(); i++) {
        std::string varName = "x" + std::to_string(i);
        ilp.add_new_variable(varName);
        ilp.add_to_objective((double) objective(i), i);
    }


    // Add constraints to ilp
    /*
        beginNewInequality   => creates new constraint
        inequalityIdentifier => name of the constraint e.g. R101
        inequalityType       => in our case always "="
    */
    long numadded = 0;
    for (int k = 0; k < A.outerSize(); ++k) {
        ilp.begin_new_inequality();
        const std::string identifier = "R" + std::to_string(numadded) + " ";
        ilp.set_inequality_identifier(identifier);
        ilp.set_inequality_type(LPMP::ILP_input::inequality_type::equal);
        for (typename Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it) {
                ilp.add_to_constraint(it.value(), it.index());
        }
        ilp.set_right_hand_side(rhs(k, 0));
        numadded++;
    }
    for (int k = 0; k < Aleq.outerSize(); ++k) {
        ilp.begin_new_inequality();
        const std::string identifier = "R" + std::to_string(numadded) + " ";
        ilp.set_inequality_identifier(identifier);
        ilp.set_inequality_type(LPMP::ILP_input::inequality_type::equal);
        for (typename Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(Aleq, k); it; ++it) {
                ilp.add_to_constraint(it.value(), it.index());
        }
        ilp.set_right_hand_side(rhsleq(k, 0));
        numadded++;
    }

    return ilp;
}
