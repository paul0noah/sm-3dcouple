//
//  ShapeMatchModelPyBinds.cpp
//  dual-decompositions
//
//  Created by Paul Rötzer on 27.06.22.
//

#include "shapeMatchModel/shapeMatchModel.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>


namespace py = pybind11;
using namespace pybind11::literals;  // NOLINT

PYBIND11_MODULE(smm_dijkstra, handle) {
    handle.doc() = "ShapeMatchModelDijkstra python bindings";

    py::class_<ShapeMatchModelDijkstra, std::shared_ptr<ShapeMatchModelDijkstra>> smm(handle, "ShapeMatchModelDijkstra");

    smm.def(py::init<Eigen::MatrixXd&, Eigen::MatrixXi&, Eigen::MatrixXd&, Eigen::MatrixXi&, Eigen::MatrixXd&, bool, bool>());
    smm.def("generate", &ShapeMatchModelDijkstra::generate);
    smm.def("getCostVector", &ShapeMatchModelDijkstra::getCostVector);
    smm.def("getAVectors", &ShapeMatchModelDijkstra::getAVectors);
    smm.def("getRHS", &ShapeMatchModelDijkstra::getRHS);
    smm.def("getProductSpace", &ShapeMatchModelDijkstra::getProductSpace);
    smm.def("getNumCouplingConstraints", &ShapeMatchModelDijkstra::getNumCouplingConstraints);
}