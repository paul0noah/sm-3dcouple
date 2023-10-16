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

    smm.def(py::init<Eigen::MatrixXd&, Eigen::MatrixXi&, Eigen::MatrixXd&, Eigen::MatrixXi&, Eigen::MatrixXd&, bool, bool, bool>());
    smm.def("generate", &ShapeMatchModelDijkstra::generate);
    smm.def("getCostVector", &ShapeMatchModelDijkstra::getCostVector);
    smm.def("getAVectors", &ShapeMatchModelDijkstra::getAVectors);
    smm.def("getRHS", &ShapeMatchModelDijkstra::getRHS);
    smm.def("getAleqVectors", &ShapeMatchModelDijkstra::getAleqVectors);
    smm.def("getRHSleq", &ShapeMatchModelDijkstra::getRHSleq);
    smm.def("getProductSpace", &ShapeMatchModelDijkstra::getProductSpace);
    smm.def("getNumCouplingConstraints", &ShapeMatchModelDijkstra::getNumCouplingConstraints);
    smm.def("getSortedMatching", &ShapeMatchModelDijkstra::getSortedMatching);
    smm.def("getIlpObj", &ShapeMatchModelDijkstra::getIlpObj);

    py::class_<LPMP::ILP_input>(handle, "ILP_instance")
            .def(py::init<>());
}
