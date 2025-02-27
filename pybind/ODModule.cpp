#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "OverdampedLangevin.h"

namespace py = pybind11;

PYBIND11_MODULE(BrownianMotion, m) {
    m.doc() = "Overdamped Langevin module"; // optional module docstring

    py::class_<vec<1>>(m, "vec1")
        .def(py::init<std::array<double, 1>>())
        .def_readwrite("x", &vec<1>::x);

    py::class_<OverdampedLangevin<1>>(m, "OverdampedLangevin")
        .def(py::init<double, OverdampedLangevin<1>::force_callable, vec<1>>())
        .def(py::init<double, vec<1>>())
        .def("eulerMaruyamaStep", &OverdampedLangevin<1>::eulerMaruyamaStep)
        .def("getTrajectory", &OverdampedLangevin<1>::getTrajectoryPy);

    py::class_<BatchedOverdampedLangevin<1>>(m, "BatchedOverdampedLangevin")
        .def(py::init<double, BatchedOverdampedLangevin<1>::force_callable, std::vector<vec<1>>>())
        .def(py::init<double, std::vector<vec<1>>>())
        .def("getBatchedTrajectory", &BatchedOverdampedLangevin<1>::getBatchedTrajectoryPy);
}