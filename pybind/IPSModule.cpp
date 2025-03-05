#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "IPS.h"
#include "Integrator.h"
#include "Particles.h"
#include "Potential.h"


namespace py = pybind11;

template struct Particles<DataLayout::SoA, 2>;

PYBIND11_MODULE(IPSModule, m) {
    m.doc() = "IPS module"; // optional module docstring

    py::class_<IPS_Simulator<DataLayout::SoA, 2, LeapFrog>>(m, "IPS_Simulator")
        // .def(py::init<IPS_Simulator<DataLayout::SoA, 2, LeapFrog>::pair_force_callable, 
        //               IPS_Simulator<DataLayout::SoA, 2, LeapFrog>::confinement_force_callable,
        //               Particles<DataLayout::SoA, 2>&>(),
        //               LeapFrog())
        .def(py::init<Particles<DataLayout::SoA, 2>&>())
        // .def(py::init<SpringConfig, RadialConfinementConfig, Particles<DataLayout::SoA, 2>&>())
        // .def(py::init<Particles<DataLayout::SoA, 2>&>())
        .def("integrate", &IPS_Simulator<DataLayout::SoA, 2, LeapFrog>::integrate)
        .def("integrate_n_steps", &IPS_Simulator<DataLayout::SoA, 2, LeapFrog>::integrate_n_steps);

    py::class_<LeapFrog>(m, "LeapFrog")
        .def(py::init<>());
        // .def("integrate_2d",
        //      [](LeapFrog &self,
        //         Particles<DataLayout::SoA, 2>& particles,
        //         std::function<double(double)> pair_force,
        //         std::function<vec<2>(const vec<2>&)> confinement_force,
        //         double step_size) {
        //             self.template integrate_on_soa<2>(particles, pair_force, confinement_force, step_size);
        //      },
        //      py::arg("particles"),
        //      py::arg("pair_force"),
        //      py::arg("confinement_force"),
        //      py::arg("step_size")
        // )

    py::class_<Particles<DataLayout::SoA, 2>>(m, "Particles")
        .def(py::init([](size_t n_particles) {
            Particles<DataLayout::SoA, 2> p;
            for (size_t d = 0; d < 2; ++d) {
                p.positions[d].resize(n_particles);
                p.velocities[d].resize(n_particles);
                p.forces[d].resize(n_particles);
            }
            return p;
        }), py::arg("n_particles"))
        .def_readwrite("positions", &Particles<DataLayout::SoA, 2>::positions)
        .def_readwrite("velocities", &Particles<DataLayout::SoA, 2>::velocities)
        .def_readwrite("forces", &Particles<DataLayout::SoA, 2>::forces)
        .def("get_positions", &Particles<DataLayout::SoA, 2>::get_positions)
        .def("get_velocities", &Particles<DataLayout::SoA, 2>::get_velocities)
        .def("get_forces", &Particles<DataLayout::SoA, 2>::get_forces);

    py::class_<SpringConfig>(m, "SpringConfig")
        .def(py::init<double, double>(), py::arg("k"), py::arg("r0"));
    
    py::class_<RadialConfinementConfig>(m, "RadialConfinementConfig")
        .def(py::init<double>(), py::arg("rad"));

    // NOTE: useless, if it is bind to python, when C++ code call operator(), it will call the python function
    // py::class_<Spring>(m, "Spring")
    //     .def(py::init<>())
    //     .def(py::init<double, double>())
    //     .def("__call__", &Spring::operator());
    
    // py::class_<RadialConfinement<2>>(m, "RadialConfinement")
    //     .def(py::init<>())
    //     .def(py::init<double>())
    //     .def("__call__", &RadialConfinement<2>::operator());
}