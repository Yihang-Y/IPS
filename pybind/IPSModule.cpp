#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "IPS.h"
#include "Integrator.h"
#include "Particles.h"
#include "Potential.h"


namespace py = pybind11;

template struct Particles<DataLayout::SoA, 2>;

ConfigValue cast_config_value(py::handle value) {
    if (py::isinstance<py::int_>(value)) {
        return value.cast<int>();
    } else if (py::isinstance<py::float_>(value)) {
        return value.cast<double>();
    } else if (py::isinstance<py::bool_>(value)) {
        return value.cast<bool>();
    } else if (py::isinstance<py::str>(value)) {
        return value.cast<std::string>();
    } else {
        throw std::runtime_error("Unsupported config value type");
    }
}

PYBIND11_MODULE(IPSModule, m) {
    m.doc() = "IPS module"; // optional module docstring

    py::class_<IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>>(m, "IPS_Simulator")
        .def(py::init<Particles<DataLayout::SoA, 2>&>())
        .def("init", [](IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>& self,
                        py::dict pair_force_config,
                        py::dict confinement_force_config) {
            Config pair_force_config_cpp;
            Config confinement_force_config_cpp;
            for (auto item : pair_force_config) {
                pair_force_config_cpp[item.first.cast<std::string>()] = cast_config_value(item.second);
            }
            for (auto item : confinement_force_config) {
                confinement_force_config_cpp[item.first.cast<std::string>()] = cast_config_value(item.second);
            }
            self.pair_force = make_pair_force(pair_force_config_cpp);
            self.confinement_force = make_confinement_force(confinement_force_config_cpp);
        })
        .def("integrate", &IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>::integrate)
        .def("integrate_n_steps", &IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>::integrate_n_steps);

    // py::class_<IPS_Simulator<Particles<DataLayout::SoA, 2>, BAOAB>>(m, "IPS_Simulator")
    //     .def(py::init<Particles<DataLayout::SoA, 2>&>())
    //     .def("init", [](IPS_Simulator<Particles<DataLayout::SoA, 2>, BAOAB>& self,
    //                     py::dict pair_force_config,
    //                     py::dict confinement_force_config) {
    //         Config pair_force_config_cpp;
    //         Config confinement_force_config_cpp;
    //         for (auto item : pair_force_config) {
    //             pair_force_config_cpp[item.first.cast<std::string>()] = cast_config_value(item.second);
    //         }
    //         for (auto item : confinement_force_config) {
    //             confinement_force_config_cpp[item.first.cast<std::string>()] = cast_config_value(item.second);
    //         }
    //         self.pair_force = make_pair_force(pair_force_config_cpp);
    //         self.confinement_force = make_confinement_force(confinement_force_config_cpp);
    //     })
    //     .def("integrate", &IPS_Simulator<Particles<DataLayout::SoA, 2>, BAOAB>::integrate)
    //     .def("integrate_n_steps", &IPS_Simulator<Particles<DataLayout::SoA, 2>, BAOAB>::integrate_n_steps);

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
        .def("get_forces", &Particles<DataLayout::SoA, 2>::get_forces)
        .def("get_temperature", &Particles<DataLayout::SoA, 2>::get_temperature);


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