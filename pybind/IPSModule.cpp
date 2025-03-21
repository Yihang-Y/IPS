#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "IPS.h"
#include "Integrator.h"
#include "Particles.h"
#include "Potential.h"


namespace py = pybind11;

template struct Particles<DataLayout::SoA, 2>;
template struct LangevinSystem<2>;

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
    m.doc() = "IPS module: Provides interfaces for the IPS simulator.\n\n"
              "This module supports several simulator classes:\n"
              "  - IPS_Simulator (using the LeapFrog integrator and Particles system)\n"
              "  - IPS_Simulator_Langevin (using the ABOBA integrator and LangevinSystem)\n"
              "  - IPS_Simulator_NoseHoover (using the NoseHoover integrator and NoseHooverSystem)\n"
              "  - IPS_Simulator_NoseHooverLangevin (using the NoseHooverLangevin integrator and NoseHooverLangevinSystem)\n\n"
              "Corresponding particle systems include: Particles, LangevinSystem, NoseHooverSystem, and NoseHooverLangevinSystem.\n\n"
              "Example usage:\n"
              "    >>> from IPSModule import IPS_Simulator\n"
              "    >>> simulator = IPS_Simulator(particles)\n"
              "    >>> simulator.init({'force_type': 'pair'}, {'confinement': 'type1'})\n"
              "    >>> simulator.integrate_n_steps(100)\n";

    // IPS_Simulator using Particles and LeapFrog
    py::class_<IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>>(m, "IPS_Simulator",
        "IPS_Simulator class for performing simulations using the LeapFrog integrator on a Particles system.\n\n"
        "This simulator supports setting up pair forces and confinement forces.")
        .def(py::init<Particles<DataLayout::SoA, 2>&>(),
             "Construct an IPS_Simulator object.\n\n"
             "Parameters:\n"
             "    particles (Particles): A collection of particles using a Structure of Arrays (SoA) layout in 2D.")
        .def("init",
             [](IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>& self,
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
             },
             "Initialize the simulator with configurations for inter-particle forces and confinement forces.\n\n"
             "Parameters:\n"
             "    pair_force_config (dict): Dictionary specifying the configuration for the pair force.\n"
             "    confinement_force_config (dict): Dictionary specifying the configuration for the confinement force.\n\n"
             "Example:\n"
             "    >>> simulator.init({'force_type': 'pair'}, {'confinement': 'type1'})")
        .def("integrate",
             &IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>::integrate,
             "Perform a single integration step of the simulation.")
        .def("integrate_n_steps",
             &IPS_Simulator<Particles<DataLayout::SoA, 2>, LeapFrog>::integrate_n_steps,
             "Perform a specified number of integration steps.\n\n"
             "Parameters:\n"
             "    steps (int): The number of integration steps to execute.\n\n"
             "Example:\n"
             "    >>> simulator.integrate_n_steps(100)");

    // IPS_Simulator_Langevin using LangevinSystem and ABOBA
    py::class_<IPS_Simulator<LangevinSystem<2>, ABOBA>>(m, "IPS_Simulator_Langevin",
        "IPS_Simulator_Langevin class for performing simulations using the ABOBA integrator on a LangevinSystem.\n\n"
        "This simulator supports Langevin dynamics with random noise and friction.")
        .def(py::init<LangevinSystem<2>&>(),
             "Construct an IPS_Simulator_Langevin object.\n\n"
             "Parameters:\n"
             "    system (LangevinSystem): A system of particles with Langevin dynamics.")
        .def("init",
             [](IPS_Simulator<LangevinSystem<2>, ABOBA>& self,
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
             },
             "Initialize the Langevin simulator with configurations for inter-particle forces and confinement forces.\n\n"
             "Parameters:\n"
             "    pair_force_config (dict): Configuration for the pair force.\n"
             "    confinement_force_config (dict): Configuration for the confinement force.")
        .def("integrate",
             &IPS_Simulator<LangevinSystem<2>, ABOBA>::integrate,
             "Perform a single integration step using the ABOBA integrator.")
        .def("integrate_n_steps",
             &IPS_Simulator<LangevinSystem<2>, ABOBA>::integrate_n_steps,
             "Perform a specified number of integration steps.");

    // IPS_Simulator_NoseHoover using NoseHooverSystem and NoseHoover
    py::class_<IPS_Simulator<NoseHooverSystem<2>, NoseHoover>>(m, "IPS_Simulator_NoseHoover",
        "IPS_Simulator_NoseHoover class for performing simulations using the NoseHoover integrator on a NoseHooverSystem.\n\n"
        "This simulator implements Nose-Hoover dynamics for thermostatting the system.")
        .def(py::init<NoseHooverSystem<2>&>(),
             "Construct an IPS_Simulator_NoseHoover object.\n\n"
             "Parameters:\n"
             "    system (NoseHooverSystem): A particle system configured for Nose-Hoover dynamics.")
        .def("init",
             [](IPS_Simulator<NoseHooverSystem<2>, NoseHoover>& self,
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
             },
             "Initialize the NoseHoover simulator with configurations for inter-particle and confinement forces.\n\n"
             "Parameters:\n"
             "    pair_force_config (dict): Configuration for the pair force.\n"
             "    confinement_force_config (dict): Configuration for the confinement force.")
        .def("integrate",
             &IPS_Simulator<NoseHooverSystem<2>, NoseHoover>::integrate,
             "Perform a single integration step using the NoseHoover integrator.")
        .def("integrate_n_steps",
             &IPS_Simulator<NoseHooverSystem<2>, NoseHoover>::integrate_n_steps,
             "Perform a specified number of integration steps.");

    // IPS_Simulator_NoseHooverLangevin using NoseHooverLangevinSystem and NoseHooverLangevin
    py::class_<IPS_Simulator<NoseHooverLangevinSystem<2>, NoseHooverLangevin>>(m, "IPS_Simulator_NoseHooverLangevin",
        "IPS_Simulator_NoseHooverLangevin class for performing simulations using the NoseHooverLangevin integrator on a NoseHooverLangevinSystem.\n\n"
        "This simulator combines Nose-Hoover dynamics with Langevin thermostatting.")
        .def(py::init<NoseHooverLangevinSystem<2>&>(),
             "Construct an IPS_Simulator_NoseHooverLangevin object.\n\n"
             "Parameters:\n"
             "    system (NoseHooverLangevinSystem): A particle system configured for Nose-Hoover-Langevin dynamics.")
        .def("init",
             [](IPS_Simulator<NoseHooverLangevinSystem<2>, NoseHooverLangevin>& self,
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
             },
             "Initialize the NoseHooverLangevin simulator with configurations for inter-particle and confinement forces.\n\n"
             "Parameters:\n"
             "    pair_force_config (dict): Configuration for the pair force.\n"
             "    confinement_force_config (dict): Configuration for the confinement force.")
        .def("integrate",
             &IPS_Simulator<NoseHooverLangevinSystem<2>, NoseHooverLangevin>::integrate,
             "Perform a single integration step using the NoseHooverLangevin integrator.")
        .def("integrate_n_steps",
             &IPS_Simulator<NoseHooverLangevinSystem<2>, NoseHooverLangevin>::integrate_n_steps,
             "Perform a specified number of integration steps.");

    // Particles class
    py::class_<Particles<DataLayout::SoA, 2>>(m, "Particles",
        "Particles class for storing particle data in 2 dimensions using a Structure of Arrays (SoA) layout.\n\n"
        "Provides arrays for positions, velocities, and forces.")
        .def(py::init([](size_t n_particles) {
            Particles<DataLayout::SoA, 2> p;
            for (size_t d = 0; d < 2; ++d) {
                p.positions[d].resize(n_particles);
                p.velocities[d].resize(n_particles);
                p.forces[d].resize(n_particles);
            }
            return p;
        }), py::arg("n_particles"),
        "Construct a Particles object.\n\n"
        "Parameters:\n"
        "    n_particles (int): The number of particles.")
        .def_readwrite("positions", &Particles<DataLayout::SoA, 2>::positions,
             "2D vector of particle positions.")
        .def_readwrite("velocities", &Particles<DataLayout::SoA, 2>::velocities,
             "2D vector of particle velocities.")
        .def_readwrite("forces", &Particles<DataLayout::SoA, 2>::forces,
             "2D vector of particle forces.")
        .def("get_positions", &Particles<DataLayout::SoA, 2>::get_positions,
             "Return the positions of all particles in the shape of (2, N).")
        .def("get_velocities", &Particles<DataLayout::SoA, 2>::get_velocities,
             "Return the velocities of all particles in the shape of (2, N).")
        .def("get_forces", &Particles<DataLayout::SoA, 2>::get_forces,
             "Return the forces on all particles in the shape of (2, N).");

    // LangevinSystem class, derived from Particles
    py::class_<LangevinSystem<2>, Particles<DataLayout::SoA, 2>>(m, "LangevinSystem",
        "LangevinSystem class: A particle system with Langevin dynamics.\n\n"
        "Includes additional parameters for friction (gamma) and temperature.")
        .def(py::init([](size_t n_particles, double gamma, double temperature) {
                LangevinSystem<2> p;
                p.gamma = gamma;
                p.temperature = temperature;
                for (size_t d = 0; d < 2; ++d) {
                    p.positions[d].resize(n_particles);
                    p.velocities[d].resize(n_particles);
                    p.forces[d].resize(n_particles);
                }
                return p;
            }), py::arg("n_particles"), py::arg("gamma"), py::arg("temperature"),
            "Construct a LangevinSystem object.\n\n"
            "Parameters:\n"
            "    n_particles (int): Number of particles.\n"
            "    gamma (float): Friction coefficient.\n"
            "    temperature (float): Target temperature.")
        .def_readwrite("gamma", &LangevinSystem<2>::gamma,
             "Friction coefficient for the Langevin dynamics.")
        .def_readwrite("temperature", &LangevinSystem<2>::temperature,
             "Target temperature for the system.");

    // NoseHooverSystem class, derived from Particles
    py::class_<NoseHooverSystem<2>, Particles<DataLayout::SoA, 2>>(m, "NoseHooverSystem",
        "NoseHooverSystem class: A particle system configured for Nose-Hoover dynamics.\n\n"
        "Includes parameters for temperature control, such as temperature, Q, and eta.")
        .def(py::init([](size_t n_particles, double temperature, double Q, double eta) {
                NoseHooverSystem<2> p;
                p.temperature = temperature;
                p.Q = Q;
                p.eta = eta;
                for (size_t d = 0; d < 2; ++d) {
                    p.positions[d].resize(n_particles);
                    p.velocities[d].resize(n_particles);
                    p.forces[d].resize(n_particles);
                }
                return p;
            }), py::arg("n_particles"), py::arg("temperature"), py::arg("Q"), py::arg("eta"),
            "Construct a NoseHooverSystem object.\n\n"
            "Parameters:\n"
            "    n_particles (int): Number of particles.\n"
            "    temperature (float): Target temperature.\n"
            "    Q (float): Nose-Hoover thermostat mass parameter.\n"
            "    eta (float): Initial thermostat variable.")
        .def_readwrite("temperature", &NoseHooverSystem<2>::temperature,
             "Target temperature for the system.")
        .def_readwrite("Q", &NoseHooverSystem<2>::Q,
             "Thermostat mass parameter.")
        .def_readwrite("eta", &NoseHooverSystem<2>::eta,
             "Thermostat variable.");

    // NoseHooverLangevinSystem class, derived from Particles
    py::class_<NoseHooverLangevinSystem<2>, Particles<DataLayout::SoA, 2>>(m, "NoseHooverLangevinSystem",
        "NoseHooverLangevinSystem class: A particle system that combines Nose-Hoover and Langevin dynamics.\n\n"
        "Includes parameters for friction (gamma), temperature control, thermostat mass (Q), and the thermostat variable (eta).")
        .def(py::init([](size_t n_particles, double gamma, double temperature, double Q, double eta) {
                NoseHooverLangevinSystem<2> p;
                p.gamma = gamma;
                p.temperature = temperature;
                p.Q = Q;
                p.eta = eta;
                for (size_t d = 0; d < 2; ++d) {
                    p.positions[d].resize(n_particles);
                    p.velocities[d].resize(n_particles);
                    p.forces[d].resize(n_particles);
                }
                return p;
            }), py::arg("n_particles"), py::arg("gamma"), py::arg("temperature"), py::arg("Q"), py::arg("eta"),
            "Construct a NoseHooverLangevinSystem object.\n\n"
            "Parameters:\n"
            "    n_particles (int): Number of particles.\n"
            "    gamma (float): Friction coefficient for Langevin dynamics.\n"
            "    temperature (float): Target temperature.\n"
            "    Q (float): Nose-Hoover thermostat mass parameter.\n"
            "    eta (float): Initial thermostat variable.")
        .def_readwrite("gamma", &NoseHooverLangevinSystem<2>::gamma,
             "Friction coefficient for the Langevin part of the dynamics.")
        .def_readwrite("temperature", &NoseHooverLangevinSystem<2>::temperature,
             "Target temperature for the system.")
        .def_readwrite("Q", &NoseHooverLangevinSystem<2>::Q,
             "Thermostat mass parameter.")
        .def_readwrite("eta", &NoseHooverLangevinSystem<2>::eta,
             "Thermostat variable.");

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