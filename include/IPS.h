#pragma once
#include "Integrator.h"
#include "Potential.h"
#include <memory>

template<typename System, typename Integrator>
class IPS_Simulator
{
public:
    constexpr static size_t Dim = System::DimVal;
    // using Layout = System::Layout;
    
    using pair_force_callable = std::function<double(double)>;
    using confinement_force_callable = std::function<vec<Dim>(const vec<Dim>&)>;
    /**
     * @brief Construct a new ips simulator object,
     * use pair potential, confinement potential, and initial positions and velocities to initialize the system.
     * 
     */
    IPS_Simulator(System& _particles)
        : particles(_particles)
    {
        pair_force = Spring(1.0, 1.0);
        confinement_force = RadialConfinement<Dim>(2.0);
        integrator = Integrator();
    }

    void integrate(double step_size) {
        integrator.integrate(particles, pair_force, confinement_force, step_size);
    }

    void integrate_n_steps(double step_size, size_t n_steps) {
        for (size_t i = 0; i < n_steps; i++)
        {
            integrate(step_size);
        }
    }

public:
    pair_force_callable pair_force;
    confinement_force_callable confinement_force;
    System& particles;
    Integrator integrator;
};