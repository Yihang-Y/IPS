#pragma once
#include "Integrator.h"
#include "Potential.h"
#include <memory>

template<DataLayout Layout, size_t Dim, typename Integrator>
class IPS_Simulator
{
public:
    using pair_force_callable = std::function<double(double)>;
    using confinement_force_callable = std::function<vec<Dim>(const vec<Dim>&)>;
    /**
     * @brief Construct a new ips simulator object,
     * use pair potential, confinement potential, and initial positions and velocities to initialize the system.
     * 
     */
    IPS_Simulator(pair_force_callable _pair_force, confinement_force_callable _confinement_force,
                  Particles<Layout, Dim>& _particles,
                  Integrator _integrator)
        : pair_force(_pair_force), confinement_force(_confinement_force), particles(_particles), integrator(std::move(_integrator))
    {}

    IPS_Simulator(Particles<Layout, Dim>& _particles)
        : particles(_particles)
    {
        pair_force = Spring(1.0, 1.0);
        confinement_force = RadialConfinement<Dim>(2.0);
        integrator = LeapFrog();
    }

    // template<typename Config1, typename Config2>
    // IPS_Simulator(Particles<Layout, Dim>& _particles, const Config1& config1, const Config2& config2)
    //     : particles(_particles)
    // {
    //     pair_force = make_potential(config1);
    //     confinement_force = make_potential(config2);
    //     integrator = LeapFrog();
    // }

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
    Particles<Layout, Dim>& particles;
    Integrator integrator;

};