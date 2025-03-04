#pragma once
#include "Integrator.h"
#include <memory>

template<DataLayout Layout, size_t Dim, size_t NumParticles, typename Integrator>
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
                  Particles<Layout, Dim, NumParticles> _particles,
                  std::unique_ptr<Integrator> _integrator)
        : pair_force(_pair_force), confinement_force(_confinement_force), particles(_particles), integrator(std::move(_integrator))
    {}

    void integrate(double step_size) {
        integrator->integrate(particles, pair_force, confinement_force, step_size);
    }

public:
    pair_force_callable pair_force;
    confinement_force_callable confinement_force;
    Particles<Layout, Dim, NumParticles> particles;
    std::unique_ptr<Integrator> integrator;

};