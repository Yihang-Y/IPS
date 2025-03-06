
#include <iostream>
#include <string>
#include "IPS.h"
#include "Integrator.h"
#include <chrono>
#include <random>
#include "Potential.h"

// arguments n_steps, step_size, num_particles, dim, output_interval
int main(int argc, char const *argv[])
{
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " n_steps step_size num_particles dim output_interval" << std::endl;
        return 1;
    }

    size_t n_steps = std::stoi(argv[1]);
    double step_size = std::stod(argv[2]);
    size_t num_particles = std::stoi(argv[3]);
    size_t dim_ = std::stoi(argv[4]);
    size_t output_interval = std::stoi(argv[5]);

    if(dim_ == 2)
    {
        // std::cerr << "Only support 2D for now." << std::endl;
        // return 1;
        constexpr size_t dim = 2;

        auto pair_force = Spring(1.0, 1.0);
        auto confinement_force = RadialConfinement<dim>(2.0);

        auto particles = generate_random_init<DataLayout::SoA, dim>(num_particles, -1, 1, -2, 2);
        auto integrator = LeapFrog();
        IPS_Simulator<DataLayout::SoA, dim, LeapFrog> ips(pair_force, confinement_force,
                                            particles, std::move(integrator)); 

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (size_t i = 0; i < n_steps; i++)
        {
            if( i % output_interval == 0)
            {
                for (size_t j = 0; j < num_particles; j++)
                {
                    std::cout << ips.particles.positions[0][j] << " " << ips.particles.positions[1][j] << " " <<
                    ips.particles.velocities[0][j] << " " << ips.particles.velocities[1][j] << " " <<
                    ips.particles.forces[0][j] << " " << ips.particles.forces[1][j] << std::endl;
                }
                std::cout << std::endl;
            }
            ips.integrate(step_size);
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto one_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        std::cout << "# Time difference = " << one_time << "[ms]" << std::endl;

        return 0;
    }
    else {
        constexpr size_t dim = 3;

        auto pair_force = LennardJones(1.0, 1.0);
        auto confinement_force = RadialConfinement<dim>(2.0);

        auto particles = generate_random_init<DataLayout::SoA, dim>(num_particles, -1, 1, -2, 2);
        auto integrator = LeapFrog();
        IPS_Simulator<DataLayout::SoA, dim, LeapFrog> ips(pair_force, confinement_force,
                                            particles, std::move(integrator));

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (size_t i = 0; i < n_steps; i++)
        {
            if( i % output_interval == 0)
            {
                for (size_t j = 0; j < num_particles; j++)
                {
                    for (size_t d = 0; d < dim; d++)
                    {
                        std::cout << ips.particles.positions[d][j] << " ";
                    }
                    for (size_t d = 0; d < dim; d++)
                    {
                        std::cout << ips.particles.velocities[d][j] << " ";
                    }
                    for (size_t d = 0; d < dim-1; d++)
                    {
                        std::cout << ips.particles.forces[d][j] << " ";
                    }
                    // we do this separately so that there is no " " at the end of the line, which causes problems when reading the particle data
                    std::cout << ips.particles.forces[2][j]; // TODO: make this more generalisable
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            ips.integrate(step_size);
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto one_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        std::cout << "# Time difference = " << one_time << "[ms]" << std::endl;                       

    }

}