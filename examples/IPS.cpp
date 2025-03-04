
#include <iostream>
#include "IPS.h"
#include "Integrator.h"
#include <chrono>

// template<size_t d, size_t n>
// auto generate_lattice_positions(double lattice_constant) -> Particles<DataLayout::SoA, d, n>
// {
//     Particles<DataLayout::SoA, d, n> particles;
//     for (size_t i = 0; i < n; i++)
//     {
//         for (size_t j = 0; j < d; j++)
//         {
//             particles.positions[j][i] = i * lattice_constant;
//         }
//     }
//     return particles;
// }
#include <random>

template<size_t d, size_t n>
auto generate_random_positions(double min, double max) -> Particles<DataLayout::SoA, d, n>
{
    Particles<DataLayout::SoA, d, n> particles;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    std::uniform_real_distribution<> dis2(-1, 1);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            particles.positions[j][i] = dis(gen);
            particles.velocities[j][i] = dis2(gen);
            // particles.velocities[j][i] = dis(gen);
        }
    }
    return particles;
}

int main() {

    constexpr size_t n_steps = 1000;
    double step_size = 0.01;
    constexpr size_t num_particles = 500;
    constexpr size_t dim = 2;

    double eps = 1;
    // distance for the pair potential to be zero
    double sigma = 1.0;

    // use Lennard-Jones 12-6 potential:
    // f_r = 4 * eps * (-12 * sigma^12 / r^14 + 6 * sigma^6 / r^8)
    // auto pair_force = [eps, sigma](double r) {
    //     double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    //     double sigma12 = sigma6 * sigma6;

    //     double r6 = r * r * r * r * r * r;
    //     double r8 = r6 * r * r;
    //     return 4 * eps * (-12 * sigma12 / (r8 * r6) + 6 * sigma6 / (r8));
    // };
    const double k = 1.0;
    const double r_0 = 1.0;
    auto pair_force = [k, r_0](double r) {
        return -1 * k * (r - r_0) / r;
    };

    const double coff = 5.0;
    // use a harmonic potential for the confinement force
    auto confinement_force = [coff](const vec<dim>& pos) { 
        return vec<dim>{-coff * pos[0], -coff * pos[1]};
    };

    auto particles = generate_random_positions<dim, num_particles>(-1, 1);
    auto integrator = std::make_unique<LeapFrog>();
    IPS_Simulator<DataLayout::SoA, dim, num_particles, LeapFrog> ips(pair_force, confinement_force,
                                         particles, std::move(integrator)); 

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (size_t i = 0; i < n_steps; i++)
    {
        // output the positions of the particles
        for (size_t j = 0; j < num_particles; j++)
        {
            std::cout << ips.particles.positions[0][j] << " " << ips.particles.positions[1][j] << " " <<
            ips.particles.velocities[0][j] << " " << ips.particles.velocities[1][j] << " " <<
            ips.particles.forces[0][j] << " " << ips.particles.forces[1][j] << std::endl;
        }
        ips.integrate(step_size);
        std::cout << std::endl;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    auto one_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    // std::cout << "# Time difference = " << one_time << "[ms]" << std::endl;

    return 0;

}