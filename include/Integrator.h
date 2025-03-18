#pragma once
#include <array>
#include <functional>
#include <cmath>
#include "Particles.h"
#include <iostream>

using pair_force_callable = std::function<double(double)>;
template <size_t Dim>
using confinement_force_callable = std::function<vec<Dim>(const vec<Dim>&)>;

template<typename Derived>
class IntegratorBase {
public:

    template<typename Particles>
    void integrate(Particles& particles, 
                    const pair_force_callable& pair_force, 
                    const confinement_force_callable<Particles::DimVal>& confinement_force,
                    double step_size){
        if constexpr (Particles::Layout == DataLayout::SoA)
        {
            static_cast<Derived*>(this)->integrate_on_soa(particles, pair_force, confinement_force, step_size);
        }
        else if constexpr (Particles::Layout == DataLayout::AoS)
        {
            static_cast<Derived*>(this)->integrate_on_aos(particles, pair_force, confinement_force, step_size);
        }
    }
};

template<typename Particles>
static inline void compute_force(Particles& particles,
                            std::array<std::vector<double>, Particles::DimVal>& forces,
                            const pair_force_callable& pair_force,
                            const confinement_force_callable<Particles::DimVal>& confinement_force) {
    const size_t N = forces[0].size();
    // 1. use new positions to calculate the forces
    constexpr size_t Dim = Particles::DimVal;

    // 2.1 caculate the confinement forces
    for (size_t i = 0; i < N; i++)
    {
        auto pos = particles.get_position(i);
        auto conF = confinement_force(pos);
        for (size_t d = 0; d < Dim; d++)
        {
            forces[d][i] += conF[d];
        }
    }

    // 2.2. calculate the pair forces f_ij = - (d\phi(r_ij) / r_ij) * (q_i - q_j)

    // 2.2.1 calculate the scalar coff: f_r = - (d\phi(r_ij) / r_ij)
    // store them into a 2D array: pair_force_r[N][N]
    std::vector<std::vector<double>> pair_force_r(N, std::vector<double>(N, 0));
    for (size_t i = 0; i < N; i++)
    {
        auto pos_i = particles.get_position(i);
        for (size_t j = i + 1; j < N; j++)
        {
            auto pos_j = particles.get_position(j);
            double r = 0;
            for (size_t d = 0; d < Dim; d++)
            {
                r += (pos_i[d] - pos_j[d]) * (pos_i[d] - pos_j[d]);
            }
            r = std::sqrt(r);
            
            double pairF = pair_force(r);
            pair_force_r[i][j] = pairF;
            // pair_force_r[j][i] = pairF;
        }
    }
    // 2.2.2 calculate the pair forces f_ij = f_r * (q_i - q_j) and update the forces
    for (size_t d = 0; d < Dim; d++)
    {
        for (size_t i = 0; i < N; i++)
        {
            auto pos_i = particles.get_position(i);
            for (size_t j = i + 1; j < N; j++)
            {
                auto pos_j = particles.get_position(j);
                double f_ij = pair_force_r[i][j] * (pos_i[d] - pos_j[d]);
                forces[d][i] -= f_ij;
                forces[d][j] += f_ij;
            }
        }
    }
}

class LeapFrog : public IntegratorBase<LeapFrog> {
public:

    template<typename Particles>
    void integrate_on_soa(Particles& particles,
                            const pair_force_callable& pair_force,
                            const confinement_force_callable<Particles::DimVal>& confinement_force,
                            double step_size){
        // does one time step
        
        const size_t N = particles.positions[0].size();
        constexpr size_t Dim = Particles::DimVal;

        // leap frog integration
        // 1. update the half step velocities, use the forces from last step
        double half_step = 0.5 * step_size;
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += half_step * particles.forces[d][i];
            }
        }
        // 2. update the positions
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.positions[d][i] += step_size * particles.velocities[d][i];
            }
        }
        // 3. use new positions to calculate the 
        std::array<std::vector<double>, Dim> forces = {};
        for (size_t d = 0; d < Dim; d++)
        {
            forces[d].resize(N);
        }
        compute_force(particles, forces, pair_force, confinement_force);


        // 4. update the velocities, using new forces
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += half_step * forces[d][i];
            }
        }

        // 5. update the forces
        particles.forces = forces;
    }

    template<typename Particles>
    void integrate_on_aos(Particles& particles,
                            const pair_force_callable& pair_force,
                            const confinement_force_callable<Particles::Dim>& confinement_force,
                            double step_size){}
       
};

class ABOBA : public IntegratorBase<ABOBA> {
public:

    /**
     * @brief Implement the ABOBA integrator for the system
     * 
     * @tparam Dim 
     * @param particles: with positions, velocities, forces in SOA and temperature 
     * @param pair_force: callable object, pair_force(r_ij) return the force between particle i and j, need to * (q_i - q_j) to get the real force
     * @param confinement_force: callable object, confinement_force(pos) return the force on the particle at pos
     * @param step_size 
     */
    template<typename Particles>
    void integrate_on_soa(Particles& particles,
                            const pair_force_callable& pair_force, 
                            const confinement_force_callable<Particles::DimVal>& confinement_force,
                            double step_size){
        // does one time step
        constexpr size_t Dim = Particles::DimVal;
        const size_t N = particles.positions[0].size();

        double half_step = 0.5 * step_size;

        // ABOBA integration
        // 1. A-step update the half step positions, use the velocities from last step
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.positions[d][i] += half_step * particles.velocities[d][i];
            }
        }

        // 2. use new positions to calculate the forces
        std::array<std::vector<double>, Dim> forces = {};
        for (size_t d = 0; d < Dim; d++)
        {
            forces[d].resize(N);
            std::fill(forces[d].begin(), forces[d].end(), 0);
        }

        compute_force(particles, forces, pair_force, confinement_force);
        for (size_t d = 0; d < Dim; d++)
        {
            std::swap(particles.forces[d], forces[d]);
        }

        // 3. B-step update the velocities, use the forces from last step
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += step_size * particles.forces[d][i];
            }
        }

        // 4. O-step update the positions, use the velocities from the B-step
        static std::mt19937 gen(std::random_device{}());
        static std::normal_distribution<double> normal_dist(0, 1);
        double c1 = std::exp(-step_size * particles.gamma);
        double c2 = std::sqrt( (1 - exp(-2 * step_size * particles.gamma)) * particles.temperature);

        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                double noise = normal_dist(gen);
                particles.velocities[d][i] = c1 * particles.velocities[d][i] + c2 * noise;
            }
        }

        // 5. B-step update the velocities, use the forces from last step
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += half_step * particles.forces[d][i];
            }
        }

        // 6. A-step update the half step positions, use the velocities from last step
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.positions[d][i] += half_step * particles.velocities[d][i];
            }
        }
    }

    template<typename Particles>
    void integrate_on_aos(Particles& particles,
                            const pair_force_callable& pair_force, 
                            const confinement_force_callable<Particles::DimVal>& confinement_force,
                            double step_size){}
};

class NoseHoover : public IntegratorBase<NoseHoover> {
public:
    template<typename Particles>
    void integrate_on_soa(Particles& particles,
                            const pair_force_callable& pair_force, 
                            const confinement_force_callable<Particles::DimVal>& confinement_force,
                            double step_size){
        constexpr size_t Dim = Particles::DimVal;
        const size_t N = particles.positions[0].size();
        
        const double half_step = 0.5 * step_size;
        // 1. A step, update the positions
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.positions[d][i] += half_step * particles.velocities[d][i];
            }
        }

        // 2. calculate the forces
        std::array<std::vector<double>, Dim> forces = {};
        for (size_t d = 0; d < Dim; d++)
        {
            forces[d].resize(N);
            std::fill(forces[d].begin(), forces[d].end(), 0);
        }
        compute_force(particles, forces, pair_force, confinement_force);
        for (size_t d = 0; d < Dim; d++)
        {
            std::swap(particles.forces[d], forces[d]);
        }

        // 3. B step, update the velocities
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += step_size * particles.forces[d][i];
            }
        }

        // 4. C step, uodate the thermostat, using p_0 * exp(-\eta * dt / 2)?
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] *= std::exp(-particles.eta * half_step);
            }
        }

        // 5. D step, update the thermostat variable eta

        // 5.1 calculate the kinetic energy
        double kinetic_energy = 0;
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                kinetic_energy += 0.5 * particles.velocities[d][i] * particles.velocities[d][i];
            }
        }
        particles.eta += step_size * (kinetic_energy - Dim * particles.temperature) / particles.Q;

        // 6. C step, uodate the thermostat, using p_0 * exp(-\eta * dt / 2)?
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] *= std::exp(-particles.eta * half_step);
            }
        }

        // 7. B step, update the velocities
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.velocities[d][i] += half_step * particles.forces[d][i];
            }
        }

        // 8. A step, update the positions
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                particles.positions[d][i] += half_step * particles.velocities[d][i];
            }
        }
    }

    template<typename Particles>
    void integrate_on_aos(Particles& particles,
                            const pair_force_callable& pair_force, 
                            const confinement_force_callable<Particles::DimVal>& confinement_force,
                            double step_size){}
};

template<size_t Dim>
static inline double LJ_potential(double r, double epsilon, double sigma) {
    double r6 = std::pow(sigma / r, 6);
    return 4 * epsilon * (r6 * r6 - r6);
}

template<size_t Dim>
static inline double total_LJ_potential(const std::array<std::vector<double>, Dim>& positions, double epsilon, double sigma) {
    const size_t N = positions[0].size();
    double total_potential = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i + 1; j < N; j++)
        {
            double r = 0;
            for (size_t d = 0; d < Dim; d++)
            {
                r += (positions[d][i] - positions[d][j]) * (positions[d][i] - positions[d][j]);
            }
            r = std::sqrt(r);
            total_potential += LJ_potential(r, epsilon, sigma);
        }
    }
    return total_potential;
}

// class NosePoincare : public IntegratorBase<NosePoincare> {
// public:
//     template<typename Particles>
//     void integrate_on_soa(Particles& particles,
//                             const pair_force_callable& pair_force, 
//                             const confinement_force_callable<Particles::DimVal>& confinement_force,
//                             double step_size){
//         constexpr size_t Dim = Particles::DimVal;
//         const size_t N = particles.positions[0].size();

//         const double half_step = 0.5 * step_size;
//         // 1. update velocities, p_1/2 = p_0 - dt / 2 * s_0 * F(q_0)

//         // 1.1 calculate the forces
//         std::array<std::vector<double>, Dim> forces = {};
//         for (size_t d = 0; d < Dim; d++)
//         {
//             forces[d].resize(N);
//             std::fill(forces[d].begin(), forces[d].end(), 0);
//         }
//         compute_force(particles, forces, pair_force, confinement_force);
//         for (size_t d = 0; d < Dim; d++)
//         {
//             std::swap(particles.forces[d], forces[d]);
//         }

//         // 1.2 update the velocities
//         for (size_t d = 0; d < Dim; d++)
//         {
//             for (size_t i = 0; i < N; i++)
//             {
//                 particles.velocities[d][i] -= half_step * particles.s * particles.forces[d][i];
//             }
//         }

//         // 2. update pi_1/2, using explicit solution to quadratic equation

//         // 2.1 calculate the kinetic energy, using 1/2 * \sum_i m_i v_i^2 / s^2
//         double kinetic_energy = 0;
//         for (size_t d = 0; d < Dim; d++)
//         {
//             for (size_t i = 0; i < N; i++)
//             {
//                 kinetic_energy += 0.5 * particles.velocities[d][i] * particles.velocities[d][i];
//             }
//         }
//         kinetic_energy /= particles.s * particles.s;

//         // 2.2 calculate the potential energy, 
//         double eplison = pair_force.target<LennardJones>()->eps;
//         double sigma = pair_force.target<LennardJones>()->sigma;
//         double potential_energy = total_LJ_potential(particles.positions, eplison, sigma);

//         // 2.3 calculate the H_0, using 1/2 * \sum_i m_i v_i^2 + U(q) + 1/2 * s^2 pi^2
//         // static double H_0 = 

//         double C = half_step * (Dim * particles.temperature * (1 + log(particles.s));
//                                 - kinetic_energy + potential_energy - H_0) - particles.pi;


//     }

//     template<typename Particles>
//     void integrate_on_aos(Particles& particles,
//                             const pair_force_callable& pair_force, 
//                             const confinement_force_callable<Particles::DimVal>& confinement_force,
//                             double step_size){}
// };