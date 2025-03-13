#pragma once
#include <array>
#include <functional>
#include <cmath>
#include "Particles.h"
#include <iostream>


template<typename Derived>
class IntegratorBase {
public:
    using pair_force_callable = std::function<double(double)>;
    template <size_t Dim>
    using confinement_force_callable = std::function<vec<Dim>(const vec<Dim>&)>;

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

        // 3. use new positions to calculate the forces
        // std::array<std::array<double, N>, Dim> forces = {};
        std::array<std::vector<double>, Dim> forces = {};
        for (size_t d = 0; d < Dim; d++)
        {
            forces[d].resize(N);
        }

        // 3.1 caculate the confinement forces
        for (size_t i = 0; i < N; i++)
        {
            auto pos = particles.get_position(i);
            auto conF = confinement_force(pos);
            for (size_t d = 0; d < Dim; d++)
            {
                forces[d][i] += conF[d];
            }
        }

        // 3.2. calculate the pair forces f_ij = - (d\phi(r_ij) / r_ij) * (q_i - q_j)

        // 3.2.1 calculate the scalar coff: f_r = - (d\phi(r_ij) / r_ij)
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
        // 3.2.2 calculate the pair forces f_ij = f_r * (q_i - q_j) and update the forces
        for (size_t d = 0; d < Dim; d++)
        {
            for (size_t i = 0; i < N; i++)
            {
                auto pos_i = particles.get_position(i);
                for (size_t j = i + 1; j < N; j++)
                {
                    auto pos_j = particles.get_position(j);
                    double f_ij = pair_force_r[i][j] * (pos_i[d] - pos_j[d]);
                    forces[d][i] += f_ij;
                    forces[d][j] -= f_ij;
                }
            }
        }

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
        // std::array<std::array<double, N>, Dim> forces = {};
        // std::array<std::vector<double>, Dim> forces = {};
        for (size_t d = 0; d < Dim; d++)
        {
            std::fill(particles.forces[d].begin(), particles.forces[d].end(), 0.0);
        }

        // 2.1 caculate the confinement forces
        for (size_t i = 0; i < N; i++)
        {
            auto pos = particles.get_position(i);
            auto conF = confinement_force(pos);
            for (size_t d = 0; d < Dim; d++)
            {
                particles.forces[d][i] += conF[d];
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
                    particles.forces[d][i] -= f_ij;
                    particles.forces[d][j] += f_ij;
                }
            }
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