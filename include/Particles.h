#pragma once
#include "vec.h"
#include <vector>
#include <array>
#include <random>
#include <stdexcept>

#include "utils.h"

enum class DataLayout
{
    SoA,
    AoS
};

/**
 * @brief Stores info on positions, velocities, and for
 * Uses SoA to store the particle system, which is more cache-friendly and might be helpful for the design of the integrator.
 * 
 * @tparam Dim 
 */
template<DataLayout Layout, size_t Dim>
struct Particles;

template<size_t Dim>
struct Particles<DataLayout::SoA, Dim>
{
    constexpr static size_t DimVal = Dim;
    constexpr static DataLayout Layout = DataLayout::SoA;

    std::array<std::vector<double>, Dim> positions = {};
    std::array<std::vector<double>, Dim> velocities = {};
    std::array<std::vector<double>, Dim> forces = {};

    // function to get the position of the i'th particle.
    vec<Dim> get_position(size_t i) const
    {
        size_t NumParticles = positions[0].size();
        if (i >= NumParticles) {
            throw std::invalid_argument("Requested particle " + std::to_string(i+1) + " but there are only " + std::to_string(NumParticles) + " particles. Note zero indexing.");
        } 
        
        vec<Dim> pos;
        for (size_t d = 0; d < Dim; d++)
        {
            pos[d] = positions[d][i];
        }
        return pos;
    }

// NOTE: the code below is for pybind11, do not call them in C++ code
#ifdef USE_PYBIND11
    py::list get_positions() {
        py::list result;
        for (size_t d = 0; d < Dim; ++d) {
            result.append(vector_to_array(positions[d]));
        }
        return result;
    }

    py::list get_velocities() {
        py::list result;
        for (size_t d = 0; d < Dim; ++d) {
            result.append(vector_to_array(velocities[d]));
        }
        return result;
    }

    py::list get_forces() {
        py::list result;
        for (size_t d = 0; d < Dim; ++d) {
            result.append(vector_to_array(forces[d]));
        }
        return result;
    }
#endif
};


template<size_t Dim>
struct LangevinSystem : public Particles<DataLayout::SoA, Dim>
{
    double gamma = 1.0;
    double temperature = 1.0;

#ifdef USE_PYBIND11
    double* get_gamma() {
        return &gamma;
    }

    double* get_temperature() {
        return &temperature;
    }
#endif
};

template<DataLayout Layout, size_t Dim>

/**
 * @brief Generates initial values for particle positions and velocities by uniformly sampling between
 pmin and pmax for position and vmin and vmax for velocity. Initial forces set to 0.
 * 
 * @param num number of particles
 * @param pmin minimum position value (in each coordinate)
 * @param pmax max position value
 * @param vmin min velocity value
 * @param vmax max velocity value
 * @return auto Particles object
 */
auto generate_random_init(size_t num, double pmin, double pmax, double vmin, double vmax)
{
    Particles<Layout, Dim> particles;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_p(pmin, pmax);
    std::uniform_real_distribution<> dis_v(vmin, vmax);
    for (size_t d = 0; d < Dim; d++) {
        particles.positions[d].reserve(num);
        particles.velocities[d].reserve(num);
        particles.forces[d].reserve(num);
    }

    for (size_t i = 0; i < num; i++)
    {
        for (size_t j = 0; j < Dim; j++)
        {
            particles.positions[j].push_back(dis_p(gen));
            particles.velocities[j].push_back(dis_v(gen));
            particles.forces[j].push_back(0);
        }
    }
    return particles;
}