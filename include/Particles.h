#include "vec.h"

enum class DataLayout
{
    SoA,
    AoS
};

/**
 * @brief Use SoA to store the particle system, which is more cache-friendly and might be helpful for the design of the integrator.
 * 
 * @tparam Dim 
 * @tparam NumParticles 
 */
template<DataLayout Layout, size_t Dim, size_t NumParticles>
struct Particles;

template<size_t Dim, size_t NumParticles>
struct Particles<DataLayout::SoA, Dim, NumParticles>
{
    std::array<std::array<double, NumParticles>, Dim> positions = {};
    std::array<std::array<double, NumParticles>, Dim> velocities = {};
    std::array<std::array<double, NumParticles>, Dim> forces = {};

    vec<Dim> get_position(size_t i) const
    {
        vec<Dim> pos;
        for (size_t d = 0; d < Dim; d++)
        {
            pos[d] = positions[d][i];
        }
        return pos;
    }
};
