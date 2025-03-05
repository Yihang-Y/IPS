#pragma once
#include "vec.h"

struct LJConfig {
    double eps;
    double sigma;
};

struct SpringConfig {
    double k;
    double r_0;
};

struct HarmonicConfinementConfig {
    double coff;
};

struct RadialConfinementConfig {
    double rad;
};

class LennardJones {
public:
    LennardJones() : eps(1.0), sigma(1.0) {}
    LennardJones(const LJConfig& config) : eps(config.eps), sigma(config.sigma) {}
    LennardJones(double eps, double sigma) : eps(eps), sigma(sigma) {}
    double operator()(double r) const {
        double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sigma12 = sigma6 * sigma6;

        double r6 = r * r * r * r * r * r;
        double r12 = r6 * r6;

        return 4 * eps * (sigma12 / r12 - sigma6 / r6);
    }

private:
    double eps;
    double sigma;
};

class Spring {
public:
    Spring() : k(1.0), r_0(1.0) {}
    Spring(const SpringConfig& config) : k(config.k), r_0(config.r_0) {}
    Spring(double k, double r_0) : k(k), r_0(r_0) {}
    double operator()(double r) const {
        return -1 * k * (r - r_0) / r;
    }
private:
    double k;
    double r_0;
};

template<size_t dim>
class HarmonicConfinement {
public:
    HarmonicConfinement() : coff(1.0) {}
    HarmonicConfinement(const HarmonicConfinementConfig& config) : coff(config.coff) {}
    HarmonicConfinement(double coff) : coff(coff) {}
    vec<dim> operator()(const vec<dim>& pos) const {
        return vec<dim>{-coff * pos[0], -coff * pos[1]};
    }
private:
    double coff;
};

template<size_t dim>
class RadialConfinement {
public:
    RadialConfinement() : rad(1.0){}
    RadialConfinement(const RadialConfinementConfig& config) : rad(config.rad){}
    RadialConfinement(double rad) : rad(rad){}
    vec<dim> operator()(const vec<dim>& pos) const {
        double r = 0;
        for (size_t d = 0; d < dim; d++){
            r += pos[d] * pos[d];
        }
        r = std::sqrt(r);
        if (r > rad)
        {
            double coff = (r - rad) / (1 - (r - rad));
            return vec<dim>{-coff * pos[0], -coff * pos[1]};
        }
        else 
        {
            return vec<dim>{0, 0};
        }
    }

private:
    double rad;
};


template<typename Config>
struct config_traits;

template<>
struct config_traits<LJConfig> {
    using potential_type = LennardJones;
};

template<>
struct config_traits<SpringConfig> {
    using potential_type = Spring;
};

template<>
struct config_traits<HarmonicConfinementConfig> {
    using potential_type = HarmonicConfinement<2>;
};

template<>
struct config_traits<RadialConfinementConfig> {
    using potential_type = RadialConfinement<2>;
};

template<typename Config>
auto make_potential(const Config& config) {
    using potential_type = typename config_traits<Config>::potential_type;
    // FIXME: maybe we can use other ways, don't need to add a new init function in each potential class
    return potential_type(config);
}