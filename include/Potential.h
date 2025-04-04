#pragma once
#include "vec.h"
#include "utils.h"
#include <map>
#include <variant>
#include <memory>
#include <iostream>

class LennardJones {
public:
    LennardJones() : eps(1.0), sigma(1.0) {}
    LennardJones(double eps, double sigma) : eps(eps), sigma(sigma) {}
    double operator()(double r) const {
        double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sigma12 = sigma6 * sigma6;

        double r6 = r * r * r * r * r * r;
        double r8 = r6 * r * r;
        double r14 = r6 * r8;

        return 24 * eps * ( -2 * sigma12 / r14 + sigma6 / r8);
    }

    REFLECT(eps, sigma);
private:
    double eps;
    double sigma;
};
REGISTER_PAIR_FORCE(LennardJones);

class Spring {
public:
    Spring() : k(1.0), r_0(1.0) {}
    Spring(double k, double r_0) : k(k), r_0(r_0) {}
    double operator()(double r) const {
        return k * (r - r_0) / r;
    }

    REFLECT(k, r_0);
private:
    double k;
    double r_0;
};
REGISTER_PAIR_FORCE(Spring);

template<size_t dim>
class HarmonicConfinement {
public:
    HarmonicConfinement() : coff(1.0) {}
    HarmonicConfinement(double coff) : coff(coff) {}
    vec<dim> operator()(const vec<dim>& pos) const {
        return vec<dim>{-coff * pos[0], -coff * pos[1]};
    }

    REFLECT(coff);
private:
    double coff;
};
// REGISTER_CONFINEMENT_FORCE(HarmonicConfinement(dim));

template<size_t dim>
class RadialConfinement {
public:
    RadialConfinement() : rad(1.0){}
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

    REFLECT(rad);
private:
    double rad;
};
// REGISTER_CONFINEMENT_FORCE(RadialConfinement(dim));

auto make_pair_force(const Config& config) -> std::function<double(double)> {
    auto type = config.get<std::string>("type", "LennardJones");
    auto it = get_PFF_map().find(type);
    if (it != get_PFF_map().end()) {
        return it->second(config);
    }
    throw std::invalid_argument("Unknown pair force type: " + type);
}

auto make_confinement_force(const Config& config) -> std::function<vec<2>(const vec<2>&)> {
    auto type = config.get<std::string>("type", "Harmonic");
    if (type == "Harmonic") {
        HarmonicConfinement<2> instance;
        instance.from_config(config);
        return instance;
    } else if (type == "Radial") {
        RadialConfinement<2> instance;
        instance.from_config(config);
        return instance;
    } else {
        throw std::invalid_argument("Unknown confinement force type: " + type);
    }
}