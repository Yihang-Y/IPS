#pragma once
#include "vec.h"
#include "utils.h"
#include <map>
#include <variant>
#include <memory>
#include <iostream>

using ConfigValue = std::variant<int, double, std::string, bool>;
using ConfigMap = std::map<std::string, ConfigValue>;

struct Config {
    ConfigMap config_map;
    Config() = default;

    ConfigValue& operator[](const std::string& key) {
        return config_map[key];
    }

    const ConfigValue& operator[](const std::string& key) const {
        return config_map.at(key);
    }

    template<typename T>
    T get(const std::string& key, const T& default_value = T{}) const {
        auto it = config_map.find(key);
        if (it == config_map.end()) {
            return default_value;
        }
        return std::get<T>(it->second);
    }
};

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

class Spring {
public:
    Spring() : k(1.0), r_0(1.0) {}
    Spring(double k, double r_0) : k(k), r_0(r_0) {}
    double operator()(double r) const {
        return -1 * k * (r - r_0) / r;
    }

    REFLECT(k, r_0);
private:
    double k;
    double r_0;
};

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

// FIXME: should use more reflection to make the code more generic, so we don't need to write the same code again and again
// but for now, it works.
auto make_pair_force(const Config& config) -> std::function<double(double)> {
    auto type = config.get<std::string>("type", "LennardJones");
    if (type == "LennardJones") {
        LennardJones instance;
        instance.from_config(config);
        return instance;
    } else if (type == "Spring") {
        Spring instance;
        instance.from_config(config);
        return instance;
    } else {
        throw std::invalid_argument("Unknown pair force type: " + type);
    }
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