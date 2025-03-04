#pragma once
#include "vec.h"

class LennardJones {
public:
    LennardJones() : eps(1.0), sigma(1.0) {}
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
