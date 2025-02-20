/**
 * @file OverdampedLangevin.h
 * @author Yihang Yin (yihangyin@hotmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-02-20
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once
#include <cstddef>
#include <array>
#include <cmath>
#include <random>
#include <functional>
#include <vector>


template<int Dim>
struct position
{
    std::array<double, Dim> x;
};

template<int Dim>
class OverdampedLangevin
{
public:
    using force_callable = std::function<double(position<Dim>&)>;
    OverdampedLangevin(double _kbT, force_callable _force, position<Dim> _init): kbT(_kbT), force(_force), current_position(_init){
        generator.seed(std::random_device()());
    }
    ~OverdampedLangevin(){
        
    }
    void eulerMaruyamaStep(double step_size){
        // update the position
        for (size_t i = 0; i < Dim; i++)
        {
            current_position.x[i] += force(current_position) * step_size + std::sqrt(2 * kbT * step_size) * std::normal_distribution<double>(0, 1)(generator);
        }
    }
    const position<Dim>& getCurrentPosition(){
        return current_position;
    }
    std::vector<std::array<double, Dim>> getTrajectory(size_t n_steps, double step_size){
        std::vector<std::array<double, Dim>> trajectory;
        for (size_t i = 0; i < n_steps; i++)
        {
            eulerMaruyamaStep(step_size);
            trajectory.push_back(current_position.x);
        }
        return trajectory;
    }

private:
    double kbT;
    force_callable force;
    position<Dim> current_position;
    std::mt19937 generator;
};