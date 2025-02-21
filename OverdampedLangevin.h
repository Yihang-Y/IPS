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
#include <thread>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <memory>

namespace py = pybind11;

template<int Dim>
py::array convert_from_shared_ptr(std::shared_ptr<std::vector<std::array<double, Dim>>> ptr)
{
    double* data_ptr = &((*ptr)[0][0]);
    std::vector<ssize_t> shape = { static_cast<ssize_t>(ptr->size()), static_cast<ssize_t>(Dim) };
    std::vector<ssize_t> strides = { static_cast<ssize_t>(Dim * sizeof(double)), static_cast<ssize_t>(sizeof(double)) };
    return py::array_t<double>(py::buffer_info(
        data_ptr,                            // 数据指针
        sizeof(double),                      // 每个元素大小
        py::format_descriptor<double>::format(), // 数据格式
        2,                                   // 维度
        shape,                               // 形状
        strides                              // 步长
    ));
}

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
    // 修改：返回 NumPy 数组对象
    auto getTrajectory(size_t n_steps, double step_size)
    {
        std::vector<std::array<double, Dim>> trajectory;
        trajectory.reserve(n_steps);
        for (size_t i = 0; i < n_steps; i++)
        {
            eulerMaruyamaStep(step_size);
            trajectory.push_back(current_position.x);
        }

        return std::make_shared<std::vector<std::array<double, Dim>>>(trajectory);
    }

    py::array getTrajectoryPy(size_t n_steps, double step_size)
    {
        return convert_from_shared_ptr<Dim>(getTrajectory(n_steps, step_size));
    }
    

private:
    double kbT;
    force_callable force;
    position<Dim> current_position;
    std::mt19937 generator;
};



template<int Dim>
class BatchedOverdampedLangevin
{
public:
    using force_callable = std::function<double(position<Dim>&)>;
    
    BatchedOverdampedLangevin(double _kbT, force_callable _force, const std::vector<position<Dim>>& _init)
        : kbT(_kbT), force(_force), initial_positions(_init)
    {
    }

    ~BatchedOverdampedLangevin() = default;


    auto getBatchedTrajectory(size_t n_steps, double step_size)
    {
        size_t batch_size = initial_positions.size();
        std::vector<std::shared_ptr<std::vector<std::array<double, Dim>>>> trajectories;
        trajectories.resize(batch_size);
        
        size_t n_threads = std::thread::hardware_concurrency();
        if (n_threads == 0) {
            n_threads = 1; 
        }
        std::cout << "system has " << n_threads << " threads" << std::endl;

        size_t total_size = batch_size;
        size_t chunk = total_size / n_threads;
        size_t remainder = total_size % n_threads;

        std::vector<std::thread> threads;
        threads.reserve(n_threads);

        size_t start_index = 0;
        for (size_t i = 0; i < n_threads; i++)
        {
            size_t work_size = chunk + ((i < remainder) ? 1 : 0);
            size_t end_index = start_index + work_size;

            if (work_size > 0)
            {
                threads.emplace_back(
                    [this, &trajectories, n_steps, step_size, start_index, end_index]()
                    {
                        for (size_t j = start_index; j < end_index; j++)
                        {
                            OverdampedLangevin<Dim> odl(kbT, force, initial_positions[j]);
                            trajectories[j] = odl.getTrajectory(n_steps, step_size);
                        }
                    }
                );
            }
            start_index = end_index;
        }

        for (auto& thread : threads)
        {
            if (thread.joinable()) {
                thread.join();
            }
        }

        return trajectories;
    }

    py::list getBatchedTrajectoryPy(size_t n_steps, double step_size)
    {
        py::gil_scoped_release release;
        auto trajectories = getBatchedTrajectory(n_steps, step_size);
        py::gil_scoped_acquire acquire;
        
        py::list out_list;
        for (auto& arr : trajectories)
        {
            out_list.append(convert_from_shared_ptr<Dim>(arr));
        }
        return py::list(out_list);
    }
private:
    double kbT;
    force_callable force;
    std::vector<position<Dim>> initial_positions;
};
