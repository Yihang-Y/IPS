/**
 * @brief use Euler-Maruyama method to solve the overdamped Langevin equation, and show some results regarding the step-size.
 * 
 */

#include "OverdampedLangevin.h"
#include <cstddef>
#include <iostream>
#include <vector>
#include <chrono>

#define NUM_TASKS 16


auto generate_random_init(size_t num) -> std::vector<vec<1>>{
    std::vector<vec<1>> inits;
    // std::cout << "# Generating " << num << " random initial positions." << std::endl;
    inits.reserve(num);
    for (size_t i = 0; i < num; i++)
    {
        // just generate a random number between -1 and 1
        inits.push_back({{2 * (double)rand() / RAND_MAX - 1}});
    }
    return inits;
}

/**
 * @brief 
 * 
 * @param argc
 * @param argv: <number of steps> <step size> <number of tasks>
 */
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <number of steps> <step size> <number of tasks>" << std::endl;
        return 1;
    }

    size_t n_steps = std::stoul(argv[1]);
    double step_size = std::stod(argv[2]);
    size_t num_tasks = std::stoul(argv[3]);

    // double well force term:  F(x) = -dU/dx = -4x(x^2 - 1)
    auto force_func = [](vec<1>& pos) {
        return -4 * pos.x[0] * (pos.x[0] * pos.x[0] - 1);
    };

    if (num_tasks == 1)
    {
        auto init_pos = generate_random_init(1);
        OverdampedLangevin<1> odLangevin(1, force_func, init_pos[0]);

        // measure the time
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto traj = odLangevin.getTrajectory(n_steps, step_size);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto one_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        // std::cout << "# Time difference = " << one_time << "[ms]" << std::endl;

        // add a new line to separate the time and the trajectories
        // std::cout << std::endl;

        // output the trajectory
        // std::cout << "# Trajectory 0:" << std::endl;
        std::cout << "# NUM_STEPS: " << n_steps << std::endl;
        std::cout << "# STEP_SIZE: " << step_size << std::endl;
        for (size_t i = 0; i < n_steps; i++)
        {
            std::cout << traj->at(i)[0] << " ";
        }
    }
    else
    {
        auto inits = generate_random_init(num_tasks);
        auto b_odl = BatchedOverdampedLangevin<1>(1, force_func, inits);
        // record the time
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto batched_traj = b_odl.getBatchedTrajectory(n_steps, step_size);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // std::cout << "# Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        
        // add a new line to separate the time and the trajectories
        // std::cout << std::endl;

        // output the trajectories
        for (size_t i = 0; i < num_tasks; i++)
        {
            // std::cout << "# Trajectory " << i << ":" << std::endl;
            std::cout << "# NUM_STEPS: " << n_steps << std::endl;
            std::cout << "# STEP_SIZE: " << step_size << std::endl;
            for (size_t j = 0; j < n_steps; j++)
            {
                std::cout << batched_traj[i]->at(j)[0] << " ";
            }
            std::cout << std::endl;
        }
    }
    return 0;
}