/**
 * @brief use Euler-Maruyama method to solve the overdamped Langevin equation, and show some results regarding the step-size.
 * 
 */

#include "OverdampedLangevin.h"
#include <cstddef>
#include <iostream>
#include <vector>
#include <chrono>



auto generate_random_init(size_t num) -> std::vector<position<1>>{
    std::vector<position<1>> inits;
    std::cout << "Generating " << num << " random initial positions." << std::endl;
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
 * @param argv 
 * @return int 
 */
int main(int argc, char *argv[])
{

    // double well force term:  F(x) = -dU/dx = -4x(x^2 - 1)
    auto force_func = [](position<1>& pos){return -4 * pos.x[0] * (pos.x[0] * pos.x[0] - 1);};

    auto init_pos = generate_random_init(1);
    OverdampedLangevin<1> odLangevin(1, force_func, init_pos[0]);

    // measure the time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    auto traj = odLangevin.getTrajectory(1000000, 0.01);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    auto inits = generate_random_init(512);
    std::cout << inits.size() << std::endl;
    auto b_odl = BatchedOverdampedLangevin<1>(1, force_func, inits);
    // record the time
    begin = std::chrono::steady_clock::now();
    auto batched_traj = b_odl.getBatchedTrajectory(1000000, 0.01);
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}