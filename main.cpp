/**
 * @brief use Euler-Maruyama method to solve the overdamped Langevin equation, and show some results regarding the step-size.
 * 
 */

#include "OverdampedLangevin.h"
#include <cstddef>
#include <iostream>

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
    OverdampedLangevin<1> odLangevin(1, force_func, position<1>{{1}});

    auto traj = odLangevin.getTrajectory(1000, 0.01);
    for (size_t i = 0; i < traj.size(); i++)
    {
        std::cout << traj[i][0] << std::endl;
    }

    return 0;
}