#include "IPS.h"
#include "Integrator.h"
#include <chrono>
#include <random>
#include "Potential.h"

int main() {
    return 0;
}
// int main() {
//     constexpr size_t dim = 2;

//     auto pair_force = Spring(1.0, 1.0);
//     auto confinement_force = RadialConfinement<dim>(2.0);

//     auto particles = generate_random_init<DataLayout::SoA, dim>(10, -1, 1, -2, 2);
//     auto integrator = ABOBA();
//     IPS_Simulator<decltype(particles), decltype(integrator)> ips(particles);
//     ips.pair_force = pair_force;
//     ips.confinement_force = confinement_force;
//     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//     for (size_t i = 0; i < 100; i++)
//     {
//         ips.integrate(0.01);
//     }

//     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//     auto one_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
//     std::cout << "# Time difference = " << one_time << "[ms]" << std::endl;

//     return 0;

// }