#include  "read_PLOT3D.h"
#include "cell.h"
#include "SpatialDiscretization.h"
#include "TemporalDiscretization.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <omp.h>
#include <chrono>

// Function to calculate sum of squares serially
long long serialSumOfSquares(const std::vector<int>& data) {
    long long sum = 0;
    for (int i = 0; i < data.size(); ++i) {
        sum += data[i] * data[i];
    }
    return sum;
}

// Function to calculate sum of squares using OpenMP parallelization
long long parallelSumOfSquares(const std::vector<int>& data, int numThreads) {
    long long sum = 0;

    // Set the number of threads to be used
    omp_set_num_threads(numThreads);

    // OpenMP parallelization with reduction for sum
    #pragma omp parallel
    {
        // Print the thread number for each iteration (can be useful to see thread usage)
        int threadID = omp_get_thread_num();

        // This block will only be executed by the master thread
        #pragma omp master
        {
            std::cout << "Number of threads used: " << omp_get_num_threads() << "\n";
        }

        // Parallelized loop with reduction
        #pragma omp for reduction(+:sum)
        for (int i = 0; i < data.size(); ++i) {
            sum += data[i] * data[i];
        }
    }
    return sum;
}



int main() {

    // // Read the PLOT3D mesh from a file
    // auto [x, y] = read_PLOT3D_mesh("../mesh/x.8");
    //
    // // Output the dimensions and some values for verification
    // std::cout << "Grid dimensions: " << x.size() << " x " << x[0].size() << std::endl;
    //
    // constexpr double Mach = 0.8;
    // constexpr double alpha = 1.25*M_PI/180;
    // constexpr double p_inf = 1E5;
    // constexpr double T_inf = 215.0;
    // constexpr double rho_inf = p_inf/(T_inf*287);
    //
    // constexpr double a = std::sqrt(1.4*287*T_inf);
    // constexpr double Vitesse = Mach*a;
    // constexpr double u_inf = Vitesse*std::cos(alpha);
    // constexpr double v_inf = Vitesse*std::sin(alpha);
    // constexpr double E_inf = p_inf/((1.4-1)*rho_inf) + 0.5*std::pow(Vitesse, 2);
    //
    // constexpr double l_ref = 1.0;
    // constexpr double U_ref = std::sqrt(p_inf/rho_inf);
    //
    // constexpr double rho = 1.0;
    // constexpr double u = u_inf/U_ref;
    // constexpr double v = v_inf/U_ref;
    // constexpr double E = E_inf/(U_ref*U_ref);
    // constexpr double T = 1.0;
    // constexpr double p = 1.0;
    //
    // // // SpatialDiscretization current_state(x, y, rho, u, v, E, T, p, T_inf, U_ref);
    // // SpatialDiscretization current_state(x, y, rho_inf, u_inf, v_inf, E_inf, T_inf, p_inf, 1, 1);
    // // current_state.run_even();
    //
    // TemporalDiscretization FVM(x, y, rho, u, v, E, T, p, T_inf, U_ref);
    // auto[q, q_vertex, Residuals] = FVM.RungeKutta(50000);
    //
    // TemporalDiscretization::save_checkpoint(q, {static_cast<int>(Residuals.size())}, Residuals, "checkpoint_test.txt");
    // write_plot3d_2d(x, y, q_vertex, Mach, alpha, 0, 0, rho_inf, U_ref,"test.xy", "test.q");
    // std::cout << "PLOT3D files written successfully." << std::endl;
    //
    //
    // return 0;

}














