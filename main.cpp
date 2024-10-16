#include  "read_PLOT3D.h"
#include "cell.h"
#include "SpatialDiscretization.h"
#include "TemporalDiscretization.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <csignal>

std::tuple<double, double, double, double, double, double> conservativeVariableFromW(const std::vector<double>& W, const double gamma=1.4) {
    double variable[4]; // Assuming 4 variables in W
    for (int i = 0; i < 4; ++i) {
        variable[i] = W[i] / W[0]; // Normalize W
    }

    double rho = W[0];
    double u = variable[1];
    double v = variable[2];
    double E = variable[3];
    double p = (gamma - 1) * rho * (E - (u * u + v * v) / 2);
    double T = p / (rho * 287); // Assume gas constant R = 287

    return std::make_tuple(rho, u, v, E, T, p);
}

void signal_handler_1(int signal_num)
{
    std::cout << "The interrupt signal is (" << signal_num
         << "). \n";

    // It terminates the  program
    exit(signal_num);
}

int main() {
    // Read the PLOT3D mesh from a file
    auto [x, y] = read_PLOT3D_mesh("../mesh/x.6");

    // Output the dimensions and some values for verification
    std::cout << "Grid dimensions: " << x.size() << " x " << x[0].size() << std::endl;
    std::cout << "Sample value at (0,0): " << x[0][0] << ", " << y[0][0] << std::endl;

    // auto [ni, nj, mach, alpha, reyn, time, q] = read_PLOT3D_solution("solution_file/test_M05_alpha0.q");
    //
    // // Output read values for verification
    // std::cout << "Grid dimensions: " << ni << " x " << nj << std::endl;
    // std::cout << "Mach: " << mach << ", Alpha: " << alpha << ", Reynolds: " << reyn << ", Time: " << time << std::endl;
    //
    // // Example of accessing flow variable at (j, i, n)
    // std::cout << "Flow variable (nj=0, ni=0, n=0): " << q[0][0][0] << std::endl;
    //
    // write_plot3d_2d(x, y, q, mach, alpha, reyn, time, "test.xy", "test.q");
    // std::cout << "PLOT3D files written successfully." << std::endl;

    constexpr double Mach = 0.8;
    constexpr double alpha = 1.25*M_PI/180;
    constexpr double p_inf = 1E5;
    constexpr double T_inf = 288;

    constexpr double a = std::sqrt(1.4*287*T_inf);
    constexpr double Vitesse = Mach*a;
    constexpr double u = Vitesse*std::cos(alpha);
    constexpr double v = Vitesse*std::sin(alpha);
    constexpr double rho = p_inf/(T_inf*287);
    constexpr double E = p_inf/((1.4-1)*rho) + 0.5*std::pow(Vitesse, 2);

    // SpatialDiscretization current_state(x, y, rho, u, v, E, T_inf, p_inf);
    // current_state.run();
    // std::vector<std::vector<cell>> cells = current_state.cells;
    // cell cell_test = cells[2][2];

    TemporalDiscretization FVM(x, y, rho, u, v, E, T_inf, p_inf);
    auto[q, q_vertex, Residuals] = FVM.RungeKutta(1000);

    TemporalDiscretization::save_checkpoint(q, {static_cast<int>(Residuals.size())}, Residuals, "checkpoint_test_3.txt");
    write_plot3d_2d(x, y, q_vertex, Mach, alpha, 0, 0, "test_M08_alpha125_x6.xy", "test_M08_alpha125_x6.q");
    std::cout << "PLOT3D files written successfully." << std::endl;


    return 0;

}
