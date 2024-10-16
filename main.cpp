#include  "read_PLOT3D.h"
#include "cell.h"
#include "SpatialDiscretization.h"
#include "TemporalDiscretization.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>

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

int main() {
    // Read the PLOT3D mesh from a file
    auto [x, y] = read_PLOT3D_mesh("../mesh/x.9");

    // Output the dimensions and some values for verification
    std::cout << "Grid dimensions: " << x.size() << " x " << x[0].size() << std::endl;

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
    auto[q, q_vertex, Residuals] = FVM.RungeKutta(50000);

    TemporalDiscretization::save_checkpoint(q, {static_cast<int>(Residuals.size())}, Residuals, "checkpoint_test.txt");
    write_plot3d_2d(x, y, q_vertex, Mach, alpha, 0, 0, "test.xy", "test.q");
    std::cout << "PLOT3D files written successfully." << std::endl;


    return 0;

}
