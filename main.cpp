#include  "read_PLOT3D.h"
#include "cell.h"
#include "SpatialDiscretization.h"
#include "TemporalDiscretization.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>

std::vector<std::vector<double>> thomasAlgorithm1(const std::vector<double>& a, // subdiagonal
                                                 const std::vector<double>& b, // main diagonal
                                                 const std::vector<double>& c, // superdiagonal
                                                 const std::vector<std::vector<double>>& d) {  // right hand side
    int n = b.size();
    int numRHS = d.size();

    // Initialize modified vectors
    std::vector<double> cp(n, 0.0);                          // Modified super-diagonal
    std::vector<std::vector<double>> dp(numRHS, std::vector<double>(n, 0.0));  // Modified right-hand side
    std::vector<std::vector<double>> x(numRHS, std::vector<double>(n, 0.0));   // Solution vector

    // Forward elimination for each right-hand side
    cp[0] = c[0] / b[0];
    for (int j = 0; j < numRHS; j++) {
        dp[j][0] = d[j][0] / b[0];
    }

    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / m;
        for (int j = 0; j < numRHS; j++) {
            dp[j][i] = (d[j][i] - a[i] * dp[j][i - 1]) / m;
        }
    }

    // Back substitution for each right-hand side
    for (int j = 0; j < numRHS; j++) {
        x[j][n - 1] = dp[j][n - 1];
        for (int i = n - 2; i >= 0; i--) {
            x[j][i] = dp[j][i] - cp[i] * x[j][i + 1];
        }
    }

    return x;
}

int main() {

    // Read the PLOT3D mesh from a file
    auto [x, y] = read_PLOT3D_mesh("../mesh/x.9");

    // Output the dimensions and some values for verification
    std::cout << "Grid dimensions: " << x.size() << " x " << x[0].size() << std::endl;

    constexpr double Mach = 0.8;
    constexpr double alpha = 1.25*M_PI/180;
    constexpr double p_inf = 1E5;
    constexpr double T_inf = 215.0;
    constexpr double rho_inf = p_inf/(T_inf*287);

    constexpr double a = std::sqrt(1.4*287*T_inf);
    constexpr double Vitesse = Mach*a;
    constexpr double u_inf = Vitesse*std::cos(alpha);
    constexpr double v_inf = Vitesse*std::sin(alpha);
    constexpr double E_inf = p_inf/((1.4-1)*rho_inf) + 0.5*std::pow(Vitesse, 2);

    constexpr double l_ref = 1.0;
    constexpr double U_ref = std::sqrt(p_inf/rho_inf);

    constexpr double rho = 1.0;
    constexpr double u = u_inf/U_ref;
    constexpr double v = v_inf/U_ref;
    constexpr double E = E_inf/(U_ref*U_ref);
    constexpr double T = 1.0;
    constexpr double p = 1.0;

    // // SpatialDiscretization current_state(x, y, rho, u, v, E, T, p, T_inf, U_ref);
    // SpatialDiscretization current_state(x, y, rho_inf, u_inf, v_inf, E_inf, T_inf, p_inf, 1, 1);
    // current_state.run_even();

    TemporalDiscretization FVM(x, y, rho, u, v, E, T, p, T_inf, U_ref);
    auto[q, q_vertex, Residuals] = FVM.RungeKutta(50000);

    TemporalDiscretization::save_checkpoint(q, {static_cast<int>(Residuals.size())}, Residuals, "checkpoint_test.txt");
    write_plot3d_2d(x, y, q_vertex, Mach, alpha, 0, 0, rho_inf, U_ref,"test.xy", "test.q");
    std::cout << "PLOT3D files written successfully." << std::endl;


    return 0;

}














