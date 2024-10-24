#include "TemporalDiscretization.h"
#include "vector_helper.h"
#include "read_PLOT3D.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

TemporalDiscretization::TemporalDiscretization(const std::vector<std::vector<double>>& x,
                                               const std::vector<std::vector<double>>& y,
                                               const double& rho,
                                               const double& u,
                                               const double& v,
                                               const double& E,
                                               const double& T,
                                               const double& p,
                                               const double& T_ref,
                                               const double& U_ref)
    : x(x), y(y), rho(rho), u(u), v(v), E(E), T(T), p(p), T_ref(T_ref), U_ref(U_ref),
      current_state(x, y, rho, u, v, E, T, p, T_ref, U_ref) {}

double TemporalDiscretization::compute_dt(const cell& cell_IJ, const double sigma) const {
    // Extract conservative variables from the cell
    auto [rho_IJ, u_IJ, v_IJ, E_IJ, T_IJ, p_IJ] = current_state.SpatialDiscretization::conservative_variable_from_W(cell_IJ.W);
    double c_IJ = std::sqrt(1.4 * 287 * T_IJ * T_ref)/U_ref;  // Speed of sound

    // Calculate normal vectors and Ds
    const std::vector<double> n_I = vector_scale(0.5, vector_subtract(cell_IJ.n2, cell_IJ.n4));
    const std::vector<double> n_J = vector_scale(0.5, vector_subtract(cell_IJ.n1, cell_IJ.n3));
    double Ds_I = 0.5 * (cell_IJ.Ds2 + cell_IJ.Ds4);
    double Ds_J = 0.5 * (cell_IJ.Ds1 + cell_IJ.Ds3);

    // Calculate lambda_I and lambda_J
    double lambda_I = (std::abs(u_IJ * n_I[0] + v_IJ * n_I[1]) + c_IJ) * Ds_I;
    double lambda_J = (std::abs(u_IJ * n_J[0] + v_IJ * n_J[1]) + c_IJ) * Ds_J;

    // Compute time step
    const double dt = sigma * cell_IJ.OMEGA / (lambda_I + lambda_J);
    // std::cout << dt << std::endl;

    return dt;
}

std::vector<double> TemporalDiscretization::compute_L2_norm(const std::vector<std::vector<std::vector<double>>> &residuals) {
    const auto m = residuals.size();          // Number of rows
    const auto n = residuals[0].size();       // Number of columns
    const auto num_components = residuals[0][0].size(); // Number of components (4)

    std::vector<double> l2_norms(num_components, 0.0);
    const auto N = static_cast<double>(m * n); // Total number of cells

    for (int k = 0; k < num_components; ++k) {
        double sum = 0.0;
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                sum += residuals[j][i][k] * residuals[j][i][k]; // Squared component
            }
        }
        l2_norms[k] = std::sqrt(sum / N);
    }

    return l2_norms;
}

void TemporalDiscretization::save_checkpoint(const std::vector<std::vector<std::vector<double>>>& q,
                                             const std::vector<int>& iteration,
                                             const std::vector<std::vector<double>>& Residuals,
                                             const std::string& file_name) {
    // Open the file for writing
    std::ofstream file(file_name);

    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file for saving checkpoint." << std::endl;
        return;
    }

    // Write the current iteration numbers
    file << "Iterations:\n";
    for (const auto& it : iteration) {
        file << it << " ";
    }
    file << "\n";

    // Get dimensions for q
    const auto nj = q.size();               // Number of rows (j)
    const auto ni = q[0].size();            // Number of columns (i)
    const auto nVars = q[0][0].size();      // Number of variables (e.g., density, momentum, energy)

    // Write the solution vector q
    file << "Solution vector (q):\n";
    file << ni << " " << nj << " " << nVars << "\n";  // Write dimensions for q

    // Write flow variables (density, x-momentum, y-momentum, energy) in row-major order
    for (int n = 0; n < nVars; ++n) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                  // Iterate over the number of variables
                    file << std::scientific << std::setprecision(16) << q[j][i][n] << "\n";
                }
            }
    }

    // Write the residuals
    file << "Residuals (4D vectors):\n";
    for (const auto& res : Residuals) {
        for (const auto& val : res) {
            file << std::scientific << std::setprecision(16) << val << " ";  // Write each value in the residual vector
        }
        file << "\n";  // New line after each row of residuals
    }

    // Close the file
    file.close();
    std::cout << "Checkpoint saved to " << file_name << std::endl;
}


std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<int>,
           std::vector<std::vector<double>>> TemporalDiscretization::load_checkpoint(const std::string& file_name) {
    // Open the file for reading
    std::ifstream file(file_name);

    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file for loading checkpoint." << std::endl;
        return std::make_tuple(std::vector<std::vector<std::vector<double>>>(), std::vector<int>(), std::vector<std::vector<double>>());
    }

    std::string line;
    std::vector<int> iteration;
    std::vector<std::vector<std::vector<double>>> q;
    std::vector<std::vector<double>> Residuals;

    // Read iteration numbers
    std::getline(file, line);  // Read the header line
    std::getline(file, line);  // Read the iteration line
    std::istringstream iss(line);
    int it;
    while (iss >> it) {
        iteration.push_back(it);
    }

    // Read dimensions for q
    std::getline(file, line);  // Read the header for the solution vector
    std::getline(file, line);  // Read the dimensions line
    std::istringstream dimStream(line);
    int ni, nj, nVars;
    dimStream >> ni >> nj >> nVars;

    // Resize the solution vector q
    q.resize(nj, std::vector<std::vector<double>>(ni, std::vector<double>(nVars)));

    // Read the flow variables (density, x-momentum, y-momentum, energy) in row-major order
    for (int n = 0; n < nVars; ++n) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                double value;
                file >> value;  // Read the value directly
                q[j][i][n] = value;
            }
        }
    }

    // Read the residuals
    std::getline(file, line);  // Read the header for residuals
    std::getline(file, line);  // Read the residuals line
    while (std::getline(file, line)) {
        std::istringstream resStream(line);
        std::vector<double> resRow;
        double val;
        while (resStream >> val) {
            resRow.push_back(val);
        }
        Residuals.push_back(resRow);
    }

    // Close the file
    file.close();
    std::cout << "Checkpoint loaded from " << file_name << std::endl;

    return std::make_tuple(q, iteration, Residuals);
}

std::tuple<std::vector<std::vector<std::vector<double>>>,
           std::vector<std::vector<std::vector<double>>>,
           std::vector<std::vector<double>>> TemporalDiscretization::RungeKutta(int it_max) {

    double a1 = 0.25; double b1 = 1.0;
    double a2 = 0.1667; double b2 = 0.0;
    double a3 = 0.3750; double b3 = 0.56;
    double a4 = 0.5; double b4 = 0.0;
    double a5 = 1.0; double b5 = 0.44;

    auto ny = current_state.cells.size();
    auto nx = current_state.cells[0].size();
    std::cout << ny << " " << nx << std::endl;

    current_state.run_even();
    // Initialize R_d0
    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {
            current_state.cells[j][i].R_d0 = current_state.cells[j][i].R_d;
        }
    }

    std::vector<std::vector<double>> Residuals;
    std::vector<int> iteration;


    Residuals = std::vector<std::vector<double>>{};
    iteration = std::vector<int>{};


    std::vector<std::vector<std::vector<double>>> all_Res(ny - 4, std::vector<std::vector<double>>(nx, std::vector<double>(4, 1.0)));
    std::vector<std::vector<std::vector<double>>> all_dw(ny - 4, std::vector<std::vector<double>>(nx, std::vector<double>(4, 1.0)));
    std::vector<std::vector<std::vector<double>>> q(ny - 4, std::vector<std::vector<double>>(nx, std::vector<double>(4, 1.0)));

    // Fill q array
    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {
            q[j - 2][i][0] = current_state.cells[j][i].W[0];
            q[j - 2][i][1] = current_state.cells[j][i].W[1];
            q[j - 2][i][2] = current_state.cells[j][i].W[2];
            q[j - 2][i][3] = current_state.cells[j][i].W[3];
        }
    }

    try {
        std::vector<double> first_residual;
        int it = 0;
        std::vector<double> normalized_residuals = {1, 1, 1, 1};

        while (it < it_max) {
            // Stage 1
            for (int j = 2; j < ny - 2; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double dt = compute_dt(current_state.cells[j][i]);
                    const std::vector<double>& Rd0 = current_state.cells[j][i].R_d0;
                    std::vector<double> dW = vector_scale(-a1 * dt / current_state.cells[j][i].OMEGA, vector_subtract(current_state.cells[j][i].R_c, Rd0));
                    current_state.cells[j][i].W = vector_add(current_state.cells[j][i].W, dW) ;
                }
            }
            current_state.run_odd();

            // Stage 2
            for (int j = 2; j < ny - 2; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double dt = compute_dt(current_state.cells[j][i]);
                    const std::vector<double>& Rd0 = current_state.cells[j][i].R_d0;
                    std::vector<double> dW = vector_scale(-a2 * dt / current_state.cells[j][i].OMEGA, vector_subtract(current_state.cells[j][i].R_c, Rd0));
                    current_state.cells[j][i].W = vector_add(current_state.cells[j][i].W, dW) ;
                }
            }
            current_state.run_even();

            // Stage 3
            for (int j = 2; j < ny - 2; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double dt = compute_dt(current_state.cells[j][i]);
                    const std::vector<double> Rd20 = vector_add(vector_scale(b3, current_state.cells[j][i].R_d), vector_scale(1-b3, current_state.cells[j][i].R_d0));
                    current_state.cells[j][i].R_d0 = Rd20;
                    std::vector<double> dW = vector_scale(-a3 * dt / current_state.cells[j][i].OMEGA,vector_subtract(current_state.cells[j][i].R_c, Rd20));
                    current_state.cells[j][i].W = vector_add(current_state.cells[j][i].W, dW) ;
                }
            }
            current_state.run_odd();

            // Stage 4
            for (int j = 2; j < ny - 2; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double dt = compute_dt(current_state.cells[j][i]);
                    const std::vector<double>& Rd20 = current_state.cells[j][i].R_d0;
                    std::vector<double> dW = vector_scale(-a4 * dt / current_state.cells[j][i].OMEGA, vector_subtract(current_state.cells[j][i].R_c, Rd20));
                    current_state.cells[j][i].W = vector_add(current_state.cells[j][i].W, dW) ;
                }
            }
            current_state.run_even();

            // Stage 5, Final update
            for (int j = 2; j < ny - 2; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double dt = compute_dt(current_state.cells[j][i]);
                    const std::vector<double> Rd42 = vector_add(vector_scale(b5, current_state.cells[j][i].R_d), vector_scale(1-b5, current_state.cells[j][i].R_d0));
                    current_state.cells[j][i].R_d0 = Rd42;
                    std::vector<double> Res = vector_subtract(current_state.cells[j][i].R_c, Rd42);
                    std::vector<double> dW = vector_scale(-a5 * dt / current_state.cells[j][i].OMEGA, Res);
                    current_state.cells[j][i].W = vector_add(current_state.cells[j][i].W, dW) ;

                    all_Res[j - 2][i] = Res;
                    all_dw[j - 2][i] = dW;
                    q[j - 2][i][0] = current_state.cells[j][i].W[0];
                    q[j - 2][i][1] = current_state.cells[j][i].W[1];
                    q[j - 2][i][2] = current_state.cells[j][i].W[2];
                    q[j - 2][i][3] = current_state.cells[j][i].W[3];
                }
            }
            current_state.run_odd();


            // Compute L2 norm (placeholder logic)
            std::vector<double> l2_norm = compute_L2_norm(all_dw);

            if (it == 0) {
              first_residual = l2_norm;
            }
            normalized_residuals = vector_divide(l2_norm, first_residual);
            iteration.push_back(it);
            Residuals.push_back(l2_norm);

            std::cout << "Iteration " << it << ": L2 Norms = ";
            for (const auto &norm : normalized_residuals) {
                std::cout << norm << " ";
            }
            std::cout << std::endl;

            // Check for convergence
            if (*std::ranges::max_element(normalized_residuals) <= 1e-11) {
                break; // Exit the loop if convergence criterion is met
            }

            // Save checkpoint at each 1000 iteration
            if (it%1000 == 0) {
                save_checkpoint(q, iteration, Residuals);
            }



            it++;

        }



        std::vector<std::vector<std::vector<double>>> q_cell_dummy(ny - 2, std::vector<std::vector<double>>(nx, std::vector<double>(4, 1.0)));

        // Compute q_vertex
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 0; i < nx; ++i) {
                q_cell_dummy[j - 1][i][0] = current_state.cells[j][i].W[0];
                q_cell_dummy[j - 1][i][1] = current_state.cells[j][i].W[1];
                q_cell_dummy[j - 1][i][2] = current_state.cells[j][i].W[2];
                q_cell_dummy[j - 1][i][3] = current_state.cells[j][i].W[3];
            }
        }
        std::vector<std::vector<std::vector<double>>> q_vertex = cell_dummy_to_vertex_centered_airfoil(q_cell_dummy);

        return std::make_tuple(q, q_vertex, Residuals);
    }

    catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }

    std::vector<std::vector<std::vector<double>>> q_vertex = cell_dummy_to_vertex_centered_airfoil(q);
    return std::make_tuple(q, q_vertex, Residuals);
}

void TemporalDiscretization::run() {
    auto[q, q_vertex, Residuals] = TemporalDiscretization::RungeKutta();
}
