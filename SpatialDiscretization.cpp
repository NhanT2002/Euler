#include "SpatialDiscretization.h"
#include "vector_helper.h"

#include <array>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>



SpatialDiscretization::SpatialDiscretization(const std::vector<std::vector<double>>& x,
                          const std::vector<std::vector<double>>& y,
                          const double& rho,
                          const double& u,
                          const double& v,
                          const double& E,
                          const double& T,
                          const double& p)
                              :x(x), y(y), rho(rho), u(u), v(v), E(E), T(T), p(p) {
    ny = y.size();
    nx = x[0].size();
    // Cells generation
    domain_cells.resize(ny - 1, std::vector<cell>(nx - 1));

    for (size_t j = 0; j < ny - 1; ++j) {
        for (size_t i = 0; i < nx - 1; ++i) {
            domain_cells[j][i] = cell(x[j][i], y[j][i],
                                       x[j][i + 1], y[j][i + 1],
                                       x[j + 1][i + 1], y[j + 1][i + 1],
                                       x[j + 1][i], y[j + 1][i],
                                       rho, u, v, E);
        }
    }

    // Initialize solid wall cells (dummy cells)
    solid_wall_cells.resize(2, std::vector<cell>(nx - 1));

    // Initialize farfield cells (dummy cells)
    farfield_cells.resize(2, std::vector<cell>(nx - 1));

    // Combine all cells
    cells.resize(solid_wall_cells.size() + domain_cells.size() + farfield_cells.size());
    std::copy(solid_wall_cells.begin(), solid_wall_cells.end(), cells.begin());
    std::copy(domain_cells.begin(), domain_cells.end(), cells.begin() + solid_wall_cells.size());
    std::copy(farfield_cells.begin(), farfield_cells.end(), cells.begin() + solid_wall_cells.size() + domain_cells.size());
}

void SpatialDiscretization::compute_dummy_cells() {
    // Solid wall
    for (size_t i = 0; i < nx - 1; ++i) {
        double p3, p4;
        auto [rho_val, u_val, v_val, E_val, T_val, p2] = conservative_variable_from_W(cells[2][i].W);
        std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p3) = conservative_variable_from_W(cells[3][i].W);
        std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p4) = conservative_variable_from_W(cells[4][i].W);

        double pw = (15 * p2 - 10 * p3 + 3 * p4) / 8.0; // Blazek
        double p1 = 2 * pw - p2;
        std::vector<double> vel = {u_val, v_val};

        std::vector<double> n = cells[2][i].n1;

        std::vector<std::vector<double>> R = { {-n[1], n[0]}, {n[0], n[1]} };
        double q_t = -R[0][0] * vel[0] - R[0][1] * vel[1];
        double q_n = -R[1][0] * vel[0] - R[1][1] * vel[1];

        double y_eta = cells[2][i].s1[0] / cells[2][i].Ds1;
        double x_eta = cells[2][i].s1[1] / cells[2][i].Ds1;

        // Swanson Turkel
        double u_dummy = x_eta * q_t + y_eta * q_n;
        double v_dummy = -y_eta * q_t + x_eta * q_n;


        E_val = p1 / (1.4 - 1) / rho_val + 0.5 * (u_dummy * u_dummy + v_dummy * v_dummy);

        cells[0][i] = cell(x[0][i], y[0][i], x[0][i + 1], y[0][i + 1],
                                       x[1][i + 1], y[1][i + 1],
                                       x[1][i], y[1][i],
                                       rho_val, u_dummy, v_dummy, E_val);
        cells[1][i] = cell(x[0][i], y[0][i], x[0][i + 1], y[0][i + 1],
                                       x[1][i + 1], y[1][i + 1],
                                       x[1][i], y[1][i],
                                       rho_val, u_dummy, v_dummy, E_val);
    }

    // Farfield
    for (size_t i = 0; i < nx - 1; ++i) {
        auto [rho_val, u_val, v_val, E_val, T_val, p_val] = conservative_variable_from_W(cells[cells.size() - 3][i].W);
        double c = std::sqrt(1.4 * 287 * T_val);
        double M = std::sqrt(u_val * u_val + v_val * v_val) / c;
        std::vector<double> n = cells[cells.size() - 3][i].n3;

        if (u_val * n[0] + v_val * n[1] > 0) { // Out of cell
            if (M >= 1) {
                cells[cells.size()-1][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                             rho_val, u_val, v_val, E_val);
                cells[cells.size()-2][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                             rho_val, u_val, v_val, E_val);
            }
            else {  // Subsonic
                double p_b = this->p;  // Boundary pressure
                double rho_b = rho_val + (p_b - p_val) / (c * c);
                double u_b = u_val + n[0] * (p_val - p_b) / (rho_val * c);
                double v_b = v_val + n[1] * (p_val - p_b) / (rho_val * c);
                double E_b = p_b / ((1.4 - 1) * rho_b) + 0.5 * (u_b * u_b + v_b * v_b);

                std::vector<double> W_b = {rho_b, rho_b * u_b, rho_b * v_b, rho_b * E_b};
                std::vector<double> W_a = {2 * W_b[0] - cells[cells.size() - 3][i].W[0],
                                       2 * W_b[1] - cells[cells.size() - 3][i].W[1],
                                       2 * W_b[2] - cells[cells.size() - 3][i].W[2],
                                       2 * W_b[3] - cells[cells.size() - 3][i].W[3]};

                auto [rho_a, u_a, v_a, E_a, T_a, p_a] = conservative_variable_from_W(W_a);


                cells[cells.size()-1][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            rho_a, u_a, v_a, E_a);
                cells[cells.size()-2][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            rho_a, u_a, v_a, E_a);
            }
        }
        else {  // Moving into the cell
            if (M >= 1) {  // Supersonic
                cells[cells.size()-1][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            this->rho, this->u, this->v, this->E);
                cells[cells.size()-2][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            this->rho, this->u, this->v, this->E);
            } else {  // Subsonic
                double p_b = 0.5 * (this->p + p_val - rho_val * c * (n[0] * (this->u - u_val) + n[1] * (this->v - v_val)));
                double rho_b = this->rho + (p_b - this->p) / (c * c);
                double u_b = this->u - n[0] * (this->p - p_b) / (rho_val * c);
                double v_b = this->v - n[1] * (this->p - p_b) / (rho_val * c);
                double E_b = p_b / ((1.4 - 1) * rho_b) + 0.5 * (u_b * u_b + v_b * v_b);

                std::vector<double> W_b = {rho_b, rho_b * u_b, rho_b * v_b, rho_b * E_b};
                std::vector<double> W_a = {2 * W_b[0] - cells[cells.size() - 3][i].W[0],
                                       2 * W_b[1] - cells[cells.size() - 3][i].W[1],
                                       2 * W_b[2] - cells[cells.size() - 3][i].W[2],
                                       2 * W_b[3] - cells[cells.size() - 3][i].W[3]};

                auto [rho_a, u_a, v_a, E_a, T_a, p_a] = conservative_variable_from_W(W_a);

                cells[cells.size()-1][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            rho_a, u_a, v_a, E_a);
                cells[cells.size()-2][i] = cell(x[ny - 2][i], y[ny - 2][i],
                                             x[ny - 2][i + 1], y[ny - 2][i + 1],
                                             x[ny - 1][i + 1], y[ny - 1][i + 1],
                                             x[ny - 1][i], y[ny - 1][i],
                                            rho_a, u_a, v_a, E_a);
            }
        }
    }
}

// Define the conservative_variable_from_W function as per your requirements
std::tuple<double, double, double, double, double, double> SpatialDiscretization::conservative_variable_from_W(const std::vector<double>& W) {
    // Implement the conversion from W to (rho, u, v, E)
    double rho = W[0];
    double u = W[1] / rho;
    double v = W[2] / rho;
    double E = W[3] / rho;
    double p = (1.4-1)*rho*(E-(u*u+v*v)/2);
    double T = p/(rho*287);
    return std::make_tuple(rho, u, v, E, T, p);
}

std::vector<double> SpatialDiscretization::FcDs(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) {
    auto [rho, u, v, E, T, p] = conservative_variable_from_W(W);
    double V = n[0]*u + n[1]*v;
    double H = E + p/rho;

    return {rho*V*Ds, (rho*u*V + n[0]*p)*Ds, (rho*v*V + n[1]*p)*Ds, rho*H*V*Ds};
}

double SpatialDiscretization::Lambdac(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) {
    auto [rho, u, v, E, T, p] = conservative_variable_from_W(W);
    double c = std::sqrt(1.4*287*T);
    double V = n[0]*u + n[1]*v;
    double lambda = (std::abs(V) + c)*Ds;

    return lambda;
}

void SpatialDiscretization::compute_Fc_DeltaS_Lambdac() {
    int ny = cells.size();
    int nx = cells[0].size();

    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {
            std::vector<double> avg_W1 = vector_scale(0.5, vector_add(cells[j][i].W, cells[j - 1][i].W));
            std::vector<double> avg_W2 = vector_scale(0.5, vector_add(cells[j][i].W, cells[j][(i + 1) % nx].W));
            std::vector<double> avg_W3 = vector_scale(0.5, vector_add(cells[j][i].W, cells[j + 1][i].W));
            std::vector<double> avg_W4 = vector_scale(0.5, vector_add(cells[j][i].W, cells[j][(i - 1 + nx) % nx].W));

            std::vector<double> FcDs_1 = FcDs(avg_W1, cells[j][i].n1, cells[j][i].Ds1);
            std::vector<double> FcDs_2 = FcDs(avg_W2, cells[j][i].n2, cells[j][i].Ds2);
            std::vector<double> FcDs_3 = FcDs(avg_W3, cells[j][i].n3, cells[j][i].Ds3);
            std::vector<double> FcDs_4 = FcDs(avg_W4, cells[j][i].n4, cells[j][i].Ds4);

            cells[j][i].FcDS_1 = FcDs_1;
            cells[j][i].FcDS_2 = FcDs_2;
            cells[j][i].FcDS_3 = FcDs_3;
            cells[j][i].FcDS_4 = FcDs_4;
        }
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // Calculate Lambda values
            std::vector<double> n1_n3 = vector_scale(0.5, vector_subtract(cells[j][i].n1, cells[j][i].n3));
            std::vector<double> n2_n4 = vector_scale(0.5, vector_subtract(cells[j][i].n2, cells[j][i].n4));
            std::vector<double> n3_n1 = vector_scale(0.5, vector_subtract(cells[j][i].n3, cells[j][i].n1));
            std::vector<double> n4_n2 = vector_scale(0.5, vector_subtract(cells[j][i].n4, cells[j][i].n2));

            double ds2_plus_ds4 = 0.5 * (cells[j][i].Ds2 + cells[j][i].Ds4);
            double ds1_plus_ds3 = 0.5 * (cells[j][i].Ds1 + cells[j][i].Ds3);

            // Compute Lambda values
            cells[j][i].Lambda_1_I = Lambdac(cells[j][i].W, n1_n3, ds2_plus_ds4);
            cells[j][i].Lambda_1_J = Lambdac(cells[j][i].W, n1_n3, ds1_plus_ds3);
            cells[j][i].Lambda_2_I = Lambdac(cells[j][i].W, n2_n4, ds2_plus_ds4);
            cells[j][i].Lambda_2_J = Lambdac(cells[j][i].W, n2_n4, ds1_plus_ds3);
            cells[j][i].Lambda_3_I = Lambdac(cells[j][i].W, n3_n1, ds2_plus_ds4);
            cells[j][i].Lambda_3_J = Lambdac(cells[j][i].W, n3_n1, ds1_plus_ds3);
            cells[j][i].Lambda_4_I = Lambdac(cells[j][i].W, n4_n2, ds2_plus_ds4);
            cells[j][i].Lambda_4_J = Lambdac(cells[j][i].W, n4_n2, ds1_plus_ds3);
        }
    }

}

std::tuple<double, double> SpatialDiscretization::compute_epsilon(const cell& cell_Im1, const cell& cell_I,
                                                                 const cell& cell_Ip1, const cell& cell_Ip2,
                                                                 double k2, double k4) {
    // Retrieve pressure from the conservative variables (assuming the last element is pressure)
    double p_Im1, p_I, p_Ip1, p_Ip2;
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Im1) = conservative_variable_from_W(cell_Im1.W);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_I) = conservative_variable_from_W(cell_I.W);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Ip1) = conservative_variable_from_W(cell_Ip1.W);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Ip2) = conservative_variable_from_W(cell_Ip2.W);

    // Calculate Gamma_I and Gamma_Ip1
    double Gamma_I = std::abs(p_Ip1 - 2.0 * p_I + p_Im1) / (p_Ip1 + 2.0 * p_I + p_Im1);
    double Gamma_Ip1 = std::abs(p_Ip2 - 2.0 * p_Ip1 + p_I) / (p_Ip2 + 2.0 * p_Ip1 + p_I);

    // Compute eps2 and eps4
    double eps2 = k2 * std::max(Gamma_I, Gamma_Ip1);
    double eps4 = std::max(0.0, k4 - eps2);

    // Return the results as a pair
    return std::make_tuple(eps2, eps4);
}

void SpatialDiscretization::compute_dissipation() {
    int ny = cells.size();
    int nx = cells[0].size();

    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {
            cell& cell_IJ = cells[j][i];
            cell& cell_Ip1J = cells[j][(i + 1) % nx];
            cell& cell_IJp1 = cells[j + 1][i];
            cell& cell_Im1J = cells[j][(i - 1 + nx) % nx];
            cell& cell_IJm1 = cells[j - 1][i];
            cell& cell_Ip2J = cells[j][(i + 2) % nx];
            cell& cell_IJp2 = cells[j + 2][i];
            cell& cell_Im2J = cells[j][(i - 2 + nx) % nx];
            cell& cell_IJm2 = cells[j - 2][i];

            // Lambda calculations
            double Lambda_2_I = 0.5 * (cell_IJ.Lambda_2_I + cell_Ip1J.Lambda_2_I);
            double Lambda_2_J = 0.5 * (cell_IJ.Lambda_2_J + cell_Ip1J.Lambda_2_J);
            double Lambda_2_S = Lambda_2_I + Lambda_2_J;

            double Lambda_3_J = 0.5 * (cell_IJ.Lambda_3_J + cell_IJp1.Lambda_3_J);
            double Lambda_3_I = 0.5 * (cell_IJ.Lambda_3_I + cell_IJp1.Lambda_3_I);
            double Lambda_3_S = Lambda_3_I + Lambda_3_J;

            double Lambda_4_I = 0.5 * (cell_IJ.Lambda_4_I + cell_Im1J.Lambda_4_I);
            double Lambda_4_J = 0.5 * (cell_IJ.Lambda_4_J + cell_Im1J.Lambda_4_J);
            double Lambda_4_S = Lambda_4_I + Lambda_4_J;

            double Lambda_1_J = 0.5 * (cell_IJ.Lambda_1_J + cell_IJm1.Lambda_1_J);
            double Lambda_1_I = 0.5 * (cell_IJ.Lambda_1_I + cell_IJm1.Lambda_1_I);
            double Lambda_1_S = Lambda_1_I + Lambda_1_J;

            cell_IJ.Lambda_1_S = Lambda_1_S;
            cell_IJ.Lambda_2_S = Lambda_2_S;
            cell_IJ.Lambda_3_S = Lambda_3_S;
            cell_IJ.Lambda_4_S = Lambda_4_S;

            // Epsilon calculations
            auto[eps2_2, eps4_2] = compute_epsilon(cell_Im1J, cell_IJ, cell_Ip1J, cell_Ip2J);
            auto[eps2_3, eps4_3] = compute_epsilon(cell_IJm1, cell_IJ, cell_IJp1, cell_IJp2);
            auto[eps2_4, eps4_4] = compute_epsilon(cell_Ip1J, cell_IJ, cell_Im1J, cell_Im2J);
            auto[eps2_1, eps4_1] = compute_epsilon(cell_IJp1, cell_IJ, cell_IJm1, cell_IJm2);

            cell_IJ.eps2_2 = eps2_2;
            cell_IJ.eps4_2 = eps4_2;
            cell_IJ.eps2_3 = eps2_3;
            cell_IJ.eps4_3 = eps4_3;
            cell_IJ.eps2_4 = eps2_4;
            cell_IJ.eps4_4 = eps4_4;
            cell_IJ.eps2_1 = eps2_1;
            cell_IJ.eps4_1 = eps4_1;

            // Dissipation terms
            std::vector<double> D_2 = vector_scale(Lambda_2_S,
                vector_subtract(
                    vector_scale(eps2_2, vector_subtract(cell_Ip1J.W, cell_IJ.W)),
                    vector_scale(eps4_2, vector_add(
                        vector_subtract(cell_Ip2J.W, vector_scale(3, cell_Ip1J.W)),
                        vector_subtract(
                            vector_scale(3, cell_IJ.W),
                            cell_Im1J.W
                        )
                    ))
                )
            );

            std::vector<double> D_3 = vector_scale(Lambda_3_S,
                vector_subtract(
                    vector_scale(eps2_3, vector_subtract(cell_IJp1.W, cell_IJ.W)),
                    vector_scale(eps4_3, vector_add(
                        vector_subtract(cell_IJp2.W, vector_scale(3, cell_IJp1.W)),
                        vector_subtract(
                            vector_scale(3, cell_IJ.W),
                            cell_IJm1.W
                        )
                    ))
                )
            );

            std::vector<double> D_4 = vector_scale(Lambda_4_S,
                vector_subtract(
                    vector_scale(eps2_4, vector_subtract(cell_Im1J.W, cell_IJ.W)),
                    vector_scale(eps4_4, vector_add(
                        vector_subtract(cell_Im2J.W, vector_scale(3, cell_Im1J.W)),
                        vector_subtract(
                            vector_scale(3, cell_IJ.W),
                            cell_Ip1J.W
                        )
                    ))
                )
            );

            std::vector<double> D_1 = vector_scale(Lambda_1_S,
                vector_subtract(
                    vector_scale(eps2_1, vector_subtract(cell_IJm1.W, cell_IJ.W)),
                    vector_scale(eps4_1, vector_add(
                        vector_subtract(cell_IJm2.W, vector_scale(3, cell_IJm1.W)),
                        vector_subtract(
                            vector_scale(3, cell_IJ.W),
                            cell_IJp1.W
                        )
                    ))
                )
            );

            cell_IJ.D_1 = D_1;
            cell_IJ.D_2 = D_2;
            cell_IJ.D_3 = D_3;
            cell_IJ.D_4 = D_4;
        }
    }

    // Boundary conditions
    for (int i = 0; i < nx; ++i) {

        // Calculate D_1 for cell (3, i)
        cells[3][i].D_1 = vector_scale(cells[3][i].Lambda_1_S,
            vector_subtract(
                vector_scale(cells[3][i].eps2_1, vector_subtract(cells[2][i].W, cells[3][i].W)),

                vector_scale(cells[4][i].eps4_1, vector_subtract(vector_subtract(
                    vector_scale(2.0, cells[3][i].W), cells[2][i].W), cells[4][i].W))
            )
        );

        // Calculate D_1 for cell (2, i)
        cells[2][i].D_1 = vector_scale(cells[2][i].Lambda_1_S,
            vector_subtract(
                vector_scale(cells[3][i].eps2_1, vector_subtract(cells[2][i].W, cells[3][i].W)),

                vector_scale(cells[4][i].eps4_1, vector_subtract(vector_subtract(
                    vector_scale(2.0, cells[3][i].W), cells[2][i].W), cells[4][i].W))
            )
        );
    }
}

void SpatialDiscretization::compute_R() {
    int ny = cells.size();
    int nx = cells[0].size();

    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {
            cell& current_cell = cells[j][i];

            // Extract the flux vectors
            const std::vector<double>& FcDS_1 = current_cell.FcDS_1;
            const std::vector<double>& FcDS_2 = current_cell.FcDS_2;
            const std::vector<double>& FcDS_3 = current_cell.FcDS_3;
            const std::vector<double>& FcDS_4 = current_cell.FcDS_4;

            // Extract the dissipation vectors
            const std::vector<double>& D_1 = current_cell.D_1;
            const std::vector<double>& D_2 = current_cell.D_2;
            const std::vector<double>& D_3 = current_cell.D_3;
            const std::vector<double>& D_4 = current_cell.D_4;

            std::vector<double> F1 = vector_subtract(FcDS_1, D_1);
            std::vector<double> F2 = vector_subtract(FcDS_2, D_2);
            std::vector<double> F3 = vector_subtract(FcDS_3, D_3);
            std::vector<double> F4 = vector_subtract(FcDS_4, D_4);

            current_cell.R = vector_add(vector_add(F1, F2), vector_add(F3, F4));

        }
    }
}


void SpatialDiscretization::run() {
    SpatialDiscretization::compute_dummy_cells();
    SpatialDiscretization::compute_Fc_DeltaS_Lambdac();
    SpatialDiscretization::compute_dissipation();
    SpatialDiscretization::compute_R();
}

