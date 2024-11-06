#include "SpatialDiscretization.h"
#include "vector_helper.h"

#include <array>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

template <typename T>
T combineBoundaryValues(const T& solidWall, const T& interior, const T& farfield) {
    T combined;
    combined.reserve(2 * solidWall.size() + interior.size() + 2 * farfield.size());

    // Insert `solidWall` twice
    combined.insert(combined.end(), solidWall.begin(), solidWall.end());
    combined.insert(combined.end(), solidWall.begin(), solidWall.end());

    // Insert `interior`
    combined.insert(combined.end(), interior.begin(), interior.end());

    // Insert `farfield` twice
    combined.insert(combined.end(), farfield.begin(), farfield.end());
    combined.insert(combined.end(), farfield.begin(), farfield.end());

    return combined;
}


SpatialDiscretization::SpatialDiscretization(const std::vector<std::vector<double>>& x,
                          const std::vector<std::vector<double>>& y,
                          const double& rho,
                          const double& u,
                          const double& v,
                          const double& E,
                          const double& T,
                          const double& p,
                          const double& T_ref,
                          const double& U_ref)
                              :x(x), y(y), rho(rho), u(u), v(v), E(E), T(T), p(p), T_ref(T_ref), U_ref(U_ref){
    ny = static_cast<int>(y.size());
    nx = static_cast<int>(x[0].size());

    std::vector OMEGA_domain(ny - 1, std::vector<double>(nx - 1));
    std::vector s_domain(ny - 1, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    std::vector Ds_domain(ny - 1, std::vector(nx - 1, std::vector<double>(2)));
    std::vector n_domain(ny - 1, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    std::vector W_domain(ny - 1, std::vector(nx - 1, std::vector<double>(4)));

    R_c.resize(ny - 1, std::vector(nx - 1, std::vector<double>(4)));
    R_d.resize(ny - 1, std::vector(nx - 1, std::vector<double>(4)));
    R_d0.resize(ny - 1, std::vector(nx - 1, std::vector<double>(4)));
    flux.resize(ny - 1 + 4, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    D.resize(ny - 1 + 4, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    eps_2.resize(ny - 1 + 4, std::vector(nx - 1,std::vector<double>(2)));
    eps_4.resize(ny - 1 + 4, std::vector(nx - 1,std::vector<double>(2)));
    Lambda_I.resize(ny - 1 + 4, std::vector<double>(nx - 1));
    Lambda_J.resize(ny - 1 + 4, std::vector<double>(nx - 1));
    Lambda_S.resize(ny - 1 + 4, std::vector(nx - 1, std::vector<double>(4)));
    for (size_t j = 0; j < ny - 1 ; ++j) {
        for (size_t i = 0; i < nx - 1; ++i) {
            const double& x1 = x[j][i];
            const double& x2 = x[j][i+1];
            const double& x3 = x[j+1][i+1];
            const double& x4 = x[j+1][i];
            const double& y1 = y[j][i];
            const double& y2 = y[j][i+1];
            const double& y3 = y[j+1][i+1];
            const double& y4 = y[j+1][i];

            // Calculate OMEGA
            OMEGA_domain[j][i] = 0.5 * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));

            // Set s and compute Ds using s values
            s_domain[j][i][0] = {y2 - y1, x1 - x2};
            s_domain[j][i][1] = {y1 - y4, x4 - x1};

            // Length of s vectors
            Ds_domain[j][i][0] = std::hypot(s_domain[j][i][0][0], s_domain[j][i][0][1]);
            Ds_domain[j][i][1] = std::hypot(s_domain[j][i][1][0], s_domain[j][i][1][1]);

            // Normal vectors
            n_domain[j][i][0] = {s_domain[j][i][0][0] / Ds_domain[j][i][0], s_domain[j][i][0][1] / Ds_domain[j][i][0]};
            n_domain[j][i][1] = {s_domain[j][i][1][0] / Ds_domain[j][i][1], s_domain[j][i][1][1] / Ds_domain[j][i][1]};

            // Compute W
            W_domain[j][i] = {rho, rho * u, rho * v, rho * E};
        }
    }

    std::vector OMEGA_solidWall(OMEGA_domain.begin(), OMEGA_domain.begin() + 1);
    std::vector OMEGA_farfield(OMEGA_domain.end() - 1, OMEGA_domain.end());

    std::vector s_solidWall(s_domain.begin(), s_domain.begin() + 1);
    std::vector s_farfield(s_domain.end() - 1, s_domain.end());
    std::vector Ds_solidWall(Ds_domain.begin(), Ds_domain.begin() + 1);
    std::vector Ds_farfield(Ds_domain.end() - 1, Ds_domain.end());
    std::vector n_solidWall(n_domain.begin(), n_domain.begin() + 1);
    std::vector n_farfield(n_domain.end() - 1, n_domain.end());
    std::vector W_solidWall(W_domain.begin(), W_domain.begin() + 1);
    std::vector W_farfield(W_domain.end() - 1, W_domain.end());

    for (size_t i = 0; i < nx - 1; ++i) {
        const double& x1 = x[ny-1][i];
        const double& x2 = x[ny-1][i+1];
        const double& y1 = y[ny-1][i];
        const double& y2 = y[ny-1][i+1];

        // Set s and compute Ds using s values
        s_farfield[0][i][0] = {y2 - y1, x1 - x2};

        // Length of s vectors
        Ds_farfield[0][i][0] = std::hypot(s_farfield[0][i][0][0], s_farfield[0][i][0][1]);

        // Normal vectors
        n_farfield[0][i][0] = {s_farfield[0][i][0][0] / Ds_farfield[0][i][0], s_farfield[0][i][0][1] / Ds_farfield[0][i][0]};
    }

    OMEGA.resize(ny -1 + 4, std::vector<double>(nx - 1));
    s.resize(ny - 1 + 4, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    Ds.resize(ny - 1 + 4, std::vector(nx - 1, std::vector<double>(2)));
    n.resize(ny - 1 + 4, std::vector(nx - 1, std::vector(2, std::vector<double>(2))));
    W.resize(ny - 1 + 4, std::vector(nx - 1, std::vector<double>(4)));

    // Combine using the helper function
    OMEGA = combineBoundaryValues(OMEGA_solidWall, OMEGA_domain, OMEGA_farfield);
    s = combineBoundaryValues(s_solidWall, s_domain, s_farfield);
    Ds = combineBoundaryValues(Ds_solidWall, Ds_domain, Ds_farfield);
    n = combineBoundaryValues(n_solidWall, n_domain, n_farfield);
    W = combineBoundaryValues(W_solidWall, W_domain, W_farfield);

    std::cout << "end" << std::endl;
}

void SpatialDiscretization::compute_dummy_cells() {
    // Solid wall
    for (int i = 0; i < nx - 1; ++i) {
        double p3, p4;
        auto [rho_val, u_val, v_val, E_val, T_val, p2] = conservative_variable_from_W(W[2][i]);
        std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p3) = conservative_variable_from_W(W[3][i]);
        std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p4) = conservative_variable_from_W(W[4][i]);

        const double pw = (15 * p2 - 10 * p3 + 3 * p4) / 8.0; // Blazek
        const double p1 = 2 * pw - p2;
        std::vector<double> vel = {u_val, v_val};

        std::vector<double> n1 = n[2][i][0];

        std::vector<std::vector<double>> R = { {-n1[1], n1[0]}, {n1[0], n1[1]} };
        const double q_t = -R[0][0] * vel[0] - R[0][1] * vel[1];
        const double q_n = -R[1][0] * vel[0] - R[1][1] * vel[1];

        const double y_eta = s[2][i][0][0] / Ds[2][i][0];
        const double x_eta = s[2][i][0][1] / Ds[2][i][0];

        // Swanson Turkel
        const double u_dummy = x_eta * q_t + y_eta * q_n;
        const double v_dummy = -y_eta * q_t + x_eta * q_n;


        E_val = p1 / (1.4 - 1) / rho_val + 0.5 * (u_dummy * u_dummy + v_dummy * v_dummy);

        W[0][i] = {rho_val, rho_val * u_dummy, rho_val * v_dummy, rho_val * E_val};
        W[1][i] = {rho_val, rho_val * u_dummy, rho_val * v_dummy, rho_val * E_val};
    }

    // Farfield
    for (int i = 0; i < nx - 1; ++i) {
        auto [rho_val, u_val, v_val, E_val, T_val, p_val] = conservative_variable_from_W(W[W.size()-3][i]);
        const double c = std::sqrt(1.4 * 287 * T_val*T_ref)/U_ref;
        const double M = std::sqrt(u_val * u_val + v_val * v_val) / c;
        std::vector<double> n3 = vector_scale(-1, n[n.size()-2][i][0]);

        if (u_val * n3[0] + v_val * n3[1] > 0) { // Out of cell
            if (M >= 1) {
                W[W.size()-2][i] = {rho_val, rho_val * u_val, rho_val * v_val, rho_val * E_val};
                W[W.size()-1][i] = {rho_val, rho_val * u_val, rho_val * v_val, rho_val * E_val};

            }
            else {  // Subsonic
                const double p_b = this->p;  // Boundary pressure
                const double rho_b = rho_val + (p_b - p_val) / (c * c);
                const double u_b = u_val + n3[0] * (p_val - p_b) / (rho_val * c);
                const double v_b = v_val + n3[1] * (p_val - p_b) / (rho_val * c);
                const double E_b = p_b / ((1.4 - 1) * rho_b) + 0.5 * (u_b * u_b + v_b * v_b);

                std::vector<double> W_b = {rho_b, rho_b * u_b, rho_b * v_b, rho_b * E_b};
                std::vector<double> W_a = vector_subtract(vector_scale(2, W_b), W[W.size()-3][i]);

                W[W.size()-2][i] = W_a;
                W[W.size()-1][i] = W_a;

            }
        }
        else {  // Moving into the cell
            if (M >= 1) {  // Supersonic
                W[W.size()-2][i] = {this->rho, this->rho * this->u, this->rho * this->v, this->rho * this->E};
                W[W.size()-1][i] = {this->rho, this->rho * this->u, this->rho * this->v, this->rho * this->E};

            } else {  // Subsonic
                const double p_b = 0.5 * (this->p + p_val - rho_val * c * (n3[0] * (this->u - u_val) + n3[1] * (this->v - v_val)));
                const double rho_b = this->rho + (p_b - this->p) / (c * c);
                const double u_b = this->u - n3[0] * (this->p - p_b) / (rho_val * c);
                const double v_b = this->v - n3[1] * (this->p - p_b) / (rho_val * c);
                const double E_b = p_b / ((1.4 - 1) * rho_b) + 0.5 * (u_b * u_b + v_b * v_b);

                std::vector<double> W_b = {rho_b, rho_b * u_b, rho_b * v_b, rho_b * E_b};
                std::vector<double> W_a = vector_subtract(vector_scale(2, W_b), W[W.size()-3][i]);

                W[W.size()-2][i] = W_a;
                W[W.size()-1][i] = W_a;
            }
        }
    }
}

// Define the conservative_variable_from_W function as per your requirements
std::tuple<double, double, double, double, double, double> SpatialDiscretization::conservative_variable_from_W(const std::vector<double>& W) const {
    // Implement the conversion from W to (rho, u, v, E)
    double rho = W[0];
    double u = W[1] / rho;
    double v = W[2] / rho;
    double E = W[3] / rho;
    double p = (1.4-1)*rho*(E-(u*u+v*v)/2);
    double T = p/(rho*287)*U_ref*U_ref/T_ref;
    return std::make_tuple(rho, u, v, E, T, p);
}

std::vector<double> SpatialDiscretization::FcDs(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) const {
    auto [rho, u, v, E, T, p] = conservative_variable_from_W(W);
    double V = n[0]*u + n[1]*v;
    double H = E + p/rho;

    return {rho*V*Ds, (rho*u*V + n[0]*p)*Ds, (rho*v*V + n[1]*p)*Ds, rho*H*V*Ds};
}

double SpatialDiscretization::Lambdac(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) const {
    auto [rho, u, v, E, T, p] = conservative_variable_from_W(W);
    double c = std::sqrt(1.4*287*T*T_ref)/U_ref;
    const double V = n[0]*u + n[1]*v;
    const double lambda = (std::abs(V) + c)*Ds;

    return lambda;
}

void SpatialDiscretization::compute_Fc_DeltaS() {
    const auto ny = W.size();
    const auto nx = W[0].size();

    for (int j = 2; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            std::vector<double> avg_W1 = vector_scale(0.5, vector_add(W[j][i], W[j - 1][i]));
            std::vector<double> avg_W4 = vector_scale(0.5, vector_add(W[j][i], W[j][(i - 1 + nx) % nx]));

            std::vector<double> FcDs_1 = FcDs(avg_W1, n[j][i][0], Ds[j][i][0]);
            std::vector<double> FcDs_4 = FcDs(avg_W4, n[j][i][1], Ds[j][i][1]);

            flux[j][i][0] = FcDs_1;
            flux[j][i][1] = FcDs_4;
        }
    }

}

std::tuple<double, double> SpatialDiscretization::compute_epsilon(const std::vector<double>& W_Im1, const std::vector<double>& W_I,
                                                                 const std::vector<double>& W_Ip1, const std::vector<double>& W_Ip2,
                                                                 double k2, double k4) const {
    // Retrieve pressure from the conservative variables (assuming the last element is pressure)
    double p_Im1, p_I, p_Ip1, p_Ip2;
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Im1) = conservative_variable_from_W(W_Im1);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_I) = conservative_variable_from_W(W_I);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Ip1) = conservative_variable_from_W(W_Ip1);
    std::tie(std::ignore, std::ignore,std::ignore, std::ignore, std::ignore, p_Ip2) = conservative_variable_from_W(W_Ip2);

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
    const auto ny = W.size();
    const auto nx = W[0].size();

    for (int j = 0; j < ny-1; ++j) {
        for (int i = 0; i < nx; ++i) {
            // Calculate Lambda values
            std::vector<double> n1_n3 = vector_scale(0.5, vector_add(n[j][i][0], n[j+1][i][0]));
            std::vector<double> n2_n4 = vector_scale(0.5, vector_add(n[j][(i + 1) % nx][1], n[j][i][1]));

            double ds2_plus_ds4 = 0.5 * (Ds[j][(i + 1) % nx][1] + Ds[j][i][1]);
            double ds1_plus_ds3 = 0.5 * (Ds[j][i][0] + Ds[j+1][i][1]);

            // Compute Lambda values
            Lambda_I[j][i] = Lambdac(W[j][i], n2_n4, ds2_plus_ds4);
            Lambda_J[j][i] = Lambdac(W[j][i], n1_n3, ds1_plus_ds3);
        }
    }

    for (int j = 2; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            std::vector<double>& W_IJ = W[j][i];
            std::vector<double>& W_Ip1J = W[j][(i + 1) % nx];
            std::vector<double>& W_IJp1 = W[j + 1][i];
            std::vector<double>& W_Im1J = W[j][(i - 1 + nx) % nx];
            std::vector<double>& W_IJm1 = W[j - 1][i];
            std::vector<double>& W_Im2J = W[j][(i - 2 + nx) % nx];
            std::vector<double>& W_IJm2 = W[j - 2][i];

            // Lambda calculations

            const double Lambda_4_I = 0.5 * (Lambda_I[j][i] + Lambda_I[j][(i - 1 + nx) % nx]);
            const double Lambda_4_J = 0.5 * (Lambda_J[j][i] + Lambda_J[j][(i - 1 + nx) % nx]);
            const double Lambda_4_S = Lambda_4_I + Lambda_4_J;

            const double Lambda_1_J = 0.5 * (Lambda_J[j][i] + Lambda_J[j - 1][i]);
            const double Lambda_1_I = 0.5 * (Lambda_I[j][i] + Lambda_I[j - 1][i]);
            const double Lambda_1_S = Lambda_1_I + Lambda_1_J;

            Lambda_S[j][i][0] = Lambda_1_S;
            Lambda_S[j][i][3] = Lambda_4_S;


            // Epsilon calculations
            auto[eps2_4, eps4_4] = compute_epsilon(W_Ip1J, W_IJ, W_Im1J, W_Im2J);
            auto[eps2_1, eps4_1] = compute_epsilon(W_IJp1, W_IJ, W_IJm1, W_IJm2);

            eps_2[j-2][i][0] = eps2_1;
            eps_2[j-2][i][1] = eps2_4;
            eps_4[j-2][i][0] = eps4_1;
            eps_4[j-2][i][1] = eps4_4;

            // Dissipation terms
            std::vector<double> D_1 = vector_scale(Lambda_1_S,
               vector_subtract(
                   vector_scale(eps2_1, vector_subtract(W_IJm1, W_IJ)),
                   vector_scale(eps4_1,
                       vector_add(
                            vector_subtract(W_IJm2, vector_scale(3, W_IJm1)),
                            vector_subtract(vector_scale(3, W_IJ),W_IJp1
                       )
                   ))
               )
            );

            std::vector<double> D_4 = vector_scale(Lambda_4_S,
                vector_subtract(
                    vector_scale(eps2_4, vector_subtract(W_Im1J, W_IJ)),
                    vector_scale(eps4_4,
                        vector_add(
                            vector_subtract(W_Im2J, vector_scale(3, W_Im1J)),
                            vector_subtract(vector_scale(3, W_IJ),W_Ip1J
                        )
                    ))
                )
            );

            D[j][i][0] = D_1;
            D[j][i][1] = D_4;

        }
    }

    // Boundary conditions
    for (int i = 0; i < nx; ++i) {

        // Calculate D_1 for cell (3, i)
        D[3][i][0] = vector_scale(Lambda_S[3][i][0],
            vector_subtract(
                vector_scale(eps_2[3][i][0], vector_subtract(W[2][i], W[3][i])),

                vector_scale(eps_4[4][i][0], vector_subtract(vector_subtract(
                    vector_scale(2.0, W[3][i]), W[2][i]), W[4][i]))
            )
        );

        // Calculate D_1 for cell (2, i)
        D[2][i][0] = vector_scale(Lambda_S[2][i][0],
            vector_subtract(
                vector_scale(eps_2[3][i][0], vector_subtract(W[2][i], W[3][i])),

                vector_scale(eps_4[4][i][0], vector_subtract(vector_subtract(
                    vector_scale(2.0, W[3][i]), W[2][i]), W[4][i]))
            )
    );
    }
}

void SpatialDiscretization::compute_R_c() {
    const auto ny = W.size();
    const auto nx = W[0].size();

    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {

            // Extract the flux vectors
            const std::vector<double>& FcDS_1 = flux[j][i][0];
            const std::vector<double>& FcDS_2 = vector_scale(-1, flux[j][(i + 1) % nx][1]);
            const std::vector<double>& FcDS_3 = vector_scale(-1, flux[j+1][i][0]);
            const std::vector<double>& FcDS_4 = flux[j][i][1];


            R_c[j-2][i] = vector_add(vector_add(FcDS_1, FcDS_2), vector_add(FcDS_3, FcDS_4));

        }
    }
}

void SpatialDiscretization::compute_R_d() {
    const auto ny = W.size();
    const auto nx = W[0].size();

    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 0; i < nx; ++i) {

            // Extract the dissipation vectors

            const std::vector<double>& D_1 = D[j][i][0];
            const std::vector<double>& D_2 = vector_scale(-1, D[j][(i + 1) % nx][1]);
            const std::vector<double>& D_3 = vector_scale(-1, D[j+1][i][0]);
            const std::vector<double>& D_4 = D[j][i][1];


            R_d[j-2][i] = vector_add(vector_add(D_1, D_2), vector_add(D_3, D_4));

        }
    }
}

void SpatialDiscretization::run_odd() {
    SpatialDiscretization::compute_dummy_cells();
    SpatialDiscretization::compute_Fc_DeltaS();
    SpatialDiscretization::compute_R_c();
}

void SpatialDiscretization::run_even() {
    SpatialDiscretization::compute_dummy_cells();
    SpatialDiscretization::compute_Fc_DeltaS();
    SpatialDiscretization::compute_dissipation();
    SpatialDiscretization::compute_R_c();
    SpatialDiscretization::compute_R_d();
}

