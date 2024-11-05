#ifndef SPATIALDISCRETIZATION_H
#define SPATIALDISCRETIZATION_H

#include "cell.h"
#include <vector>

class SpatialDiscretization {
public:


    std::vector<std::vector<double>> OMEGA;
    std::vector<std::vector<std::vector<std::vector<double>>>> s;
    std::vector<std::vector<std::vector<double>>> Ds;
    std::vector<std::vector<std::vector<std::vector<double>>>> n;
    std::vector<std::vector<std::vector<double>>> W;
    std::vector<std::vector<std::vector<double>>> R_c;
    std::vector<std::vector<std::vector<double>>> R_d;
    std::vector<std::vector<std::vector<double>>> R_d0;
    std::vector<std::vector<std::vector<std::vector<double>>>> flux;
    std::vector<std::vector<std::vector<std::vector<double>>>> D;
    std::vector<std::vector<std::vector<std::vector<double>>>> eps_2;
    std::vector<std::vector<std::vector<std::vector<double>>>> eps_4;
    std::vector<std::vector<double>> Lambda_I;
    std::vector<std::vector<double>> Lambda_J;

    std::vector<std::vector<double>> x, y;
    double rho, u, v, E, T, p;
    double T_ref, U_ref;
    int ny, nx;

    std::vector<std::vector<cell>> domain_cells;
    std::vector<std::vector<cell>> solid_wall_cells;
    std::vector<std::vector<cell>> farfield_cells;
    std::vector<std::vector<cell>> cells;

    SpatialDiscretization(const std::vector<std::vector<double>>& x,
                          const std::vector<std::vector<double>>& y,
                          const double& rho,
                          const double& u,
                          const double& v,
                          const double& E,
                          const double& T,
                          const double& p,
                          const double& T_ref,
                          const double& U_ref);

    std::tuple<double, double, double, double, double, double> conservative_variable_from_W(const std::vector<double>& W) const;

    void compute_dummy_cells();

    std::vector<double> FcDs(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) const;

    double Lambdac(const std::vector<double>& W, const std::vector<double>& n, const double& Ds) const;

    void compute_Fc_DeltaS();
    std::tuple<double, double> compute_epsilon(const cell& cell_Im1, const cell& cell_I,
                                              const cell& cell_Ip1, const cell& cell_Ip2,
                                              double k2 = 0.5, double k4 = 1.0/32.0) const;

    void compute_dissipation();

    void compute_R_c();

    void compute_R_d();

    void run_odd();

    void run_even();
};



#endif //SPATIALDISCRETIZATION_H
