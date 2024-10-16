#ifndef SPATIALDISCRETIZATION_H
#define SPATIALDISCRETIZATION_H

#include "cell.h"
#include <vector>

class SpatialDiscretization {
public:
    std::vector<std::vector<cell>> domain_cells;
    std::vector<std::vector<cell>> solid_wall_cells;
    std::vector<std::vector<cell>> farfield_cells;
    std::vector<std::vector<cell>> cells;

    std::vector<std::vector<double>> x, y;
    double rho, u, v, E, T, p;
    int ny, nx;

    SpatialDiscretization(const std::vector<std::vector<double>>& x,
                          const std::vector<std::vector<double>>& y,
                          const double& rho,
                          const double& u,
                          const double& v,
                          const double& E,
                          const double& T,
                          const double& p);

    static std::tuple<double, double, double, double, double, double> conservative_variable_from_W(const std::vector<double>& W);

    void compute_dummy_cells();

    static std::vector<double> FcDs(const std::vector<double>& W, const std::vector<double>& n, const double& Ds);

    static double Lambdac(const std::vector<double>& W, const std::vector<double>& n, const double& Ds);

    void compute_Fc_DeltaS();

    static std::tuple<double, double> compute_epsilon(const cell& cell_Im1, const cell& cell_I,
                                              const cell& cell_Ip1, const cell& cell_Ip2,
                                              double k2 = 0.5, double k4 = 1.0/64.0);

    void compute_dissipation();

    void compute_R_c();

    void compute_R_d();

    void run_odd();

    void run_even();
};



#endif //SPATIALDISCRETIZATION_H
