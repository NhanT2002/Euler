//
// Created by hieun on 10/15/2024.
//

#ifndef TEMPORALDISCRETIZATION_H
#define TEMPORALDISCRETIZATION_H


#include <string>
#include "SpatialDiscretization.h"
#include <vector>

class TemporalDiscretization{
public:
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
    double rho, u, v, E, T, p;
    double T_ref, U_ref;

    SpatialDiscretization current_state;

    TemporalDiscretization(const std::vector<std::vector<double>>& x,
                           const std::vector<std::vector<double>>& y,
                           const double& rho,
                           const double& u,
                           const double& v,
                           const double& E,
                           const double& T,
                           const double& p,
                           const double& T_ref,
                           const double& U_ref);

    double compute_dt(const cell& cell_IJ, double sigma=0.5) const;

    static std::vector<double> compute_L2_norm(const std::vector<std::vector<std::vector<double>>> &residuals);

    static void save_checkpoint(const std::vector<std::vector<std::vector<double>>>& q,
                         const std::vector<int>& iteration,
                         const std::vector<std::vector<double>>& Residuals,
                         const std::string& file_name = "checkpoint.txt");

    static std::tuple<std::vector<std::vector<std::vector<double>>>,std::vector<int>,
           std::vector<std::vector<double>>> load_checkpoint(const std::string& file_name);

    std::tuple<std::vector<std::vector<std::vector<double>>>,
               std::vector<std::vector<std::vector<double>>>,
               std::vector<std::vector<double>>> RungeKutta(int it_max = 20000);

    void run();
};



#endif //TEMPORALDISCRETIZATION_H
