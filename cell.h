#ifndef CELL_H
#define CELL_H

#include <vector>

class cell {
public:
    // Constructor
    cell(double x1, double y1, double x2, double y2,
         double x3, double y3, double x4, double y4,
         double rho, double u, double v, double E);

    // Default constructor
    cell() : x1(0), y1(0), x2(0), y2(0), x3(0), y3(0), x4(0), y4(0),
             rho(0), u(0), v(0), E(0) {
    }
    // Member variables
    double x1, y1, x2, y2, x3, y3, x4, y4;
    double rho, u, v, E;

    double OMEGA{};

    std::vector<double> s1, s2, s3, s4;
    double Ds1{}, Ds2{}, Ds3{}, Ds4{};
    std::vector<double> n1, n2, n3, n4;

    std::vector<double> W;

    std::vector<double> FcDS_1, FcDS_2, FcDS_3, FcDS_4;

    double Lambda_1_I{}, Lambda_1_J{};
    double Lambda_2_I{}, Lambda_2_J{};
    double Lambda_3_I{}, Lambda_3_J{};
    double Lambda_4_I{}, Lambda_4_J{};

    double Lambda_1_S{}, Lambda_2_S{}, Lambda_3_S{}, Lambda_4_S{};

    double eps2_2{}, eps4_2{};
    double eps2_3{}, eps4_3{};
    double eps2_4{}, eps4_4{};
    double eps2_1{}, eps4_1{};

    std::vector<double> D_1, D_2, D_3, D_4;

    std::vector<double> R;

};



#endif //CELL_H
