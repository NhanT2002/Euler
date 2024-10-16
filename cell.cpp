#include "cell.h"
#include <iostream>
#include <vector>
#include <cmath>

cell::cell(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double rho, double u, double v, double E)
    : x1(x1), y1(y1), x2(x2), y2(y2), x3(x3), y3(y3), x4(x4), y4(y4),
      rho(rho), u(u), v(v), E(E),
      s1(2), s2(2), s3(2), s4(2),
      W(4), FcDS_1(4), FcDS_2(4), FcDS_3(4), FcDS_4(4),
      D_1(4), D_2(4), D_3(4), D_4(4),
      R(4) {

    // Cell volume
    OMEGA = 0.5*((x1-x3)*(y2-y4) + (x4-x2)*(y1-y3));

    // Initialize side vectors
    s1[0] = y2 - y1; s1[1] = x1 - x2; // s1 = [dy, -dx]
    s2[0] = y3 - y2; s2[1] = x2 - x3;
    s3[0] = y4 - y3; s3[1] = x3 - x4;
    s4[0] = y1 - y4; s4[1] = x4 - x1;

    // Compute norms
    Ds1 = std::hypot(s1[0], s1[1]);
    Ds2 = std::hypot(s2[0], s2[1]);
    Ds3 = std::hypot(s3[0], s3[1]);
    Ds4 = std::hypot(s4[0], s4[1]);

    // Normal vectors
    n1 = {s1[0] / Ds1, s1[1] / Ds1};
    n2 = {s2[0] / Ds2, s2[1] / Ds2};
    n3 = {s3[0] / Ds3, s3[1] / Ds3};
    n4 = {s4[0] / Ds4, s4[1] / Ds4};

    // Conservative variable initialization
    W = {rho, rho * u, rho * v, rho * E};

    // Initialize Lambda properties
    Lambda_1_I = Lambda_1_J = Lambda_2_I = Lambda_2_J = 0.0;
    Lambda_3_I = Lambda_3_J = Lambda_4_I = Lambda_4_J = 0.0;
    Lambda_1_S = Lambda_2_S = Lambda_3_S = Lambda_4_S = 0.0;

    // Initialize epsilon properties
    eps2_2 = eps4_2 = eps2_3 = eps4_3 = eps2_4 = eps4_4 = eps2_1 = eps4_1 = 0.0;

}

