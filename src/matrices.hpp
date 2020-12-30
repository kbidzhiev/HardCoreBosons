#pragma once

#include <complex>
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "integral_kernel.hpp"


using Eigen::MatrixXd;
using Eigen::MatrixXcd;





MatrixXcd IdMatrix(int size);

MatrixXcd V(const double eta ,const double x_coordinate, const double t_time);

