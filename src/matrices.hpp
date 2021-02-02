#pragma once

#include <complex>
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "integral_kernel.hpp"


using Eigen::MatrixXcd;




pair<MatrixXcd,MatrixXcd> OnePlusV_W (const double eta,
		const double x_coordinate,
		const double t_time,
		const double momenta_truncation);


Cplx G_inf (const double x_coordinate, const double t_time);
