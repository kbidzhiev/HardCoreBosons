//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"


#include <iostream>
#include <math.h> // for definition of pi = 3.14...    M_PI

using Eigen::MatrixXd;
using namespace boost::math::quadrature;


int main(){
//
//  MatrixXd m(2,2);
//  m(0,0) = 3;
//  m(1,0) = 2.5;
//  m(0,1) = -1;
//  m(1,1) = m(1,0) + m(0,1);
//  std::cout << m << std::endl;


	//auto f = [](const double& t) { return t * t * std::atan(t); };
	auto f = [](const double& t) { return std::exp(-t * t); };


//	boost::multiprecision::cpp_bin_float_quad Q2 = gauss<boost::multiprecision::cpp_bin_float_quad, 20>::integrate(f2, -M_PI, M_PI);

	double Q = gauss<double, 10>::integrate(f, -M_PI, M_PI);


	std::cout << std::setprecision(16)
		<< Q  << "\n"
		<< 1.7724381183457067 ;
}
