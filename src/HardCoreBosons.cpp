//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"


#include <iostream>

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


	auto f = [](const double& t) { return t * t * std::atan(t); };
	double Q = gauss<double, 100>::integrate(f, 0, 1);
	std::cout << Q ;
}
