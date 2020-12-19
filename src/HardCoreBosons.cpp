//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"


#include <iostream>
#include <cmath> 	 // pow (x,3) = x^3 = x*x*x and M_PI = pi = 3.14
#include <complex.h> // trigoniometric functions sin() cos()

using Eigen::MatrixXd;
using namespace boost::math::quadrature;

using Cplx = complex<double>; //pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0,1);


complex<double> Weight (const double momenta,
		const double beta,
		const double gamma,
		const double magnetization){

	const Cplx sqrt_term = sqrt(pow(magnetization-cos(momenta),2) + pow(gamma * sin(momenta),2));

	complex<double> result =
			magnetization - cos(momenta) - Cplx_i * gamma * sin(momenta) ;

	result /= sqrt_term ;
	result *= -tanh(beta * sqrt_term /2.);
	result += 1;
	result *= 0.5;


	return result;
}


Cplx Kernel (const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2
		){
	Cplx kernel = sin(0.5 * x * (momenta1-momenta2))/
			sin(0.5 * x * (momenta1-momenta2));
	Cplx result = -weight/M_PI;
	if (abs (momenta1 - momenta2) > 1e-14 ){
		result *= kernel;
	} else {
		result *= x;
	}
	return  result;
}

Cplx KernelFiniteRank (const Cplx kernel,
		const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2
		){
	Cplx result = kernel;
	result += weight * exp(-Cplx_i * x * 0.5* (momenta1 + momenta2))
		*exp(-Cplx_i * 0.5* (momenta1 - momenta2));
	return result;
}




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


	double Q{0};
	{
		LOG_DURATION("gauss  30");
		Q = gauss<double, 30>::integrate(f, -M_PI, M_PI);
	}
	{
		LOG_DURATION("gauss  40");
		Q = gauss<double, 40>::integrate(f, -M_PI, M_PI);
	}
	{
		LOG_DURATION("gauss 100");
		Q = gauss<double, 100>::integrate(f, -M_PI, M_PI);
	}
	{
		LOG_DURATION("gauss 500");
		Q = gauss<double, 500>::integrate(f, -M_PI, M_PI);
	}



	std::cout << std::setprecision(16)
		<< Q  << "\n"
		<< 1.7724381183457067 ;
}
