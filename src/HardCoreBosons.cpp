//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"


#include <iostream>	 // screen input output
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


MatrixXd ConstructMatrix( 	const int size,
							const double x,
							const double beta,
							const double gamma,
							const double magnetization){
	MatrixXd m(size,size);
	gauss<double, 5> g;

	auto momenta = [&](const size_t i){
		size_t middle_point = g.abscissa().size();
		bool size_parity_is_odd = g.abscissa().front() == 0 ? false : true ;
		if (size_parity_is_odd){
			return  i < middle_point ? -M_PI*g.abscissa()[middle_point - i - 1] : M_PI*g.abscissa()[i - middle_point];
		} else {
			return  i < (middle_point-1) ? -M_PI*g.abscissa()[middle_point - i - 1] : M_PI*g.abscissa()[i - middle_point + 1];
		}
	};

	auto weight = [&](const size_t i){
		size_t middle_point = g.weights().size();
		bool size_parity_is_odd = g.abscissa().front() == 0 ? false : true ;
		//cerr << size_parity_is_odd << endl;
		if (size_parity_is_odd){
			return  i < middle_point ? M_PI*g.weights()[middle_point - i - 1]
									 : M_PI*g.weights()[i - middle_point];
		} else {
			return  i < (middle_point-1) ? M_PI*g.weights()[middle_point - i - 1]
										 : M_PI*g.weights()[i - middle_point + 1];
		}

	};




	return m;
}

/*
 *  gauss<double, 10> g; in wolfram its
 *  GaussianQuadratureWeights[10, -1, 1]
 * 	g.weights() = {0.295524, 0.269267, 0.219086, 0.149451, 0.0666713} positive part
 * 	g.abscissa()= {0.148874, 0.433395, 0.67941,  0.865063, 0.973907 } positive part
 *
 *	to obtain GaussianQuadratureWeights[10, -pi, pi] one should
 *	pi*g.weights() and pi*g.abscissa()
 */


int main(){
//
//  MatrixXd m(2,2);
//  m(0,0) = 3;
//  m(1,0) = 2.5;
//  m(0,1) = -1;
//  m(1,1) = m(1,0) + m(0,1);
//  std::cout << m << std::endl;




//	std::cout << std::setprecision(16)
//		<< Q  << "\n"
//		<< 1.7724381183457067 ;
}
