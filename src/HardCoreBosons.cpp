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
using Eigen::MatrixXcf;
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
	result *= -tanh(0.5*beta * sqrt_term);
	result += 1;
	result *= 0.5;


	return result;
}


Cplx Kernel (const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2) {
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
		*exp(-Cplx_i * 0.5* (momenta1 - momenta2))/M_PI;
	return result;
}


MatrixXcf ConstructMatrix(
							const double x,
							const double beta,
							const double gamma,
							const double magnetization,
							const bool finite_rank = false){
	//MatrixXd m(5,5);
	MatrixXcf m(50,50);
	gauss<double, 50> g;
	bool size_parity_is_odd = g.abscissa().front() == 0 ? false : true ;

	auto momenta = [&](const size_t i){
		size_t middle_point = g.abscissa().size();
		if (size_parity_is_odd){
			return  i < middle_point ? -M_PI*g.abscissa()[middle_point - i - 1]
									  : M_PI*g.abscissa()[i - middle_point];
		} else {
			return  i < (middle_point-1) ? -M_PI*g.abscissa()[middle_point - i - 1]
										  : M_PI*g.abscissa()[i - middle_point + 1];
		}
	};

	auto weight = [&](const size_t i){
		size_t middle_point = g.weights().size();
		if (size_parity_is_odd){
			return  i < middle_point ? M_PI*g.weights()[middle_point - i - 1]
									 : M_PI*g.weights()[i - middle_point];
		} else {
			return  i < (middle_point-1) ? M_PI*g.weights()[middle_point - i - 1]
										 : M_PI*g.weights()[i - middle_point + 1];
		}

	};

	size_t length = size_parity_is_odd ? 2 * g.abscissa().size()
			: (2 * g.abscissa().size() -1);

	for (size_t i = 0; i<length; i++){
		for (size_t j = 0; j< length; j++){
				Cplx TD_weight  = Weight (momenta(i), beta, gamma, magnetization);
				Cplx kernel = Kernel (TD_weight, x, momenta(i), momenta(j) );
				if (finite_rank == true){
					kernel = KernelFiniteRank (kernel , TD_weight, x, momenta(i), momenta(j));
				}
				m(i,j) = sqrt(weight(i)) * kernel * sqrt(weight(j));
			}
	}
	m +=  MatrixXcf::Identity(50, 50);
	return  m;
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
	{LOG_DURATION("N = 50");

	MatrixXcf m = ConstructMatrix( 1., 100, 0, 0 );

	MatrixXcf m_finite = ConstructMatrix( 1., 100, 0, 0, true );
	//cout << m_finite << endl;
	auto det = m.determinant() - m_finite.determinant();
	cout << det << endl;
	}

}
