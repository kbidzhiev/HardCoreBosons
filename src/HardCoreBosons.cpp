//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"
#include "integral_kernel.hpp"

#include <fstream>	 // file input output
#include <iostream>	 // screen input output
#include <cmath> 	 // pow (x,3) = x^3 = x*x*x and M_PI = pi = 3.14
					 // trigoniometric functions sin() cos()
#include <complex.h>
using Eigen::MatrixXd;
using Eigen::MatrixXcf;
using Eigen::MatrixXcd;
using namespace boost::math::quadrature;
//
//using Cplx = complex<double>; //pseudoname for complex<double>
//
//const Cplx Cplx_i = Cplx(0,1);



MatrixXcd ConstructMatrix(
							const double x,
							const double beta,
							const double gamma,
							const double magnetization,
							const bool finite_rank = false){
	const int s = 100;
	gauss<double, 100> g;
	MatrixXcd m(s,s);
	auto identity = MatrixXcd::Identity(s, s);
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
	m += identity;
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
complex<double> det(double coordinate) {
	MatrixXcd m = ConstructMatrix(coordinate, 100, 0, 0);

	MatrixXcd m_finite = ConstructMatrix(coordinate, 100, 0, 0, true);

	return m_finite.determinant() - m.determinant();
}

int main(){

	{
		LOG_DURATION("total");
		ofstream correlator; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)
		correlator.open("Correlator_100.dat", mode);
		correlator.precision(15);
		correlator << "#x \t correlator \t time \n";

		// ------- Correlator profile -------
		double dx = 0.1;
		double system_size = 20.0;
		const int n_steps = system_size / dx;
		for (int n = -n_steps / 2; n <= n_steps / 2; ++n) {
			const double coordinate = n * dx; //+param.val("time_shift");
			const complex<double> determ = det(coordinate);
			correlator << coordinate << "\t" << real(determ) << "\t" << imag(determ) << "\t" << endl;
		}
	}
}
