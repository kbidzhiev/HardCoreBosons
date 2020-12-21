//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "boost/math/quadrature/gauss.hpp"
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"

#include <fstream>	 // file input output
#include <iostream>	 // screen input output
#include <cmath> 	 // pow (x,3) = x^3 = x*x*x and M_PI = pi = 3.14
#include <complex.h> // trigoniometric functions sin() cos()

using Eigen::MatrixXd;
using Eigen::MatrixXcf;
using Eigen::MatrixXcd;
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
			sin(0.5 * (momenta1-momenta2));
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
		* exp(Cplx_i * 0.5* (momenta1 - momenta2))/M_PI;
	return result;
}



MatrixXcd ConstructMatrix(
							const double x,
							const double beta,
							const double gamma,
							const double magnetization,
							const bool finite_rank = false){
	const int s = 10;
	gauss<double, 10> g;
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
					cerr << "finite rank" << endl;
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
	if (abs(coordinate + 10.) < 1e-10){
		cout << m << endl;
	}
	MatrixXcd m_finite = ConstructMatrix(coordinate, 100, 0, 0, true);

	return m_finite.determinant() - m.determinant();
}

int main(){
	double momenta1 = -2.5;
	double momenta2 = 1.;
	double beta = 2.9;
	double gamma = -0.3;
	double magnetization = 0.4;
	double coordinate = -10;

	auto weight = Weight (momenta1, beta, gamma, magnetization);
	auto kern = Kernel(weight, coordinate, momenta1, momenta2);
	cout << "weight = " << weight << endl;
	cout << "kernel = " << kern << endl;
	cout << ConstructMatrix(coordinate, beta,gamma,magnetization ) << endl;

	{
		LOG_DURATION("total");
		ofstream correlator; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)
		correlator.open("Correlator_100.dat", mode);
		correlator.precision(15);
		correlator << "#x \t correlator \t time \n";

		// ------- Correlator profile -------
		double dx = 0.01;
		double system_size = 1.0;
		const int n_steps = system_size / dx;
		for (int n = -n_steps / 2; n <= n_steps / 2; ++n) {
			const double coordinate = n * dx; //+param.val("time_shift");
		//  cout << "#x" << n << "/" << n_steps << "\ttime=" << time << endl;


		//	for (int time = 0; time <= time_total; time++) {
		//		correlator << "\"t=" << dt*time << "\"" << endl;

//				const complex<double> determ = det(coordinate);
//				correlator << coordinate << "\t" << real(determ) << "\t" << imag(determ) << "\t" << endl;
		//	}
		//	correlator << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)



		}
		//cout << det(-10.);
	}

}
