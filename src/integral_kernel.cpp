#include "integral_kernel.hpp"
#include <omp.h>
#include "profile.h"

#include <cstdlib>
#include <complex.h> // trigoniometric functions sin() cos()
#include "gauss.hpp" // Gauss-Legendre quadrature
#include "gauss_kronrod.hpp"
//#include "trapezoidal.hpp"


using namespace boost::math::quadrature;
//using boost::math::quadrature::trapezoidal;
using namespace std;

using Cplx = complex<double>; //pseudoname for complex<double>
const Cplx Cplx_i = Cplx(0,1);

const double mass = 1.0;
const double g_coupling = 999.;
const double b_beta = 1000.0;
const double rho = 0.5 ; // rho = #particles / system size;    so 0 < rho < 1
const double chem_potential = (M_PI * rho * M_PI * rho)/ (2 * mass);


const gauss<double, 10> g; // has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points
const double gauss_limits = 1.0; // gauss PV integration domain (-gauss_limits : gauss_limits)
double TRUNC = 1e-4; // controls convergence for interval from "gauss_limits" to \infty
double convergence = 1e-6; //gauss_kronrod convergence criteria

const double step = 1.0; //integration domain (i : i + step)




double Energy (const double q_momenta){
	return  q_momenta * q_momenta * 0.5 / mass;
}

double Tau (const double q_momenta, const double x_coordinate, const double t_time){
	return t_time * Energy(q_momenta) - x_coordinate * q_momenta;
}

Cplx PrincipalValue(const double q_momenta, const double x_coordinate, const double t_time) {

	auto f = [&](const double &p_momenta){
		/*
		 * Hilbert transformation
		 * https://en.wikipedia.org/wiki/Hilbert_transformation
		 */
		const Cplx exp_q_plus_p  = exp(-Cplx_i * Tau(q_momenta + p_momenta, x_coordinate, t_time));
		const Cplx exp_q_minus_p = exp(-Cplx_i * Tau(q_momenta - p_momenta, x_coordinate, t_time));
		return exp_q_plus_p - exp_q_minus_p;
	};

	const auto x = [&](const size_t i) {
		size_t middle_point = g.abscissa().size();
		return i < middle_point ?
				- g.abscissa()[middle_point - i - 1] :
				 g.abscissa()[i - middle_point];
	};

	const auto weight = [&](const size_t i) {
		size_t middle_point = g.weights().size();
		return i < middle_point ?
				g.weights()[middle_point - i - 1] :
				g.weights()[i - middle_point];
	};

	Cplx value_pole = 0;
	const int NUMBER_OF_POINTS = 2 * g.weights().size();
	for (int i = 0; i < NUMBER_OF_POINTS; i++) {
		value_pole += (weight(i) / x(i))
				* (f(gauss_limits * x(i) ) - f(0));
	}
	value_pole *= 0.5; // symmetrization of integration limits requires factor 1/2

	auto u = [&](const double &t){
		return f(t + gauss_limits) / (t + gauss_limits);
	};
	Cplx left_and_right = 0;
	complex<long double > df = 1.0 + Cplx_i;

	double trunc = TRUNC;

	for (double i = 0; abs(df) > trunc ;  ){
		//df = trapezoidal(u, i, i + step, convergence);
		df = gauss_kronrod<double, 15>::integrate(u, i, i+step, convergence);
		left_and_right += df;
		i += step;
	}

	return value_pole + left_and_right;
}



Cplx PrincipalValueDerivative(const double q_momenta, const double x_coordinate, const double t_time) {
	//G_DURATION("Principal Value");
	auto f = [&](const double &p_momenta){
		/*
		 * Hilbert transformation
		 * https://en.wikipedia.org/wiki/Hilbert_transform
		 */
		double momenta = q_momenta + p_momenta;
		Cplx exp_q_plus_p_deriv = -Cplx_i * exp(-Cplx_i * Tau(momenta, x_coordinate, t_time));
		exp_q_plus_p_deriv *= (t_time * momenta / mass - x_coordinate);

		momenta = q_momenta - p_momenta;
		Cplx exp_q_minus_p_deriv = -Cplx_i * exp(-Cplx_i * Tau(momenta, x_coordinate, t_time));
		exp_q_minus_p_deriv *= (t_time * momenta / mass - x_coordinate);

		return exp_q_plus_p_deriv - exp_q_minus_p_deriv;
	};

	const auto x = [&](const size_t i) {
		size_t middle_point = g.abscissa().size();
		return i < middle_point ?
				-g.abscissa()[middle_point - i - 1] :
				 g.abscissa()[i - middle_point];
	};

	const auto weight = [&](const size_t i) {
		size_t middle_point = g.weights().size();
		return i < middle_point ?
				g.weights()[middle_point - i - 1] :
				g.weights()[i - middle_point];
	};

	Cplx value_pole = 0;
	const int NUMBER_OF_POINTS = 2 * g.weights().size();
	for (int i = 0; i < NUMBER_OF_POINTS; i++) {
		value_pole += (weight(i) / x(i))
				* (f(gauss_limits * x(i) ) - f(0));
	}
	value_pole *= 0.5; // symmetrization of integration limits requires factor 1/2

	auto u = [&](const double &t){
		return f(t + gauss_limits) / (t + gauss_limits);
	};
	Cplx left_and_right = 0;
	complex<long double > df = 1.0 + Cplx_i;

	const double trunc = TRUNC;

	for (double i = 0; abs(df) > trunc ;  ){
		//df = trapezoidal(u, i, i + step, convergence);
		df = gauss_kronrod<double, 15>::integrate(u, i, i+step, convergence);
		left_and_right += df;
		i += step;
		cout << "q = " << i << " integr(q) = " << left_and_right <<  endl;
	}

//	while (abs(df) > trunc ){
//		df = trapezoidal(u, double(i), double(i + step));
//		left_and_right += df;
//		i += step;
//	}
	//cout << "PV " << value_pole + left_and_right << endl;
	cout << "pole = " << value_pole <<endl
			<< "rest = " << left_and_right << endl;
	return value_pole + left_and_right;
}


Cplx E_inf(const double eta, const double q_momenta, const double x_coordinate, const double t_time){
	const double tau = Tau(q_momenta, x_coordinate, t_time);
	Cplx result = PrincipalValue(q_momenta, x_coordinate, t_time)/M_PI;
	result *= sin(0.5 * eta) * sin(0.5 * eta);
	result += sin(0.5 * eta) * cos(0.5 * eta) * exp(-Cplx_i * tau);
	//cout << "result = " << result << endl;
	return result;
}

Cplx E_inf_Derivative(const double eta, const double q_momenta, const double x_coordinate, const double t_time){
	const double tau = Tau(q_momenta, x_coordinate, t_time);
	Cplx result = PrincipalValueDerivative(q_momenta, x_coordinate, t_time)/M_PI;
	result *= sin(0.5 * eta) * sin(0.5 * eta);
	result += sin(0.5 * eta) * cos(0.5 * eta) * exp(-Cplx_i * tau) * (-Cplx_i) * (t_time * q_momenta / mass - x_coordinate);
	//cout << "result = " << result << endl;
	return result;
}

double Theta(const double q_momenta){
	const double theta = 1.0/(exp(b_beta*(Energy(q_momenta) - chem_potential)) + 1.0);
	//cout << "Theta = " << theta << endl;
	return  theta;
}

double UpperMomenta(){
	// we return value there Theta < 10^15
	double result  = chem_potential;
	result += 15 * log(10) /b_beta; //15 is exponent of 10^15
	result *= 2 * mass;
	return sqrt(result);
}

Cplx E_minus(const double q_momenta, const double x_coordinate, const double t_time ){
	return exp(Cplx_i * 0.5 * Tau(q_momenta, x_coordinate, t_time)) * sqrt(Theta(q_momenta)/ M_PI) ;
}

Cplx G_0(const double x_coordinate, const double t_time){
	Cplx result = exp(-Cplx_i * M_PI /4.);
	result *= sqrt(mass / (2 * M_PI * t_time));
	result *= exp( Cplx_i * 0.5 * mass * x_coordinate * x_coordinate / t_time );
	return result;
}



// XY model :

//Cplx Weight (const double momenta,
//		const double b_beta,
//		const double gamma,
//		const double magnetization){
//
//	const Cplx sqrt_term = sqrt(pow(magnetization-cos(momenta),2) + pow(gamma * sin(momenta),2));
//
//	Cplx result =
//			magnetization - cos(momenta) - Cplx_i * gamma * sin(momenta) ;
//
//	result /= sqrt_term ;
//	result *= -tanh(0.5*b_beta * sqrt_term);
//	result += 1;
//	result *= 0.5;
//	return result;
//}
//
//
//Cplx Kernel (const Cplx weight,
//		const double x,
//		const double momenta1,
//		const double momenta2) {
//	Cplx kernel = sin(0.5 * x * (momenta1-momenta2))/
//			sin(0.5 * (momenta1-momenta2));
//	Cplx result = -weight/M_PI;
//	if (abs (momenta1 - momenta2) > 1e-14 ){
//		result *= kernel;
//	} else {
//		result *= x;
//	}
//	return  result;
//}
//
//Cplx KernelFiniteRank (const Cplx kernel,
//		const Cplx weight,
//		const double x,
//		const double momenta1,
//		const double momenta2
//		){
//	Cplx result = kernel;
//	result += weight * exp(-Cplx_i * x * 0.5* (momenta1 + momenta2))
//		* exp(Cplx_i * 0.5* (momenta1 - momenta2))/M_PI;
//	return result;
//}












