
#include "integral_kernel.hpp"


#include <cstdlib>
#include <complex.h> // trigoniometric functions sin() cos()
#include "boost/math/quadrature/gauss.hpp" // Gauss-Legendre quadrature
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>


#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>



using namespace boost::math::quadrature;
using boost::math::quadrature::trapezoidal;
using namespace std;

using Cplx = complex<double>; //pseudoname for complex<double>
const Cplx Cplx_i = Cplx(0,1);

const double mass = 1.0;
const double g_coupling = 999.;
const double b_beta = 100.;
const double chem_potential = 0;

const double Energy (const double q_momenta){

	return  q_momenta * q_momenta * 0.5 / mass;
}

const double Tau (const double q_momenta, const double x_coordinate, const double t_time){
	return t_time * Energy(q_momenta) - x_coordinate * q_momenta;
}

Cplx PrincipalValue(const double q_momenta, const double x_coordinate, const double t_time) {

	auto ExpTau = [&](const double &p_momenta) {
		return exp(-Cplx_i * Tau(p_momenta, x_coordinate, t_time));
	};

	auto f = [&](const double &p_momenta){
		return (ExpTau(q_momenta + p_momenta) - ExpTau(q_momenta - p_momenta));
	};

	auto u = [&](const double &t){
		//double z = t/(1.0 - t);
		//double Jacobian = 1.0 / ( (1. - t)*(1.0 - t));
		return f(t) / t;
	};

	const double cutoff = 1.0;
	const gauss<double, 100> g;

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
				* (f(cutoff * x(i) ) - f(0));
	}
	value_pole *= 0.5; // symmetrization of integration limits requires factor 1/2

	Cplx left_and_right = 0;
	complex<long double > df = 1.0 + Cplx_i;
	double i = cutoff;
	const double step = 1.;
	while (abs(df) > 1e-10 ){
		df = trapezoidal(u, double(i), double(i + step));
		left_and_right += df;
		//cout << "step = " << i/step << " real = " << real(df) << " imag = " << imag(df)  << endl;
		i += step;
	}
	return value_pole + left_and_right;
}


const Cplx E_inf(const double eta, const double q_momenta, const double x_coordinate, const double t_time){
	const double tau = Tau(q_momenta, x_coordinate, t_time);
	Cplx result = PrincipalValue(q_momenta, x_coordinate, t_time)/M_PI;
	result *= sin(0.5 * eta) * sin(0.5 * eta);
	result += sin(0.5 * eta) * cos(0.5 * eta) * exp(-Cplx_i * tau);
	return 0;
	//this part is wrong !! it doesnt integrate from -inf to inf
}

const double Teta(const double q_momenta){
	return exp(b_beta*(Energy(q_momenta) - chem_potential)) + 1 ;
}

const Cplx E_minus(const double q_momenta, const double x_coordinate, const double t_time ){
	return sqrt(Teta(q_momenta)/ M_PI)*exp(Cplx_i * 0.5 * Tau(q_momenta, x_coordinate, t_time));
}

const Cplx E_plus(const double eta, const double q_momenta, const double x_coordinate, const double t_time){
	return E_inf(eta, q_momenta, x_coordinate, t_time)
			* E_minus(q_momenta, x_coordinate, t_time);
}


const Cplx V_p_q_inf(const double p_momenta, const double q_momenta,
		const double eta, const double x_coordinate, const double t_time){

	const Cplx e_plus_p = E_plus(eta, p_momenta, x_coordinate, t_time);
	const Cplx e_plus_q = E_plus(eta, q_momenta, x_coordinate, t_time);

	const Cplx e_minus_p = E_minus(p_momenta, x_coordinate, t_time);
	const Cplx e_minus_q = E_minus(q_momenta, x_coordinate, t_time);

	return (e_plus_p * e_minus_q - e_plus_q * e_minus_p)/(p_momenta- q_momenta);
}

const Cplx W_p_q_inf(const double p_momenta,
		const double q_momenta,
		const double eta,
		const double x_coordinate,
		const double t_time){

	const Cplx e_plus_p = E_plus(eta, p_momenta, x_coordinate, t_time);
	const Cplx e_plus_q = E_plus(eta, q_momenta, x_coordinate, t_time);

	return e_plus_p * e_plus_q * 0.5 / (sin(0.5* eta) * sin(0.5* eta));
}




Cplx G_0(const double x_coordinate, const double t_time){
	Cplx result = exp(-Cplx_i * M_PI /4.);
	result *= sqrt(1./(2 * M_PI * t_time));
	result *= exp( (Cplx_i * x_coordinate * x_coordinate )/ (2 * t_time) );
	return result;
}

Cplx Weight (const double momenta,
		const double b_beta,
		const double gamma,
		const double magnetization){

	const Cplx sqrt_term = sqrt(pow(magnetization-cos(momenta),2) + pow(gamma * sin(momenta),2));

	Cplx result =
			magnetization - cos(momenta) - Cplx_i * gamma * sin(momenta) ;

	result /= sqrt_term ;
	result *= -tanh(0.5*b_beta * sqrt_term);
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












