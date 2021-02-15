#include "kernel.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

const gauss<double, 20> gPV;// Use only even size; g has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points


const double gauss_limits = 1.0;
double TRUNC = 1e-3; // controls convergence for interval from "gauss_limits" to \infty
double convergence = 1e-4; //gauss_kronrod convergence criteria
const double step = 1.0; //integration domain (i : i + step)

using namespace boost::math::quadrature;


double Ksi (const size_t i){
	size_t middle_point = gPV.abscissa().size();
		return i < middle_point ?
			- gPV.abscissa()[middle_point - i - 1] :
 			  gPV.abscissa()[i - middle_point];
}


double WeightPV (const size_t i) {
	size_t middle_point = gPV.weights().size();
	return i < middle_point ?
			gPV.weights()[middle_point - i - 1] :
			gPV.weights()[i - middle_point];
};

Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime st){



	auto f = [&](const double &p_momenta){
		/*
		 * Hilbert transformation
		 * https://en.wikipedia.org/wiki/Hilbert_transformation
		 */
		const Cplx exp_q_plus_p  = exp(-Cplx_i * Tau(Q_momenta(q_momenta.value + p_momenta), st));
		const Cplx exp_q_minus_p = exp(-Cplx_i * Tau(Q_momenta(q_momenta.value - p_momenta), st));
		return exp_q_plus_p - exp_q_minus_p;
	};


	Cplx value_pole = 0;
	const int NUMBER_OF_POINTS = 2 * gPV.weights().size();
	for (int i = 0; i < NUMBER_OF_POINTS; i++) {
		value_pole += (WeightPV(i) / Ksi(i))
				* (f(gauss_limits * Ksi(i) ) - f(0));
	}
	value_pole *= 0.5; // symmetrization of integration limits requires factor 1/2

	auto u = [&](const double &t){
		return f(t + gauss_limits) / (t + gauss_limits);
	};

	Cplx left_and_right = 0;
	complex<long double > df = 1.0 + Cplx_i;

	double trunc = TRUNC;


	for (double i = 0; abs(df) > trunc ;  ){
		df = gauss_kronrod<double, 31>::integrate(u, i, i + step, convergence); //15, 31, 41, 51 and 61
		left_and_right += df;
		i += step;

	}

	return value_pole + left_and_right;



}
