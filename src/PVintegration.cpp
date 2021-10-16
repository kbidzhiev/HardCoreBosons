#include "kernel.hpp"


const gauss<double, 100> gPV;// Use only even size == 2n;
//g has pre-computed tables of abscissa and weights for 20 and 30 points


const double gauss_limits = 1.0;
double TRUNC = 1e-9; // controls convergence for interval from "gauss_limits" to \infty
double convergence = 1e-9; //gauss_kronrod convergence criteria
const double step = 1.; //integration domain (i : i + step)

using namespace boost::math;
//using namespace boost::math::quadrature;



Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime st){
	Cplx result  =  M_PI * Cplx_i * exp( -Cplx_i * Tau(q_momenta, st));//was - sign here
	Cplx erf_arg = (st.x - q_momenta.value * st.t) * (1.0 - Cplx_i) / (2.0 * sqrt(st.t));
	result *= (Faddeeva::erf(erf_arg));
	return result;
}



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
}

Cplx PrincipalValue_old(Q_momenta q_momenta,  SpaceTime st){
	/* What should one learn:
	 * Gauss quadrature needs HUGE matrices to converge (~100).
	 * Calculations in mathematica suggest that also
	 *
	 * to improve a quality of the old method one should work on
	 * trapezoidal implementation: increase number of a COUNTER
	 * and decrease convergence and truncation parameters
	 *
	 * or change trapezoidal by smth better
	 */


	auto f = [&](const double &p_momenta){ //is checked
		/*
		 * part of the Hilbert transformation
		 * https://en.wikipedia.org/wiki/Hilbert_transformation
		 */
		const Cplx exp_q_plus_p  = exp(-Cplx_i * Tau(Q_momenta(q_momenta.value + p_momenta), st));
		const Cplx exp_q_minus_p = exp(-Cplx_i * Tau(Q_momenta(q_momenta.value - p_momenta), st));
		return exp_q_plus_p - exp_q_minus_p;
	};


	// First we evaluate PV integral in the [-1,1] vicinity of the pole
	Cplx value_pole = 0;
	const int NUMBER_OF_POINTS = 2 * gPV.weights().size();
	for (int i = 0; i < NUMBER_OF_POINTS; i++) {
		value_pole += (WeightPV(i) / Ksi(i))
				* (f(gauss_limits * Ksi(i) ) - f(0));
	}
	value_pole *= 0.5; // symmetrization of integration limits around pole requires factor 1/2


	auto u = [&](const double &t){
		return f(t ) / (t );
	};
	Cplx left_and_right = 0;
	complex<double > df = 1.0 + Cplx_i;
	double trunc = TRUNC;
	int counter = 0;
	for (double i = gauss_limits; counter < 100 ; i += step){
		//df = gauss_kronrod<double, 61>::integrate(u, i, i + step, convergence); //15, 31, 41, 51 and 61
		df =  trapezoidal(u, i, i+step , convergence);
		left_and_right += df;
		if (abs(df) < trunc) {
			++counter;
			cout << counter << endl;
		};
	}
	return value_pole + left_and_right;
}
