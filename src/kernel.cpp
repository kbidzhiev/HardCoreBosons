#include "kernel.hpp"

#include <vector>
#include <omp.h>
#include "profile.h"
#include <future>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "../../cpp_libs/eigen/Eigen/Dense"

using Eigen::MatrixXcd;
using boost::math::quadrature::trapezoidal;

double KF(){
	return M_PI * RHO;
}

double Mu_chempot(){
	return KF()*KF() * 0.5 / MASS;
}

double Energy(Q_momenta q_momenta){
	return (q_momenta.value) * (q_momenta.value) * 0.5 / MASS ;
}

double Tau(Q_momenta q_momenta,  SpaceTime spacetime){
	return spacetime.t * Energy(q_momenta) - spacetime.x * q_momenta.value;
}

double Theta (Q_momenta q_momenta){
	double result = exp(B_BETA * (Energy(q_momenta) - Mu_chempot())) + 1.0;
	return 1./result;
}

Cplx Eminus (Q_momenta q_momenta,  SpaceTime spacetime){
	Cplx result = sqrt(Theta (q_momenta)/M_PI);
	result *= exp(Cplx_i * 0.5 * Tau(q_momenta, spacetime));
	return result;
}

//Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime spacetime){
//	return 0;
//}


Cplx G0 (SpaceTime spacetime){
	Cplx result = exp(-Cplx_i * M_PI /4.);
	result *= sqrt(MASS / (2 * M_PI * spacetime.t));
	result *= exp( Cplx_i * 0.5 * MASS * spacetime.x * spacetime.x / spacetime.t );
	return result;
}

double Q (const size_t i) {
	size_t middle_point = g.abscissa().size();
	return i < middle_point ?
			- KF() * g.abscissa()[middle_point - i - 1] :
			  KF() * g.abscissa()[i - middle_point];
};

double Weight (const size_t i) {
	size_t middle_point = g.weights().size();
	return i < middle_point ?
			KF() * g.weights()[middle_point - i - 1] :
			KF() * g.weights()[i - middle_point];
};

pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime spacetime){

	const size_t s = 2 * g.weights().size();
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

	vector<Cplx> e_minus(s);
	vector<Cplx> e_infty(s);
	vector<Cplx> e_infty_derivative(s);

	const Cplx first_term_of_PVderivative = - Cplx_i * spacetime.t * 2.0 * M_PI * G0(spacetime) / MASS;

#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		const Cplx tau = Tau (Q_momenta(Q(i)), spacetime);
		const Cplx pv = PrincipalValue(Q_momenta(Q(i)), spacetime);
		e_infty[i] = (pv/M_PI - Lambda * exp(- Cplx_i * tau))/(Lambda * Lambda + 1);
		e_minus[i] = Eminus(Q_momenta(Q(i)), spacetime) ;
		const Cplx second_term_of_PVderivative = - Cplx_i * (Q(i) * spacetime.t / MASS - spacetime.x) * pv;
		Cplx e_inf_der = (first_term_of_PVderivative + second_term_of_PVderivative)/M_PI;
		e_inf_der += -Lambda * exp(-Cplx_i * tau) *(-Cplx_i) * (spacetime.t * Q(i) / MASS - spacetime.x);
		e_inf_der *= 1.0 /( Lambda * Lambda + 1);
		e_infty_derivative[i] = e_inf_der;
	}


#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			Cplx v;
			Cplx w = e_infty[i] * e_minus[i] //E_+ = E_infty E_-
					* e_infty[j] * e_minus[j] ;
			if(abs(w)< 1e-15) {w = Cplx(0, 0);}
			W(i, j) = sqrt(Weight(i)) * w * sqrt(Weight(j));
			if (i == j) {
				v = e_infty_derivative[i] * e_minus[i] * e_minus[i];
				if(abs(v)< 1e-15){ v = Cplx(0, 0);}
				V(i, i) = sqrt(Weight(i)) * v * sqrt(Weight(i)) + 1.0;
			} else {
				v = (e_infty[i] - e_infty[j]) * e_minus[i] * e_minus[j];
				v /= Q(i) - Q(j);
				if(abs(v)< 1e-16){ v = Cplx(0, 0);}
				V(i, j) = sqrt(Weight(i)) * v * sqrt(Weight(j));
				V(j, i) = V(i, j);
				W(j, i) = W(i, j);
			}

		}
	}

	Cplx detV = V.determinant();
	W = V - W;
	Cplx detW = W.determinant();

	return {detV, detW};
}



Cplx Grep(SpaceTime st){
	auto f = [&](double Lambda){
			auto [detv,detw] = Determinants(Lambda, st);
			Cplx coeff = G0 (st) * 2.0/(Lambda * Lambda + 1.0) ;
			detv *= (coeff - 1.0);
			return detv + detw;
	};

	Cplx result = trapezoidal(f, -150.0, 150.0);
	return result;
}



