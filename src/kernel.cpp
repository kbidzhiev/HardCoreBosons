#include "kernel.hpp"

#include <vector>
#include <omp.h>

TotalDuration elements("Elements");
TotalDuration mf("MatrixFilling");
TotalDuration det("Determinant");

const size_t s = 2 * g.weights().size() ; //Gauss
//const size_t s = 2 * g.weights().size() - 1; //GKrondrod, we always have ODD number of weights


using Eigen::MatrixXcd;
//using boost::math::quadrature::trapezoidal;

double KF(){
	return M_PI * RHO;
}

double Mu_chempot(){
	return KF()*KF() * 0.5 / MASS;
}

double Energy(Q_momenta q_momenta){
	return (q_momenta.value) * (q_momenta.value) * 0.5 / MASS ;
}

complex<double> Tau(Q_momenta q_momenta,  SpaceTime spacetime){
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

Cplx Einf (double eta, Q_momenta q_momenta,  SpaceTime spacetime){
	Cplx result = sin(eta*0.5)*sin(eta*0.5)*PrincipalValue(q_momenta, spacetime)/M_PI;
	result += sin(eta*0.5)*cos(eta*0.5)*exp(-Cplx_i*Tau(q_momenta, spacetime));
	return result;
}

Cplx Eplus (double eta, Q_momenta q_momenta,  SpaceTime spacetime){
	return Eminus(q_momenta, spacetime) * Einf(eta, q_momenta, spacetime);
}

Cplx EplusW (double eta, Q_momenta q_momenta,  SpaceTime spacetime){
	Cplx result = sin(eta*0.5)*PrincipalValue(q_momenta, spacetime)/M_PI;
	result += cos(eta*0.5)*exp(-Cplx_i*Tau(q_momenta, spacetime));
	return Eminus(q_momenta, spacetime) * result;
}





Cplx G0 (SpaceTime spacetime){
	Cplx result = exp(-Cplx_i * 0.25 * M_PI );
	result *= sqrt(MASS / (2.0 * M_PI * spacetime.t));
	result *= exp( Cplx_i * 0.5 * MASS * spacetime.x * spacetime.x / spacetime.t );
	return result;
}

double Q_G (const size_t i) {
	//size_t middle_point = g.abscissa().size();
	return i < g.abscissa().size() ?
			- KF() * g.abscissa()[g.abscissa().size() - i - 1] :
			  KF() * g.abscissa()[i - g.abscissa().size()];
}

double Weight_G (const size_t i) {
	//size_t middle_point = g.weights().size();
	return i < g.weights().size() ?
			KF() * g.weights()[g.weights().size() - i - 1] :
			KF() * g.weights()[i - g.weights().size()];
}

double Q_Kr(const size_t i) {
	//size_t middle_point = g.abscissa().size() - 1;
	return i < g.abscissa().size() - 1 ?
		   - KF() * g.abscissa()[g.abscissa().size() - i - 1] :
			 KF() * g.abscissa()[i - g.abscissa().size() + 1];
}

double Weight_Kr(const size_t i) {
	//size_t middle_point = g.weights().size() - 1;
	return i < g.weights().size() - 1?
			KF() * g.weights()[g.weights().size() -1 - i ] :
			KF() * g.weights()[i - g.weights().size() + 1];
}

pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime spacetime){
	vector<Cplx> e_minus(s,0.0);
	vector<Cplx> e_infty(s,0.0);
	vector<Cplx> e_infty_derivative(s,0.0);
	vector<Cplx> e_plus_for_W(s,0.0);
	vector<Cplx> e_infty_eps(s,0.0);
	vector<Cplx> e_minus_eps(s,0.0);
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

//	cout << "Weight_Kr\tQabscissa\n";
//	for(size_t i = 0; i< s; ++i){
//		cout << Q_Kr(i) <<"\t" << Weight_Kr(i)<< "\n";
//	}
//	cout << endl;
//	terminate();

	const Cplx first_term_of_PVderivative = -Cplx_i * spacetime.t * 2.0 * M_PI * G0(spacetime) / MASS;

	double epsilon = 1E-6;

	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(elements);
		Q_momenta q_i(Q_Kr(i));
		Cplx tau = Tau (q_i, spacetime);
		Cplx pv = PrincipalValue(q_i, spacetime);

		Q_momenta q_i_eps(Q_Kr(i) + epsilon);
		Cplx tau_eps = Tau (q_i_eps, spacetime);
		Cplx pv_eps = PrincipalValue(q_i_eps, spacetime);

		double sinsin = 1.0/(1.0 + Lambda*Lambda);
		double sincos = -Lambda/(1.0 + Lambda*Lambda);
		double sin = sqrt(sinsin);
		double cos = sincos / sin;

		e_infty[i] = sinsin * (pv/M_PI) + sincos * exp(- Cplx_i * tau);
		e_minus[i] = Eminus(q_i, spacetime) ;
		//
		e_infty_eps[i] = sinsin * (pv_eps/M_PI)	+ sincos * exp(- Cplx_i * tau_eps);
		e_minus_eps[i] = Eminus(q_i_eps, spacetime) ;
		//
		Cplx second_term_of_PVderivative = - Cplx_i * (Q_Kr(i) * spacetime.t / MASS - spacetime.x) * pv;
		Cplx e_inf_der = sinsin * (first_term_of_PVderivative + second_term_of_PVderivative)/M_PI;
		e_inf_der += sincos * exp(-Cplx_i * tau) *(-Cplx_i) * (spacetime.t * Q_Kr(i) / MASS - spacetime.x);

		e_infty_derivative[i] = e_inf_der;

		// In order to avoid having 0/0 in W eq (3.25), we
		// cancel sin^2 (eta/2) in denominator and enumerator
		e_plus_for_W[i] = sin * (pv/M_PI) + cos*exp(- Cplx_i * tau);
		e_plus_for_W[i] *= e_minus[i];
	}

	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(mf);
		for (size_t j = i; j < s; j++) {
			Cplx v = e_minus_eps[i] * e_minus[j];
			Cplx w =  0.5 * e_plus_for_W[i] * e_plus_for_W[j];
			if(i == j){
//				v *= e_infty_derivative[i];
				Cplx de =  e_infty_eps[i] - e_infty[j];
				double dq =  Q_Kr(i) + epsilon - Q_Kr(j);
				v *= (de/dq);
				V(i, j) = 1.0 + sqrt(Weight_Kr(i)) * v * sqrt(Weight_Kr(j)) ;
				W(i, j) = 1.0 + sqrt(Weight_Kr(i)) * (v-w) * sqrt(Weight_Kr(j)) ;
			} else {
				Cplx de =  e_infty[i] - e_infty[j];
				double dq =  Q_Kr(i) - Q_Kr(j);
				v *= (de/dq);
				V(i, j) = sqrt(Weight_Kr(i)) * v * sqrt(Weight_Kr(j)) ;
				W(i, j) = sqrt(Weight_Kr(i)) * (v-w) * sqrt(Weight_Kr(j)) ;
				V(j,i) = V(i,j);
				W(j,i) = W(i,j);
			}
		}
	}

	Cplx detV = V.determinant();
	Cplx detW = W.determinant();

	double Jacobian_from_eta_to_lambda = (2.0)/(Lambda * Lambda + 1.0);
	detV *= Jacobian_from_eta_to_lambda;
	detW *= Jacobian_from_eta_to_lambda;

	return {detV, detW};
}


pair <Cplx, Cplx> Determinants_eta(double eta, SpaceTime spacetime){
	vector<Cplx> e_minus(s,0.0);
	vector<Cplx> e_infty(s,0.0);
	vector<Cplx> e_minus_eps(s,0.0);
	vector<Cplx> e_infty_eps(s,0.0);
	vector<Cplx> e_plus_W(s,0.0);
	vector<Cplx> e_plus(s,0.0);
	vector<Cplx> e_plus_eps(s,0.0);

	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

//	cout << "Weight_G\tQabscissa\n";
//	for(size_t i = 0; i< s; ++i){
//		cout << Q_G(i) <<"\t" << Weight_G(i)<< "\n";
//	}
//	cout << endl;
//	terminate();


	double eps = 1E-6;

	//const Cplx first_term_of_PVderivative = -Cplx_i * spacetime.t * 2.0 * M_PI * G0(spacetime) / MASS;
	for (size_t i = 0; i < s; i++) {
		Q_momenta q_i(Q_G(i));
		e_plus[i] = Eplus(eta,q_i, spacetime);
		e_minus[i] = Eminus(q_i, spacetime);
		e_plus_W[i] = EplusW(eta, q_i, spacetime);

		Q_momenta q_i_eps(Q_G(i) + eps);
		e_plus_eps[i] = Eplus(eta, q_i_eps, spacetime);
		e_minus_eps[i] = Eminus(q_i_eps, spacetime) ;
	}


	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			Cplx v = (e_plus_eps[i] * e_minus[j] - e_plus[j] * e_minus_eps[i])/(Q_G(i)+eps-Q_G(j));
			Cplx w =  0.5 * e_plus_W[i] * e_plus_W[j];
			if(i == j){
				V(i, j) = 1.0 + sqrt(Weight_G(i)) * v * sqrt(Weight_G(j)) ;
				W(i, j) = 1.0 + sqrt(Weight_G(i)) * (v-w) * sqrt(Weight_G(j)) ;
			} else {
				V(i, j) = sqrt(Weight_G(i)) * v * sqrt(Weight_G(j)) ;
				W(i, j) = sqrt(Weight_G(i)) * (v-w) * sqrt(Weight_G(j)) ;
				V(j,i) = V(i,j);
				W(j,i) = W(i,j);
			}
		}
	}
	Cplx detV = V.determinant();
	Cplx detW = W.determinant();
	//cout << detV << endl << endl << detW << endl;
	//terminate();

return {detV, detW};
}

Cplx GrepLambda(double Lambda,  SpaceTime st){
	auto [detv,detw] = Determinants(Lambda, st);
	detv *= (G0(st) - 1.0);
	return (detv + detw);
}

Cplx GrepEta(double eta,  SpaceTime st){
	auto [detv,detw] = Determinants_eta(eta, st);
	detv *= (G0(st) - 1.0);
	return (detv + detw);
}
Cplx Grep_new(SpaceTime st){
	auto f= [&](double Eta){
		return GrepEta(Eta, st);
	};
	double error = 100;
//	Cplx result = gauss_kronrod<double, 61>::integrate(f, -100, 100,  10, 1e-9, &error);
//	Cplx result = gauss_kronrod<double, 61>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
	Cplx result = trapezoidal(f, -M_PI, M_PI);
	return 0.5*result/M_PI;
}

Cplx Grep(SpaceTime st){
	auto f= [&](double Lambda){
		return GrepLambda(Lambda, st);
	};
	double error = 100;
	Cplx result = gauss_kronrod<double, 61>::integrate(f, -100, 100,  10, 1e-9, &error);
	//Cplx result = gauss_kronrod<double, 61>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
	return 0.5*result/M_PI;
}



