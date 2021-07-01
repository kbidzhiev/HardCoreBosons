#include "kernel.hpp"

#include <vector>
#include <omp.h>



const size_t s = 2 * g_l.weights().size() ; //Gauss
//const size_t s = 2 * g_l.weights().size() - 1; //GKrondrod, we always have EVEN number of weights


using Eigen::MatrixXcd;
using boost::math::quadrature::trapezoidal;

double KF_l(){
	return M_PI * RHO;
}

double Mu_chempot_l(){
	return 0;
	//return KF()*KF() * 0.5 / MASS;
}

double Energy_l(Q_momenta q_momenta){//verified
	return -2.0 * cos(q_momenta.value);
}

double Tau_l(Q_momenta q_momenta,  SpaceTime st){//verified
	return st.t * Energy_l(q_momenta) - st.x * q_momenta.value;
}

double Theta_l(Q_momenta q_momenta){//verified

	if (MAGN_FIELD < -1){
		return 1;
	}
	// exp (710.5) gives overflow and exp(-700) is underflow!
	double result = 1.0 + exp(2 * B_BETA * MAGN_FIELD)
			+ exp(B_BETA * (Energy_l(q_momenta) - Mu_chempot_l() + MAGN_FIELD));
	return 1./result;

//	double result = exp(B_BETA * (Energy_l(q_momenta) - Mu_chempot_l()));
//	result += 2 * cosh(B_BETA * MAGN_FIELD); // cosh (710.5) gives overflow !
//	result = 1./result;
//	result *= exp(-B_BETA * MAGN_FIELD);
//	return result;
}

Cplx Eminus_l(Q_momenta q_momenta,  SpaceTime st){//verified
	return exp(Cplx_i * 0.5 * Tau_l(q_momenta, st));
}

Cplx EPV_l (Q_momenta k_external,  SpaceTime st){//verified

	auto f = [&](double q) {
		Cplx exp_external = exp(-Cplx_i * Tau_l(k_external, st));
		Cplx exp_internal = exp(-Cplx_i * Tau_l(Q_momenta(q), st));
		return (exp_internal-exp_external) * cos(0.5*(q-k_external.value))/sin(0.5*(q-k_external.value));
	};
	Cplx result = trapezoidal(f, -M_PI, M_PI);//gauss<double, 30>::integrate(f, -M_PI, M_PI);
	return result/(2.0*M_PI);
}

Cplx EPV_derivative_l (Q_momenta k_momenta,  SpaceTime st){
	//verified numerically. Indeed it is a derivative of E_PV_l
	auto f = [&](double q) {
		Cplx exp_external = Cplx_i *(-2.0 * st.t * sin(k_momenta.value) + st.x) *
				exp(-Cplx_i * Tau_l(k_momenta, st));
		Cplx exp_internal = Cplx_i *(-2.0 * st.t * sin(q) + st.x) *
				exp(-Cplx_i * Tau_l(Q_momenta(q), st));
		return (exp_internal-exp_external) *
				cos(0.5*(q-k_momenta.value)) / sin(0.5*(q-k_momenta.value));
	};
	Cplx result = trapezoidal(f, -M_PI, M_PI);//gauss<double, 30>::integrate(f, -M_PI, M_PI);
	return result/(2.0*M_PI);
}


Cplx Eplus_l (Q_momenta q_momenta,  SpaceTime st){
	return EPV_l(q_momenta, st) * Eminus_l(q_momenta, st);
}

Cplx Lplus_l (double eta, Q_momenta q_momenta,  SpaceTime st){
	Cplx result = (1.0 - cos(eta)) * Eplus_l(q_momenta, st);
	result += (sin(eta) / Eminus_l(q_momenta, st));
	result *= 0.5 * sqrt(Theta_l(q_momenta));
	return result;
}
Cplx Lminus_l (Q_momenta q_momenta,  SpaceTime st){
	return Eminus_l(q_momenta, st) * sqrt(Theta_l(q_momenta));
}

Cplx G_l (SpaceTime st){ //verified
	return exp(Cplx_i * 0.5 * M_PI * st.x) * cyl_bessel_j(abs(st.x), 2.0 * st.t);
}

double Gamma_l(){
	if (MAGN_FIELD < -1){
		return 1;
	}

	return 1.0 + exp(2 * MAGN_FIELD * B_BETA);
}

double F_l(double eta){//verified
	// I've found a simpler form
	double result = Gamma_l() * Gamma_l() - 1.0;
	result /= Gamma_l() * Gamma_l() + 1.0 - 2.0 * Gamma_l() * cos(eta);
	return result;
}


Cplx Rplus_l (double eta, Q_momenta q_momenta, Q_momenta k_momenta,  SpaceTime st){
//	Cplx first = (1.0 - cos(eta)) * Eplus_l(q_momenta, st) * Eplus_l(k_momenta, st);
//	Cplx second = sin(eta) * (
//			Eplus_l(q_momenta, st)/Eminus_l(k_momenta, st) + Eplus_l(k_momenta, st) * Eminus_l(q_momenta, st)
//			);
//	Cplx third = (1.0 + cos(eta)) / (Eminus_l(q_momenta, st) * Eminus_l(k_momenta, st));
//	Cplx result = (first + second + third )/ (4.0 * M_PI);
//	return sqrt(Theta_l(q_momenta)) * result * sqrt(Theta_l(k_momenta));
	Cplx Eq = Eplus_l(q_momenta, st);
	Cplx Ek = Eplus_l(k_momenta, st);
	Cplx Emq = Eminus_l(q_momenta, st);
	Cplx Emk = Eminus_l(k_momenta, st);
	Cplx product = Emk * Emq * Ek * Eq;

	Cplx result = 1.0 + cos(eta) + product - product * cos(eta);
	result += (Ek * Emk + Eq * Emq) * sin(eta);
	result /= 4.0 * M_PI * Emk * Emq;
	return sqrt(Theta_l(q_momenta)) * result * sqrt(Theta_l(k_momenta));
}

Cplx Nu_diagonal_l(double eta, Q_momenta q_momenta, SpaceTime st) { //verified
	Cplx cos_term = 0.5 * (1.0 - cos(eta)) * Eminus_l(q_momenta, st)
			* Eminus_l(q_momenta, st) * EPV_derivative_l(q_momenta, st) / M_PI;
	Cplx sin_term = 0.5 * sin(eta) * Cplx_i
			* (2.0 * st.t * sin(q_momenta.value) - st.x) / M_PI;
	return (cos_term - sin_term) * Theta_l(q_momenta);
}

Cplx Nu_matrix_elem(double eta, Q_momenta k_momenta, Q_momenta k_prime_momenta, SpaceTime st) {
	Cplx nu_matrix_elem  = Lplus_l(eta, k_momenta, st)*Lminus_l(k_prime_momenta, st)
			- Lplus_l(eta, k_prime_momenta, st)*Lminus_l(k_momenta, st);
	nu_matrix_elem /= 2.0 * M_PI * tan(0.5 * (k_momenta.value - k_prime_momenta.value )); // Diagonal 0 is here 0

	return nu_matrix_elem;
}


double Q_l(const size_t i) {
	//size_t middle_point = g.abscissa().size() - 1;
	return i < g.abscissa().size() - 1 ?
			-M_PI * g.abscissa()[g.abscissa().size() - i - 1] :
			M_PI * g.abscissa()[i - g.abscissa().size() + 1];
}

double Weight_l (const size_t i) {
	//size_t middle_point = g.weights().size() - 1;
	return i < g.weights().size() - 1?
			M_PI * g.weights()[g.weights().size() -1 - i ] :
			M_PI * g.weights()[i - g.weights().size() + 1];
}

double Q_G_l (const size_t i) {
	//size_t middle_point = g.abscissa().size();
	return i < g_l.abscissa().size() ?
			-M_PI * g_l.abscissa()[g_l.abscissa().size() - i - 1] :
			M_PI * g_l.abscissa()[i - g_l.abscissa().size()];
}

double Weight_G_l (const size_t i) {
	//size_t middle_point = g.weights().size();
	return i < g_l.weights().size() ?
			M_PI * g_l.weights()[g_l.weights().size() - i - 1] :
			M_PI * g_l.weights()[i - g_l.weights().size()];
}

pair <Cplx, Cplx> Determinants_l(double eta, SpaceTime st){

	vector<Cplx> l_plus;
	vector<Cplx> l_minus;
	l_plus.reserve(s);
	l_minus.reserve(s);

	MatrixXcd One_plus_gammaQ(s, s);
	MatrixXcd One_plus_gammaQ_minus_gammaR(s, s);



//#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		Q_momenta q_i(Q_G_l(i));
		l_plus[i]  = Lplus_l(eta, q_i, st);
		l_minus[i] = Lminus_l(q_i, st);
	}


//#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			Cplx q_matrix_elem;
			Q_momenta q_i(Q_G_l(i));
			Q_momenta k_j(Q_G_l(j));
			Cplx r_plus_elem = Gamma_l() * Rplus_l(eta, q_i, k_j, st) ;
			Cplx r_minus_elem = 0.5 * (1.0 - cos(eta)) * G_l(st) * l_minus[i] * l_minus[j] / (2.0 * M_PI);

			if (i == j) {
				q_matrix_elem = Nu_diagonal_l(eta, q_i, st) - r_minus_elem;
				q_matrix_elem *= Gamma_l();
				One_plus_gammaQ(i, i) = 1.0 + Weight_G_l(i) * q_matrix_elem  ;
				One_plus_gammaQ_minus_gammaR(i, i) = 1.0 + Weight_G_l(i) * ( q_matrix_elem - r_plus_elem);

			} else {
				q_matrix_elem  = l_plus[i]*l_minus[j] - l_plus[j]*l_minus[i];
				q_matrix_elem /= 2.0 * M_PI * tan(0.5 * (Q_G_l(i) - Q_G_l(j)));
				q_matrix_elem -= r_minus_elem;
				q_matrix_elem *= Gamma_l();

				One_plus_gammaQ(i, j) = sqrt(Weight_G_l(i)) * q_matrix_elem * sqrt(Weight_G_l(j));
				One_plus_gammaQ_minus_gammaR(i, j) = sqrt(Weight_G_l(i)) *
						(q_matrix_elem - r_plus_elem) * sqrt(Weight_G_l(j));
				One_plus_gammaQ(j, i) = One_plus_gammaQ(i, j);
				One_plus_gammaQ_minus_gammaR(j, i) = One_plus_gammaQ_minus_gammaR(i, j);
			}
		}
	}

	Cplx detQ = One_plus_gammaQ.determinant();
	Cplx detR = One_plus_gammaQ_minus_gammaR.determinant();

	return {detQ, detR};
}


Cplx GrepEta_l(double eta,  SpaceTime st){
	auto [detq,detr] = Determinants_l(eta, st);
	detq *= G_l(st) - 1.0;
	return exp(Cplx_i * st.t * (Mu_chempot_l() - MAGN_FIELD)) * (detq + detr);
}

Cplx Grep_l(SpaceTime st){
	auto f= [&](double eta){
		return F_l(eta) * GrepEta_l(eta, st);
	};
	double error;
	Cplx result = g_integration.integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
	//Cplx result = gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);

	return result/(2.0 * M_PI);
}


