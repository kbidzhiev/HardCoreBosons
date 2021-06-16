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

double Energy_l(Q_momenta q_momenta){
	return -2.0 * cos(q_momenta.value) ;
}

double Tau_l(Q_momenta q_momenta,  SpaceTime spacetime){
	return spacetime.t * Energy_l(q_momenta) - spacetime.x * q_momenta.value;
}

double Theta_l(Q_momenta q_momenta){
	double result = exp(B_BETA * (Energy_l(q_momenta) - Mu_chempot_l())) + 2.0;
	return 1./result;
}

Cplx Eminus_l(Q_momenta q_momenta,  SpaceTime spacetime){
	return exp(Cplx_i * 0.5 * Tau_l(q_momenta, spacetime));
}

Cplx EPV_l (Q_momenta q_momenta,  SpaceTime spacetime){

	auto f = [&](double x) {
		Cplx exp_external = exp(-Cplx_i * Tau_l(q_momenta, spacetime));
		Cplx exp_internal = exp(-Cplx_i * Tau_l(Q_momenta(x), spacetime));
		return (exp_internal-exp_external) * cos(0.5*(x-q_momenta.value))/sin(0.5*(x-q_momenta.value));
	};

	//cout << "trapezoidal" << trapezoidal(f, -M_PI, M_PI) << endl;
	//double error;
	//cout << "gauss_kronrod" << gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error) << endl;
	//cout << "gauss" << gauss<double, 30>::integrate(f, -M_PI, M_PI) << endl;

	//Cplx result = gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
	Cplx result = trapezoidal(f, -M_PI, M_PI);//gauss<double, 30>::integrate(f, -M_PI, M_PI);
	return result/(2.0*M_PI);
}

Cplx EPV_derivative_l (Q_momenta q_momenta,  SpaceTime spacetime){
	auto f = [&](double x) {
		Cplx exp_external = Cplx_i *(-2.0 * spacetime.t * sin(q_momenta.value) + spacetime.x) *
				exp(-Cplx_i * Tau_l(q_momenta, spacetime));
		Cplx exp_internal = Cplx_i *(-2.0 * spacetime.t * sin(x) + spacetime.x) *
				exp(-Cplx_i * Tau_l(Q_momenta(x), spacetime));
		return (exp_internal-exp_external) *
				cos(0.5*(x-q_momenta.value)) / sin(0.5*(x-q_momenta.value));
	};

	//cout << "trapezoidal" << trapezoidal(f, -M_PI, M_PI) << endl;
	//double error;
	//cout << "gauss_kronrod" << gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error) << endl;
	//cout << "gauss" << gauss<double, 30>::integrate(f, -M_PI, M_PI) << endl;

	//Cplx result = gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
	Cplx result = trapezoidal(f, -M_PI, M_PI);//gauss<double, 30>::integrate(f, -M_PI, M_PI);
	return result/(2.0*M_PI);
}


Cplx Eplus_l (Q_momenta q_momenta,  SpaceTime spacetime){
	return EPV_l(q_momenta, spacetime) * Eminus_l(q_momenta, spacetime);
}

Cplx Lplus_l (double eta, Q_momenta q_momenta,  SpaceTime spacetime){
	Cplx result = (1.0 - cos(eta)) * Eplus_l(q_momenta, spacetime);
	result += sin(eta) / Eminus_l(q_momenta, spacetime);
	return 0.5 * result * sqrt(Theta_l(q_momenta));
}
Cplx Lminus_l (Q_momenta q_momenta,  SpaceTime spacetime){
	return Eminus_l(q_momenta, spacetime) * sqrt(Theta_l(q_momenta));
}

Cplx G_l (SpaceTime spacetime){
	return exp(Cplx_i * 0.5 * M_PI * spacetime.x) * cyl_bessel_j(abs(spacetime.x), 2.0 * spacetime.t);
}

double Gamma_l(){
	return 2.0;
}

Cplx F_l(double eta){
	Cplx result = Gamma_l() / (Gamma_l() - exp(Cplx_i * eta));
	result += 1.0/(exp(Cplx_i * eta)*Gamma_l() -1.0);
	return result;
}


Cplx Rplus_l (double eta, Q_momenta q_momenta, Q_momenta k_momenta,  SpaceTime st){
	Cplx first = (1.0 - cos(eta)) * Eplus_l(q_momenta, st) * Eplus_l(k_momenta, st)  ;
	Cplx second = sin(eta) * (
			Eplus_l(q_momenta, st)/Eminus_l(k_momenta, st) + Eplus_l(k_momenta, st) * Eminus_l(q_momenta, st)
			);
	Cplx third = (1.0 + cos(eta)) / (Eminus_l(q_momenta, st) * Eminus_l(k_momenta, st));
	return (first + second + third )/ (4.0 * M_PI);
}

Cplx Nu_diagonal_l(double eta, Q_momenta q_momenta, SpaceTime spacetime) {
	Cplx cos_term = 0.5 * (1.0 - cos(eta)) * Eminus_l(q_momenta, spacetime)
			* Eminus_l(q_momenta, spacetime) * EPV_derivative_l(q_momenta, spacetime) / M_PI;
	Cplx sin_term = 0.5 * sin(eta) * Cplx_i
			* (2.0 * spacetime.t * sin(q_momenta.value) - spacetime.x) / M_PI;
	Cplx result = (cos_term - sin_term) * Theta_l(q_momenta);
	result -= 0.5 * (1.0 - cos(eta)) * G_l(spacetime) * Lminus_l(q_momenta, spacetime) * Lminus_l(q_momenta, spacetime)
			/ (2.0 * M_PI);
	result *= Gamma_l();

	return result;
}


double Q_l(const size_t i) {
	//size_t middle_point = g.abscissa().size() - 1;
	return i < g.abscissa().size() - 1 ?
			-KF() * g.abscissa()[g.abscissa().size() - i - 1] :
			KF() * g.abscissa()[i - g.abscissa().size() + 1];
}

double Weight_l (const size_t i) {
	//size_t middle_point = g.weights().size() - 1;
	return i < g.weights().size() - 1?
			KF() * g.weights()[g.weights().size() -1 - i ] :
			KF() * g.weights()[i - g.weights().size() + 1];
}

double Q_G_l (const size_t i) {
	//size_t middle_point = g.abscissa().size();
	return i < g_l.abscissa().size() ?
			- KF() * g_l.abscissa()[g_l.abscissa().size() - i - 1] :
			  KF() * g_l.abscissa()[i - g_l.abscissa().size()];
}

double Weight_G_l (const size_t i) {
	//size_t middle_point = g.weights().size();
	return i < g_l.weights().size() ?
			KF() * g_l.weights()[g_l.weights().size() - i - 1] :
			KF() * g_l.weights()[i - g_l.weights().size()];
}

pair <Cplx, Cplx> Determinants_l(double eta, SpaceTime spacetime){

	vector<Cplx> l_plus;
	vector<Cplx> l_minus;
	l_plus.reserve(s);
	l_minus.reserve(s);

	MatrixXcd One_plus_gammaQ(s, s);
	MatrixXcd One_plus_gammaQ_minus_gammaR(s, s);

//#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		Q_momenta q_i(Q_G_l(i));
		l_plus[i]  = Lplus_l(eta, q_i, spacetime);
		l_minus[i] = Lminus_l(q_i, spacetime);
	}

//#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			Cplx q_matrix_elem;
			Q_momenta q_i(Q_G_l(i));
			Q_momenta k_j(Q_G_l(j));
			Cplx r_matrix_elem = Gamma_l() * Rplus_l(eta, q_i, k_j, spacetime);

			if (i == j) {
//				Cplx cos_term = 0.5 * (1.0 - cos(eta)) * Eminus_l(q_i, spacetime) *
//						Eminus_l(q_i, spacetime) * EPV_derivative_l(q_i, spacetime) / M_PI;
//				Cplx sin_term =  0.5 * sin(eta) * Cplx_i *
//						(2.0 * spacetime.t * sin(q_i.value) - spacetime.x)/ M_PI;
//				q_matrix_elem = (cos_term - sin_term) * Theta_l(q_i);
//				q_matrix_elem -= 0.5 * (1.0 - cos(eta)) * G_l(spacetime)*l_minus[i]*l_minus[i]/(2.0 * M_PI);
//				q_matrix_elem *= Gamma_l();

				q_matrix_elem = Nu_diagonal_l(eta, q_i, spacetime);

				One_plus_gammaQ(i, i) = sqrt(Weight_G_l(i)) * q_matrix_elem * sqrt(Weight_G_l(i)) + 1.0;

				One_plus_gammaQ_minus_gammaR(i, i) = sqrt(Weight_G_l(i)) *
						(q_matrix_elem - r_matrix_elem) * sqrt(Weight_G_l(i)) + 1.0; //new part
			} else {
				q_matrix_elem  = l_plus[i]*l_minus[j] - l_minus[i]*l_plus[j];
				q_matrix_elem /= 2.0 * M_PI * tan(0.5 * (Q_G_l(i) - Q_G_l(j))); // Diagonal 0 is here 0
				q_matrix_elem -= 0.5 * (1.0-cos(eta)) * G_l(spacetime)*l_minus[i]*l_minus[j]/(2.0 * M_PI);
				q_matrix_elem *= Gamma_l();

				One_plus_gammaQ(i, j) = sqrt(Weight_G_l(i)) * q_matrix_elem * sqrt(Weight_G_l(j));
				One_plus_gammaQ_minus_gammaR(i, j) = sqrt(Weight_G_l(i)) *
						(q_matrix_elem - r_matrix_elem) * sqrt(Weight_G_l(j));

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
	detq *= G_l(st)  - 1.0;
	return detq + detr;
}

Cplx Grep_l(SpaceTime st){
	auto f= [&](double eta){
		return F_l(eta) * GrepEta_l(eta, st)/(2.0 * M_PI);
	};
	double error = 100;

	Cplx result = gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);

	return result;
}


