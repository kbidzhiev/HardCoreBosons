#include "kernel.hpp"

#include <vector>
#include <omp.h>



//const size_t s = 2 * g.weights().size() ; //Gauss
const size_t s = 2 * g.weights().size() - 1; //GKrondrod, we always have EVEN number of weights


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
	return exp(Cplx_i * 0.5 * Tau(q_momenta, spacetime));
}

Cplx EPV_l (Q_momenta q_momenta,  SpaceTime spacetime){

	auto f = [&](double x) {
		Cplx exp_external = exp(-Cplx_i * Tau_l(q_momenta, spacetime));
		Cplx exp_internal = exp(-Cplx_i * Tau_l(Q_momenta(x), spacetime));
		return (exp_internal-exp_external)/tan(0.5*(x-q_momenta.value));
	};
	//cout << "tau_l " << Tau_l(q_momenta, spacetime) << endl;
	//cout << "trapezoidal" << trapezoidal(f, -M_PI, M_PI) << endl;
	double error;
	//cout << "gauss" << gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error) << endl;

	Cplx result = gauss_kronrod<double, 31>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);
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


//


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

pair <Cplx, Cplx> Determinants_l(double eta, SpaceTime spacetime){

	vector<Cplx> l_plus;
	vector<Cplx> l_minus;

	l_plus.reserve(s);
	l_minus.reserve(s);
	MatrixXcd Q_plus(s, s);
	MatrixXcd R_plus(s, s);


	{
//#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		Q_momenta q_i(Q_l(i));
		l_plus[i] = Lplus_l(eta, q_i, spacetime);
		l_minus[i] = Lminus_l(q_i, spacetime);
//		cout << "i=" << i << ", Theta = " << Theta_l(q_i);
//		cout << "i=" << i << ", E+ = " << Lplus_l(eta,q_i, spacetime) << endl;
//		cout << "i=" << i << ", E- = " << Lminus_l(q_i, spacetime) << endl;
	}
	//terminate();

	}
	//time for 2 core proc
	//avr 2125209 	parallel
	//avr 31592		consequential


	{
//#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {

		for (size_t j = i; j < s; j++) {
			Cplx q;
			Cplx r_plus = Gamma_l() * l_plus[i]*l_plus[j]/

					(M_PI * (1.0 - cos(eta))) ;
			//if(abs(r_plus)< 1e-16) {w = Cplx(0, 0);}
			//W(i, j) = sqrt(Weight(i)) * w * sqrt(Weight(j));
			if (i == j) {
				q = (0.5 - 0.5 * cos(eta))/ (sin(0.5*Q_l(i))*sin(0.5*Q_l(i))) -
						Cplx_i * sin(eta) * (spacetime.x + 2.0 * spacetime.t * sin (Q_l(i))) ;
				q *= Theta_l(Q_momenta(Q_l(i)));
				q /= 2.0 * M_PI ;
				q -= 0.5 * (1.0-cos(eta)) * G_l(spacetime)*l_minus[i]*l_minus[j]/(2.0 * M_PI);
				q *= Gamma_l();

				Q_plus(i, i) = sqrt(Weight_l(i)) *  q * sqrt(Weight_l(i)) + 1.0;

				R_plus(i, i) = sqrt(Weight_l(i)) * (q - r_plus) * sqrt(Weight_l(i)) + 1.0; //new part
			} else {
				q = l_plus[i]*l_minus[j] - l_minus[i]*l_plus[j];
				q /= 2.0 * M_PI * tan(Q_l(i) - Q_l(j));
				q-= 0.5 * (1.0-cos(eta)) * G_l(spacetime)*l_minus[i]*l_minus[j]/(2.0 * M_PI);
				q *= Gamma_l();

				Q_plus(i, j) = sqrt(Weight_l(i)) * q * sqrt(Weight_l(j));
				R_plus(i, j) = sqrt(Weight_l(i)) * (q - r_plus) * sqrt(Weight_l(j)); //new part

				Q_plus(j, i) = Q_plus(i, j);
				R_plus(j, i) = R_plus(i, j);
				cout << "q = " << q << " i = " << i << ", j = " << j  << endl;
				cout << "r = " << r_plus << " i = " << i << ", j = " << j  << endl;
			}
			//cout << "Q = " << Q_plus(i, j) << " i = " << i << ", j = " << j  << endl;
			//cout << "R = " << R_plus(i, j) << " i = " << i << ", j = " << j  << endl;
			//cout << eta << endl;
		}
		terminate();
	}
	}

	Cplx detQ;
	Cplx detR;

	{

		detQ = Q_plus.determinant();
		detR = R_plus.determinant();
	}

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


