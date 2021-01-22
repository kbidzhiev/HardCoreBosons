#include "matrices.hpp"
#include "profile.h"

#include "gauss.hpp"
#include <utility>
#include <omp.h>

#include <vector>

#include <exception>

using namespace std;
using namespace boost::math::quadrature;

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::MatrixXi;
using Cplx = complex<double>;
//pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0, 1);

gauss<double, 6> g;




MatrixXcd IdMatrix(int size) {
	return MatrixXcd::Identity(size, size);
}


pair<MatrixXcd, MatrixXcd> OnePlusV_W(const double &eta,
		const double &x_coordinate, const double &t_time,
		const double & kL, const double & kR) {

	const int s = 6;
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

	//test;
//	const MatrixXi A = MatrixXi::Random(2,2);
//	cout << "automatic \n" << A << endl;
//	cout << "automatic \n" << A << endl;
//	cout << "automatic \n" << A << endl;
//
//	cout << "by hand \n";
//	for (int i = 0 ; i < 2; i++){
//		for (int j = 0 ; j < 2; j++){
//			cout << A(i,j) << " ";
//		}cout << endl;
//	}
//	terminate();


	const auto q = [&](const size_t i) {
		size_t middle_point = g.abscissa().size();
		const double jacobian = (kR - kL) * 0.5;
		const double shift = (kR + kL) * 0.5;
		return i < middle_point ?
				shift - jacobian * g.abscissa()[middle_point - i - 1] :
				shift + jacobian * g.abscissa()[i - middle_point];
	};

	const auto weight = [&](const size_t i) {
		size_t middle_point = g.weights().size();
		const double jacobian = (kR - kL) * 0.5;
		return i < middle_point ?
				jacobian * g.weights()[middle_point - i - 1] :
				jacobian * g.weights()[i - middle_point];
	};

//	for (int i = 0; i < s ; i++){
//		cout << "coord = " << x(i) << endl;
//	}
//	terminate();


	vector<Cplx> e_minus(s);
	vector<Cplx> e_infty(s);
	vector<Cplx> e_infty_derivative(s);

//	cout << "coord = " << x_coordinate << endl;
#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		e_minus[i] = E_minus(q(i),    x_coordinate, t_time) ;
		e_infty[i] = E_inf(eta, q(i), x_coordinate, t_time) ;
		e_infty_derivative[i] = E_inf_Derivative(eta, q(i), x_coordinate,
				t_time) ;
	}
//	cout << "coord and time " << x_coordinate << " " << t_time << endl;
//	cout << "kL and kR " << kL << " " << kR << endl;
//	for (int i = 0; i<s; i++){
//		cout << "(" << q(i) << ", " << weight(i) << ")" << endl;
//	}
//	cout << "E_ =" << endl;
//	for(auto elem: e_minus){
//		cout << elem << " " ;
//	}
//	cout << endl;
//	cout << "E_inf =" << endl;
//	for(auto elem: e_infty){
//		cout << elem << " " ;
//	}
//	cout << endl;
//	cout << "E_inf_der =" << endl;
//	for(auto elem: e_infty_derivative){
//		cout << elem << " " ;
//	}
//	cout << endl;

//	E_ =
//	(0.564175,0.00407237) (0.563819,0.0204357) (0.562318,0.0459185) (0.559233,0.0746248) (0.555271,0.0999219) (0.000778486,0.000163663)
//	E_inf =
//	(0.409198,0.113781) (0.41464,0.0896644) (0.420318,0.0515106) (0.422563,0.00791475) (0.420822,-0.0307861) (0.417855,-0.0555198)
//	E_inf_der =
//	(0.0264442,-0.103196) (0.0204561,-0.104705) (0.0109458,-0.106346) (3.4863e-05,-0.107145) (-0.00967775,-0.106926) (-0.0159041,-0.106306)


//#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
		//	cout << "(i,j) = (" << i <<", " << j << ") " << endl;
			Cplx v;
			Cplx w = e_infty[i] * e_minus[i] //E_+ = E_infty E_-
					* e_infty[j] * e_minus[j] * 0.5
					/ ( sin(0.5 * eta) * sin(0.5 * eta) );
			if(abs(w)< 1e-16) w = Cplx(0, 0);

			W(i, j) = sqrt(weight(i)) * w * sqrt(weight(j));
			if (i == j) {
				v = e_infty_derivative[i] * e_minus[i] * e_minus[i];
				if(abs(v)< 1e-16) v = Cplx(0, 0);
				V(i, i) = sqrt(weight(i)) * v * sqrt(weight(i));
			} else {

				v = (e_infty[i] - e_infty[j]) * e_minus[i] * e_minus[j];
				v /= q(i) - q(j);

				if(abs(v)< 1e-16) v = Cplx(0, 0);
				V(i, j) = sqrt(weight(i)) * v * sqrt(weight(j));

				V(j, i) = V(i, j);
				W(j, i) = W(i, j);
			}
			//cout << "v = " << v <<endl;
		}
	}

//	cout << "V = " << endl << V << endl;
//	cout << "W = " << endl << W << endl;
//	terminate();

	V = IdMatrix(s) + V; // Id(s,s) + V
	W = V - W;		  // Id(s,s) + V - W

	return {V, W};
}

Cplx G_inf(const double &x_coordinate, const double &t_time) {
	const double kL = LowerMomenta();
	const double kR = UpperMomenta();
	auto f = [&](double eta) {
		auto [one_plus_V, one_plus_V_minus_W] = OnePlusV_W(eta, x_coordinate,
				t_time, kL, kR);
		Cplx detV;
		Cplx detV_W;
		detV = one_plus_V.determinant();
		detV_W = one_plus_V_minus_W.determinant();
		Cplx g0 = G_0(x_coordinate, t_time) - 1.0;
	//	cout << "detV = " << detV << " detV_W = " << detV_W << " g0 = " << g0 << '\n';
		return g0 * detV + detV_W;
	};

	//const gauss<double, 20> g;
	Cplx result = g.integrate(f, -M_PI, M_PI);
//	cout << result  << endl;
//	terminate();
	return result / (2 * M_PI);
}




//MatrixXcd ConstructMatrix(const double x, const double beta, const double gamma,
//		const double magnetization, const bool finite_rank = false) {
//	const int s = 20;
//	MatrixXcd m(s, s);
//	gauss<double, 20> g;
//	auto identity = MatrixXcd::Identity(s, s);
//	bool size_parity_is_odd = g.abscissa().front() == 0 ? true : false;
//
//	auto momenta = [&](const size_t i) {
//		size_t middle_point = g.abscissa().size();
//		if (size_parity_is_odd) {
//			return i < middle_point ?
//					-M_PI * g.abscissa()[middle_point - i - 1] :
//					M_PI * g.abscissa()[i - middle_point]; //WRONG!
//		} else {
//			return i < (middle_point - 1) ?
//					-M_PI * g.abscissa()[middle_point - i - 1] :
//					M_PI * g.abscissa()[i - middle_point + 1];
//		}
//	};
//
//	auto weight = [&](const size_t i) {
//		size_t middle_point = g.weights().size();
//		if (size_parity_is_odd) {
//			return i < middle_point ?
//					M_PI * g.weights()[middle_point - i - 1] :
//										M_PI * g.weights()[i - middle_point];
//		} else {
//			return i < (middle_point - 1) ?
//					M_PI * g.weights()[middle_point - i - 1] :
//					M_PI * g.weights()[i - middle_point + 1];
//		}
//
//	};
//
//	size_t length =
//			size_parity_is_odd ?
//					2 * g.abscissa().size() : (2 * g.abscissa().size() - 1);
//
//	for (size_t i = 0; i < length; i++) {
//		for (size_t j = 0; j < length; j++) {
//			Cplx TD_weight = Weight(momenta(i), beta, gamma, magnetization);
//			Cplx kernel = Kernel(TD_weight, x, momenta(i), momenta(j));
//			if (finite_rank == true) {
//				kernel = KernelFiniteRank(kernel, TD_weight, x, momenta(i),
//						momenta(j));
//			}
//			m(i, j) = sqrt(weight(i)) * kernel * sqrt(weight(j));
//		}
//	}
//	m += identity;
//	return m;
//}
