#include "matrices.hpp"
#include "profile.h"

#include "gauss.hpp"
#include <utility>
#include <omp.h>
#include "../../cpp_libs/eigen/Eigen/Dense"

#include <vector>

#include <exception>

using namespace std;
using namespace boost::math::quadrature;


using Eigen::MatrixXcd;

using Cplx = complex<double>;
//pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0, 1);

gauss<double, 10> g;

MatrixXcd IdMatrix(int size) {
	return MatrixXcd::Identity(size, size);
}


pair<MatrixXcd, MatrixXcd> OnePlusV_W(const double eta,
		const double x_coordinate, const double t_time,
		const double  momenta_truncation) {

	const int s = 10;
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

	const auto q = [&](const size_t i) {
		size_t middle_point = g.abscissa().size();
		return i < middle_point ?
				- momenta_truncation * g.abscissa()[middle_point - i - 1] :
				 momenta_truncation * g.abscissa()[i - middle_point];
	};

	const auto weight = [&](const size_t i) {
		size_t middle_point = g.weights().size();
		return i < middle_point ?
				momenta_truncation * g.weights()[middle_point - i - 1] :
				momenta_truncation * g.weights()[i - middle_point];
	};

//	for (int i = 0; i < s ; i++){
//		cout << "coord = " << q(i) << endl;
//	}
//	terminate();


	vector<Cplx> e_minus(s);
	vector<Cplx> e_infty(s);
	vector<Cplx> e_infty_derivative(s);

	const Cplx fresnel = exp(-Cplx_i * 3. * M_PI / 4.) * sqrt(2 * M_PI * t_time)
			* exp(Cplx_i * x_coordinate * x_coordinate / (2 * t_time));


	cout << "fresn = " << fresnel << endl;
	terminate();

	#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		e_minus[i] = E_minus(q(i),    x_coordinate, t_time) ;
		e_infty[i] = E_inf(eta, q(i), x_coordinate, t_time) ;
	}


#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		e_infty_derivative[i] = fresnel - Cplx_i * (q(i) * t_time - x_coordinate) * e_infty[i];
	}

//	cout << "coord and time " << x_coordinate << " " << t_time << endl;
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
//
//	terminate();

#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
		//	cout << "(i,j) = (" << i <<", " << j << ") " << endl;
			Cplx v;
			Cplx w = e_infty[i] * e_minus[i] //E_+ = E_infty E_-
					* e_infty[j] * e_minus[j] * 0.5
					/ ( sin(0.5 * eta) * sin(0.5 * eta) );
			if(abs(w)< 1e-15) w = Cplx(0, 0);
			W(i, j) = sqrt(weight(i)) * w * sqrt(weight(j));
			if (i == j) {
				v = e_infty_derivative[i] * e_minus[i] * e_minus[i];
				if(abs(v)< 1e-15) v = Cplx(0, 0);
				V(i, i) = sqrt(weight(i)) * v * sqrt(weight(i)) ;
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

	V = IdMatrix(s) + V;
	W = V - W;		  // Id(s,s) + V - W

//	cout << "V = " << endl << V << endl;
//	cout << "W = " << endl << W << endl;
//
//	cout << "det V" << V.determinant() << endl;
//	cout << "det W" << W.determinant() << endl;
//	terminate();

	return {V, W};
}

Cplx G_inf(const double x_coordinate, const double t_time) {
	const double momenta_truncation = UpperMomenta();
	auto f = [&](double eta) {
		auto [one_plus_V, one_plus_V_minus_W] = OnePlusV_W(eta, x_coordinate,
				t_time, momenta_truncation);
		// determinants can be computer parallel. Does it worth ?
		Cplx detV = one_plus_V.determinant();
		Cplx detV_W = one_plus_V_minus_W.determinant();
		Cplx g0 = G_0(x_coordinate, t_time) - 1.0;
		//cout << "detV = " << detV << " detV_W = " << detV_W << " g0 = " << g0 << '\n';
		return g0 * detV + detV_W;
	};

	const gauss<double, 6> g2;
	Cplx result = g2.integrate(f, -M_PI, M_PI);
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
