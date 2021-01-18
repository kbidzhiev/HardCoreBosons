#include "matrices.hpp"
#include "profile.h"

#include "boost/math/quadrature/gauss.hpp"
#include <utility>
#include <omp.h>

#include <vector>

using namespace std;
using namespace boost::math::quadrature;

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Cplx = complex<double>;
//pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0, 1);

gauss<double, 20> g;

MatrixXcd IdMatrix(int size) {
	return MatrixXcd::Identity(size, size);
}

pair<MatrixXcd, MatrixXcd> OnePlusV_W(const double &eta,
		const double &x_coordinate, const double &t_time) {

	const int s = 20;
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);

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

	vector<Cplx> e_minus, e_infty, e_infty_derivative;
	e_minus.reserve(s);
	e_infty.reserve(s);
	e_infty_derivative.reserve(s);

#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		e_minus[i] = E_minus(x(i), x_coordinate, t_time);
		e_infty[i] = E_inf(eta, x(i), x_coordinate, t_time);
		e_infty_derivative[i] = E_inf_Derivative(eta, x(i), x_coordinate,
				t_time);
	}

#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			Cplx v;
			Cplx w = e_infty[i] * e_minus[i] //E_+ = E_infty E_-
					* e_infty[j] * e_minus[j] * 0.5
					/ ( sin(0.5 * eta) * sin(0.5 * eta) );
			W(i, j) = sqrt(weight(i)) * w * sqrt(weight(j));
			if (i == j) {
				v = e_infty_derivative[i] * e_minus[i] * e_minus[i];
				V(i, j) = sqrt(weight(i)) * v * sqrt(weight(i));
			} else {

				v = (e_infty[i] - e_infty[j]) * e_minus[i] * e_minus[j];
				v /= x(i) - x(j);
				V(i, j) = sqrt(weight(i)) * v * sqrt(weight(j));

				V(j, i) = V(i, j);
				W(j, i) = W(i, j);
			}
		}
	}
	V = IdMatrix(s) + V; // Id(s,s) + V
	W = V - W;		  // Id(s,s) + V - W

	return {V, W};
}

Cplx G_inf(const double &x_coordinate, const double &t_time) {
	auto f = [&](double eta) {
		auto [one_plus_V, one_plus_V_minus_W] = OnePlusV_W(eta, x_coordinate,
				t_time);
		Cplx detV;
		Cplx detV_W;
		detV = one_plus_V.determinant();
		detV_W = one_plus_V_minus_W.determinant();
		Cplx g0 = G_0(x_coordinate, t_time) - 1.0;
		//cout << "detV = " << detV << " detV_W = " << detV_W << " g0 = " << g0 << '\n';
		return g0 * detV + detV_W;
	};

	//const gauss<double, 20> g;
	Cplx result = g.integrate(f, -M_PI, M_PI);
	return result / (2 * M_PI);
}

MatrixXcd ConstructMatrix(const double x, const double beta, const double gamma,
		const double magnetization, const bool finite_rank = false) {
	const int s = 20;
	MatrixXcd m(s, s);
	gauss<double, 20> g;
	auto identity = MatrixXcd::Identity(s, s);
	bool size_parity_is_odd = g.abscissa().front() == 0 ? true : false;

	auto momenta = [&](const size_t i) {
		size_t middle_point = g.abscissa().size();
		if (size_parity_is_odd) {
			return i < middle_point ?
					-M_PI * g.abscissa()[middle_point - i - 1] :
					M_PI * g.abscissa()[i - middle_point]; //WRONG!
		} else {
			return i < (middle_point - 1) ?
					-M_PI * g.abscissa()[middle_point - i - 1] :
					M_PI * g.abscissa()[i - middle_point + 1];
		}
	};

	auto weight = [&](const size_t i) {
		size_t middle_point = g.weights().size();
		if (size_parity_is_odd) {
			return i < middle_point ?
			M_PI * g.weights()[middle_point - i - 1] :
										M_PI * g.weights()[i - middle_point];
		} else {
			return i < (middle_point - 1) ?
					M_PI * g.weights()[middle_point - i - 1] :
					M_PI * g.weights()[i - middle_point + 1];
		}

	};

	size_t length =
			size_parity_is_odd ?
					2 * g.abscissa().size() : (2 * g.abscissa().size() - 1);

	for (size_t i = 0; i < length; i++) {
		for (size_t j = 0; j < length; j++) {
			Cplx TD_weight = Weight(momenta(i), beta, gamma, magnetization);
			Cplx kernel = Kernel(TD_weight, x, momenta(i), momenta(j));
			if (finite_rank == true) {
				kernel = KernelFiniteRank(kernel, TD_weight, x, momenta(i),
						momenta(j));
			}
			m(i, j) = sqrt(weight(i)) * kernel * sqrt(weight(j));
		}
	}
	m += identity;
	return m;
}
