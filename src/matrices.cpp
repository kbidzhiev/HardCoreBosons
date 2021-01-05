#include "matrices.hpp"
#include "profile.h"

#include "boost/math/quadrature/gauss.hpp"
#include <utility>
#include <omp.h>

using namespace std;
using namespace boost::math::quadrature;

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Cplx = complex<double>;
//pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0, 1);

MatrixXcd IdMatrix(int size) {
	return MatrixXcd::Identity(size, size);
}

pair<MatrixXcd, MatrixXcd> OnePlusV_W(const double &eta,
		const double &x_coordinate, const double &t_time) {

	const gauss<double, 8> g;
	const int s = 8;
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

//omp_set_num_threads (omp_get_num_procs()); // number of threads = num of processors
//#pragma omp parallel for collapse(2)

#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		for (size_t j = i; j < s; j++) {
			//cout << "V(i,j) = (" << i << ", " << j << ") \n";

			{
				//LOG_DURATION("V");
				if (i == j) {
					Cplx v = V_diag_inf(x(i), eta, x_coordinate, t_time);
					V(i, j) = sqrt(weight(i)) * v * sqrt(weight(j));
				} else {
					Cplx v = V_p_q_inf(x(i), x(j), eta, x_coordinate, t_time);

					V(i, j) = sqrt(weight(i)) * v * sqrt(weight(j));
					V(j, i) = V(i, j);
					//cout << "V(i,j) = " << V(i,j) << "\n";
				}
			}
			{
				//LOG_DURATION("W");
				if (i == j) {
					W(i, j) = 0;
				} else {

					W(i, j) = sqrt(weight(i))
							* W_p_q_inf(x(i), x(j), eta, x_coordinate, t_time)
							* sqrt(weight(j));
				}

			}
		}
	}
	return {IdMatrix(s) + V, IdMatrix(s) + V - W};
}

Cplx G_inf(const double &x_coordinate, const double &t_time) {
	auto f = [&](double eta) {
		auto [one_plus_V, W] = OnePlusV_W(eta, x_coordinate, t_time);
		Cplx detV = one_plus_V.determinant();
		Cplx detV_W = W.determinant();
		Cplx g0 = G_0(x_coordinate, t_time) - 1.0;
		//cout << "detV = " << detV << " detV_W = " << detV_W << " g0 = " << g0 << '\n';
		return g0 * detV + detV_W;
	};

	const gauss<double, 12> g;
	Cplx result = g.integrate(f, -M_PI, M_PI);
	return result / (2 * M_PI);
}

MatrixXcd ConstructMatrix(const double x, const double beta, const double gamma,
		const double magnetization, const bool finite_rank = false) {
	const int s = 10;
	MatrixXcd m(s, s);
	gauss<double, 10> g;
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
