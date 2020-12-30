#include "matrices.hpp"
#include "profile.h"

#include "boost/math/quadrature/gauss.hpp"
#include <utility>

using namespace std;
using namespace boost::math::quadrature;

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Cplx = complex<double>; //pseudoname for complex<double>

const Cplx Cplx_i = Cplx(0,1);




MatrixXcd IdMatrix(int size){
	return MatrixXcd::Identity(size, size);
}


pair<MatrixXcd,MatrixXcd> OnePlusV_W(const double eta,const double x_coordinate, const double t_time){

	const gauss<double, 10> g;
	const int s = 10;
	MatrixXcd V(s,s);
	MatrixXcd W(s,s);

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



	for (size_t i = 0; i< s; i++){
		for (size_t j = 0; j< s; j++){
			cout << "(i,j) = (" << i << ", " << j <<")" << endl;
			{
				LOG_DURATION("V");
				V(i,j) = sqrt(weight(i))
							* V_p_q_inf(x(i), x(j), eta, x_coordinate, t_time)
							* sqrt(weight(j));
			}
			{
				LOG_DURATION("W");
				W(i,j) = sqrt(weight(i))
							* W_p_q_inf(x(i), x(j), eta, x_coordinate, t_time)
							* sqrt(weight(j));
			}
			}
	}

	return {IdMatrix(s) + V,W};
}




Cplx G_inf (const double x_coordinate, const double t_time){
	auto f = [&](double eta){
		auto [one_plus_V, W] = OnePlusV_W(eta, x_coordinate, t_time);
		Cplx detV = one_plus_V.determinant();
		Cplx detV_W = (one_plus_V - W).determinant();
		Cplx g0 = G_0(x_coordinate, t_time) - 1.0;
		return g0*detV + detV_W;
	};


	const gauss<double, 50> g;
	Cplx result = g.integrate(f, -M_PI,M_PI);
	return result/(2*M_PI);
}





MatrixXcd ConstructMatrix(
							const double x,
							const double beta,
							const double gamma,
							const double magnetization,
							const bool finite_rank = false){
	const int s = 10;
	MatrixXcd m(s,s);
	gauss<double, 10> g;
	auto identity = MatrixXcd::Identity(s, s);
	bool size_parity_is_odd = g.abscissa().front() == 0 ? true : false ;

	auto momenta = [&](const size_t i){
		size_t middle_point = g.abscissa().size();
		if (size_parity_is_odd){
			return  i < middle_point ? -M_PI*g.abscissa()[middle_point - i - 1]
									  : M_PI*g.abscissa()[i - middle_point]; //WRONG!
		} else {
			return  i < (middle_point-1) ? -M_PI*g.abscissa()[middle_point - i - 1]
										  : M_PI*g.abscissa()[i - middle_point + 1];
		}
	};

	auto weight = [&](const size_t i){
		size_t middle_point = g.weights().size();
		if (size_parity_is_odd){
			return  i < middle_point ? M_PI*g.weights()[middle_point - i - 1]
									 : M_PI*g.weights()[i - middle_point];
		} else {
			return  i < (middle_point-1) ? M_PI*g.weights()[middle_point - i - 1]
										 : M_PI*g.weights()[i - middle_point + 1];
		}

	};

	size_t length = size_parity_is_odd ? 2 * g.abscissa().size()
			: (2 * g.abscissa().size() -1);

	for (size_t i = 0; i<length; i++){
		for (size_t j = 0; j< length; j++){
				Cplx TD_weight  = Weight (momenta(i), beta, gamma, magnetization);
				Cplx kernel = Kernel (TD_weight, x, momenta(i), momenta(j) );
				if (finite_rank == true){
					kernel = KernelFiniteRank (kernel , TD_weight, x, momenta(i), momenta(j));
				}
				m(i,j) = sqrt(weight(i)) * kernel * sqrt(weight(j));
			}
	}
	m += identity;
	return  m;
}
