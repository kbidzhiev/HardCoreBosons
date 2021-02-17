#pragma once

#include <complex>
#include <utility>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

//#include "~/Programs/cpp_libs/boost/boost/math/quadrature/gauss.hpp"
//#include "~/Programs/cpp_libs/boost/boost/math/quadrature/trapezoidal.hpp"
//#include "~/Programs/cpp_libs/boost/boost/math/quadrature/gauss_kronrod.hpp"


#include "../../cpp_libs/eigen/Eigen/Dense"
#include <stdexcept>
#include <string>

using namespace std;
using namespace boost::math::quadrature;

using Cplx = complex<double>;

const Cplx Cplx_i = Cplx(0,1);
const double MASS = 1.0;
const double B_BETA = 1000;
const double RHO = 0.5;



const gauss<double, 20> g; // Use only even size; g has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points
/*
 * g  			integr
 *
 * 6			-2.42971,-3.09114
 * 20			-2.42865,-3.09458
 * 30			-2.42847,-3.09571
 * 60			-2.4291,-3.09473
 * */


// x = 10.0, t = 0.1
// 100 0.0958141,-0.0294548
// 20  0.0958088,-0.0294753

struct T_time {
	double value;
	explicit T_time(double new_value) {
		if (new_value < 0) {
			throw logic_error(
					"Time can not be negative. t = " + to_string(new_value));
		}
		value = new_value;
	}

};

struct X_coordinate{
	double value;
	explicit X_coordinate(double new_value)
		:value{new_value}{}
};

struct Q_momenta{
	double value;
	explicit Q_momenta(double new_value)
		:value{new_value}{}
};

struct SpaceTime{
	double x;
	double t;
	SpaceTime(X_coordinate x_coord, T_time t_time)
		:x(x_coord.value)
		,t(t_time.value){}
};


double KF();
double Mu_chempot();
double Energy(Q_momenta q_momenta);
double Tau(Q_momenta q_momenta,  SpaceTime st);
double Theta (Q_momenta q_momenta);
Cplx Eminus (Q_momenta q_momenta,  SpaceTime st);
Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime st);
Cplx G0 (SpaceTime st);
pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime st);
Cplx Grep(SpaceTime st);




