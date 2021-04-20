#pragma once

#include <complex>
#include <utility>
#include <stdexcept>
#include <string>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "Faddeeva.hh"
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"

using namespace std;
using namespace boost::math::quadrature;
using boost::math::quadrature::trapezoidal;

using Cplx = complex<double>;

const Cplx Cplx_i = Cplx(0,1);
const double MASS = 1.0;
const double B_BETA = 1000;
const double RHO = 0.5;


const size_t GAUSS_RANK = 40;

const gauss<double, GAUSS_RANK> g; // Use only even size; g has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points



struct T_time {
	double value;
	explicit T_time(double new_value) {
//		if (new_value < 0) {
//			throw logic_error(
//					"Time can not be negative. t = " + to_string(new_value));
//		}
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
//Cplx Erf(const Cplx z);
Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime st);
Cplx PrincipalValue_old(Q_momenta q_momenta,  SpaceTime st);
Cplx G0 (SpaceTime st);
pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime st);
Cplx Grep(SpaceTime st);

void Fourier2D();
void Fourier1D();
void Gpt();
void foo();
