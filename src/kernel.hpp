#pragma once

#include <complex>
#include <utility>
#include <boost/math/quadrature/gauss.hpp>


using namespace std;
using namespace boost::math::quadrature;

using Cplx = complex<double>;

const Cplx Cplx_i = Cplx(0,1);
const double MASS = 1.0;
const double B_BETA = 1000;
const double RHO = 0.5;



const gauss<double, 20> g; // Use only even size; g has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points


struct T_time{
	double value;
	explicit T_time(double new_value)
		:value{new_value}{}
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




