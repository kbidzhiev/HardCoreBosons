#pragma once

#include <complex>
#include <utility>
//#include <stdexcept>
#include <string>
#include "Faddeeva.hh"
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include "../../cpp_libs/eigen/Eigen/Dense"
#include "profile.h"

using namespace std;
using namespace boost::math::quadrature;

using Cplx = complex<double>;

const Cplx Cplx_i = Cplx(0,1);
const double MASS = 1.0;
const double B_BETA = 100.0;
const double RHO = 0.5;
const double MAGN_FIELD = 0.;

// product of MAGN_FIELD * B_BETA should be less than 700
// for avoiding an overflow



const size_t GAUSS_RANK = 120;		// use only ODD = 2n+1 numbers; Precomputed 15, 31, 41, 51 and 61
const size_t GAUSS_RANK_l = 10; 	 // Use only EVEN rank
/*
 * increasing of GAUSS_RANK_l improves GrepEta_l(eta, spacetime).
 * for low GAUSS_RANK_l profiles for large (x,t) looks periodic
 *
 * */
const size_t GAUSS_INTEGRATION = 31; // Use only ODD 2n+1 rank // Precomputed 15, 31, 41, 51 and 61

const gauss<double, GAUSS_RANK> g; // Use only even size; g has pre-computed tables of abscissa and weights for 7, 15, 20, 25 and 30 points


//const gauss_kronrod<double, GAUSS_RANK> g; //  Precomputed 15, 31, 41, 51 and 61


//const gauss<double, GAUSS_RANK> g; // Precomputed //7, 15, 20, 25 and 30




const gauss<double, GAUSS_RANK_l> g_l; //7, 15, 20, 25 and 30
const gauss_kronrod<double, GAUSS_INTEGRATION> g_integration; // Precomputed 15, 31, 41, 51 and 61


struct T_time {
	double value;
	explicit T_time(double new_value)
		:value{new_value}{
//		if (new_value < 0) {
//			throw logic_error(
//					"Time can not be negative. t = " + to_string(new_value));
//		}
	}
};

struct X_coordinate{
	complex<double> value;
	explicit X_coordinate(complex<double> new_value)
		:value{new_value}{}
};

struct Q_momenta{
	double value;
	explicit Q_momenta(double new_value)
		:value{new_value}{}
};

struct SpaceTime{
	complex<double> x;
	double t;
	SpaceTime(X_coordinate x_coord, T_time t_time)
		:x{x_coord.value}
		,t{t_time.value}{}
};



/////////// Continuum version ///////////
double KF();
double Mu_chempot();
double Energy(Q_momenta q_momenta);
complex<double> Tau(Q_momenta q_momenta,  SpaceTime st);
double Theta (Q_momenta q_momenta);
Cplx Eminus (Q_momenta q_momenta,  SpaceTime st);
Cplx Einf (double eta, Q_momenta q_momenta,  SpaceTime spacetime);
Cplx Eplus (double eta, Q_momenta q_momenta,  SpaceTime spacetime);
Cplx PrincipalValue(Q_momenta q_momenta,  SpaceTime st);
Cplx PrincipalValue_old(Q_momenta q_momenta,  SpaceTime st);
Cplx G0 (SpaceTime st);
pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime st);
Cplx GrepLambda(double Lambda,  SpaceTime st);
Cplx Grep(SpaceTime st);

pair <Cplx, Cplx> Determinants_eta(double eta, SpaceTime st);
Cplx GrepEta(double eta,  SpaceTime st);
Cplx Grep_new(SpaceTime st);

/////////// Lattice version ///////////
double Energy_l(Q_momenta q_momenta);
complex<double> Tau_l(Q_momenta q_momenta,  SpaceTime st);
double Theta_l(Q_momenta q_momenta);
Cplx Eminus_l (Q_momenta q_momenta,  SpaceTime st);
Cplx Eplus_l (Q_momenta q_momenta,  SpaceTime st);
Cplx EPV_l (Q_momenta q_momenta,  SpaceTime st);
Cplx EPV_derivative_l (Q_momenta q_momenta,  SpaceTime st);
Cplx Lplus_l (double eta, Q_momenta q_momenta,  SpaceTime st);
Cplx Lminus_l (Q_momenta q_momenta,  SpaceTime st);
Cplx Nu_matrix_elem(double eta, Q_momenta k_momenta, Q_momenta k_prime_momenta, SpaceTime spacetime);
Cplx Nu_diagonal_l (double eta, Q_momenta q_momenta, SpaceTime st);
Cplx Rplus_l (double eta, Q_momenta q_momenta, Q_momenta k_momenta,  SpaceTime st);
Cplx G_l (SpaceTime st);
double F_l(double eta);
double Gamma_l();
pair <Cplx, Cplx> Determinants_l(double eta, SpaceTime st);
Cplx GrepEta_l(double eta,  SpaceTime st);
Cplx Grep_l(SpaceTime st);


Cplx Asymptotics(SpaceTime st);
void Fourier2D();
void Fourier1D();
void Gpt();
void Gp0t();
void foo();



///////////// UNIT TESTS ///////////

void UnitTests();
