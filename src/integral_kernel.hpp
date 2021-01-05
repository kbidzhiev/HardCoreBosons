#pragma once

#include <cmath>
#include <complex>

using namespace std;

using Cplx = complex<double>; //pseudoname for complex<double>



double Energy (const double& q_momenta);

double Tau (const double& q_momenta,
		const double& x_coordinate,
		const double& t_time
);

Cplx PrincipalValue(const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);
// from (-1,1) I integrate with Gauss, from 1, to \infty with trapezoid
// According to "Gauss-Legendre Principal value Integration"
// by Julian V. Noble
// DOI:10.1109/MCISE.2000.970778

Cplx PrincipalValueDerivative(const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);


Cplx E_inf(const double& eta,
		const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);

Cplx E_inf_Derivative(const double& eta,
		const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);

double Teta(const double b_beta,
		const double energy,
		const double chem_potential);

Cplx E_minus(const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);

Cplx E_plus(const double& eta,
		const double& q_momenta,
		const double& x_coordinate,
		const double& t_time);

Cplx V_p_q_inf(const double& p_momenta,
		const double& q_momenta,
		const double& eta,
		const double& x_coordinate,
		const double& t_time);

Cplx V_diag_inf (const double& p_momenta,
		const double& eta,
		const double& x_coordinate,
		const double& t_time);


Cplx W_p_q_inf(const double& p_momenta,
		const double& q_momenta,
		const double& eta,
		const double& x_coordinate,
		const double& t_time);


Cplx G_0(const double& x_coordinate, const double& t_time);//Eq (3.22) arxiv.org/abs/1511.05922


Cplx Weight (const double momenta,
		const double b_beta,
		const double gamma,
		const double magnetization);

Cplx Kernel (const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2);

Cplx KernelFiniteRank (const Cplx kernel,
		const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2);



