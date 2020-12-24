#pragma once

#include <cmath>
#include <complex>

using namespace std;

using Cplx = complex<double>; //pseudoname for complex<double>



const double Energy (const double q_momenta);

const double Tau (const double q_momenta,
		const double x_coordinate,
		const double t_time
);

Cplx PrincipalValue(const double q_momenta,
		const double x_coordinate,
		const double t_time); // According to "Gauss-Legendre Principal value Integration"
// by Julian V. Noble
// DOI:10.1109/MCISE.2000.970778



const Cplx E_inf(const double eta,
		const double q_momenta,
		const double x_coordinate,
		const double t_time);

const double Teta(const double b_beta,
		const double energy,
		const double chem_potential);

const Cplx E_minus(const double q_momenta,
		const double x_coordinate,
		const double t_time);

const Cplx E_plus(const double eta,
		const double q_momenta,
		const double x_coordinate,
		const double t_time);

const Cplx V_p_q_inf(const double p_momenta,
		const double q_momenta,
		const double eta,
		const double x_coordinate,
		const double t_time);

const Cplx W_p_q_inf(const double p_momenta,
		const double q_momenta,
		const double eta,
		const double x_coordinate,
		const double t_time);



Cplx G_0(const double x_coordinate, const double t_time);//Eq (3.22) arxiv.org/abs/1511.05922


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



