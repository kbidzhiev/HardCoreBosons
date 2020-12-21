#pragma once

#include <cmath>
#include <complex>

using namespace std;

using Cplx = complex<double>; //pseudoname for complex<double>



Cplx Weight (const double momenta,
		const double beta,
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

