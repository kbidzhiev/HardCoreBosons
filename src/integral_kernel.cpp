
#include "integral_kernel.hpp"

#include <complex.h> // trigoniometric functions sin() cos()


using namespace std;

const Cplx Cplx_i = Cplx(0,1);

const double mass =  1.;

const double Energy (const double q_momenta){
	return  q_momenta * q_momenta * 0.5 / mass;
}

const double Tau (const double q_momenta,
		const double coordinate,
		const double time
){
	return time * Energy(q_momenta) - coordinate * q_momenta;
}








Cplx G_0(const double coordinate, const double time){
	Cplx result = exp(-Cplx_i * M_PI /4.);
	result *= sqrt(1./(2 * M_PI * time));
	result *= exp( (Cplx_i * coordinate * coordinate )/ (2 * time) );
	return result;
}

Cplx Weight (const double momenta,
		const double beta,
		const double gamma,
		const double magnetization){

	const Cplx sqrt_term = sqrt(pow(magnetization-cos(momenta),2) + pow(gamma * sin(momenta),2));

	Cplx result =
			magnetization - cos(momenta) - Cplx_i * gamma * sin(momenta) ;

	result /= sqrt_term ;
	result *= -tanh(0.5*beta * sqrt_term);
	result += 1;
	result *= 0.5;
	return result;
}


Cplx Kernel (const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2) {
	Cplx kernel = sin(0.5 * x * (momenta1-momenta2))/
			sin(0.5 * (momenta1-momenta2));
	Cplx result = -weight/M_PI;
	if (abs (momenta1 - momenta2) > 1e-14 ){
		result *= kernel;
	} else {
		result *= x;
	}
	return  result;
}

Cplx KernelFiniteRank (const Cplx kernel,
		const Cplx weight,
		const double x,
		const double momenta1,
		const double momenta2
		){
	Cplx result = kernel;
	result += weight * exp(-Cplx_i * x * 0.5* (momenta1 + momenta2))
		* exp(Cplx_i * 0.5* (momenta1 - momenta2))/M_PI;
	return result;
}







const Cplx E_plus();
const Cplx E_minus();

const Cplx V_p_q(const double p_momenta, const double q_momenta, const double g_coupling){
	return (E_plus(p_momenta) * E_minus(q_momenta) - E_minus(p_momenta) * E_plus (q_momenta))/(p_momenta- q_momenta);
}

const Cplx V_inf_p_q(const double p_momenta, const double q_momenta){
	return V_p_q(p_momenta, q_momenta, 999);
}


