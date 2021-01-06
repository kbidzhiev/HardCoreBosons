//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================


#include "profile.h"

#include "matrices.hpp"

#include <fstream>	 // file input output
#include <iostream>	 // screen input output
#include <iomanip>   // setprecision(15)
#include <cmath> 	 // pow (x,3) = x^3 = x*x*x and M_PI = pi = 3.14
					 // trigoniometric functions sin() cos()
#include <complex.h>


//
//using Cplx = complex<double>; //pseudoname for complex<double>
//
//const Cplx Cplx_i = Cplx(0,1);



/*
 *  gauss<double, 10> g; in wolfram its
 *  GaussianQuadratureWeights[10, -1, 1]
 * 	g.weights() = {0.295524, 0.269267, 0.219086, 0.149451, 0.0666713} always positive
 * 	g.abscissa()= {0.148874, 0.433395, 0.67941,  0.865063, 0.973907 } positive part
 *
 *	to obtain GaussianQuadratureWeights[10, -pi, pi] one should
 *	pi*g.weights() and pi*g.abscissa()
 */



struct Parameters {
	const double mass = 1.;
	const double g_coupling = 999.;
	const double b_beta = 10.;
	const double chem_potential = 0;
};

int main(){



	{

		ofstream correlator; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)
		correlator.open("Correlator.dat", mode);
		correlator.precision(15);
		correlator << "#x \t correlator \t time \n";

		// ------- Correlator profile -------
		double dx = 0.1;
		double system_size = 1.0;
		const double t_time = 0.1;
		const int n_steps = system_size / dx;
		for (int n = -n_steps / 2; n <= n_steps / 2; ++n) {
			LOG_DURATION("time step");
			const double x_coordinate = n * dx; //+param.val("time_shift");
			const complex<double> determ = G_inf (x_coordinate, t_time);
			correlator << x_coordinate << "\t" << real(determ) << "\t" << imag(determ) << "\n" ;
		}
	}
//	double pole = 10.0;
//	double x_coordinate = -3.0;
//	double t_time = 2.0;
//	{	LOG_DURATION("total");
//		//cout << setprecision(15) << endl;
//		//cout << PrincipalValue ( pole,  x_coordinate,  t_time) << endl;
//		//cout << PrincipalValueDerivative( pole,  x_coordinate,  t_time) << endl;
//	//	cout << G_inf (x_coordinate, t_time) << endl;
//	}
}
