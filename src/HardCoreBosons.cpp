//============================================================================
// Name        : HardCoreBosons.cpp
// Author      : BK
// Description : HardCoreBosons fredholm determinant
//============================================================================

#include "test_runner.h"
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



void TestPVIntegration(){
	LOG_DURATION("time");
	const double momenta = 1.2;
	const double coordinate = -2.1;
	const double time = 0.2;
	const complex<double> PV = PrincipalValue(momenta, coordinate, time);
	const complex<double> PV_deriv = PrincipalValueDerivative(momenta, coordinate, time);
	const complex<double> PV_wolfram = -1.10414 + 3.12551 * complex<double>(0,1);
	const complex<double> PV_deriv_wolfram = 6.49805 + 3.35269 * complex<double>(0,1);
	cout << PV << '\n' << PV_deriv << endl ;
	ASSERT(abs(PV - PV_wolfram) < 1e-3);
	ASSERT_EQUAL(abs(PV_deriv - PV_deriv_wolfram) < 1e-3, true);
}

int main(){
//	TestRunner tr;
//	RUN_TEST(tr, TestPVIntegration);


	 LOG_DURATION("Total");

		ofstream correlator; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)
		correlator.open("Correlator.dat", mode);
		correlator.precision(15);
		correlator << "#x \t correlator \t time \n";

		// ------- Correlator profile -------

		//int(param.val("Sz");


		const double dt = 0.1;
		const double system_size = 0.5;
		const double dx = 0.01;
		const int n_steps = system_size / dx;
		const double time_total = 0.1;
		const int T =  time_total / dt;
		const int total_steps = (n_steps + 1) * T ;
		cout << n_steps << " " << T << endl;
		int counter = 0;
		for (double time = 0.001; time < time_total; time += dt) {
			LOG_DURATION("time step");
			correlator << "\"t=" << time << "\"" << endl;
			for (int n = -n_steps/2; n <= n_steps / 2; ++n) {
				const double x_coordinate = n * dx; //+param.val("time_shift");
				const complex<double> g_inf = G_inf(x_coordinate, time);
				correlator << x_coordinate << "\t" << real(g_inf) << "\t"
				<< imag(g_inf) << "\n";
				cout << "step " << counter << "/" << total_steps << '\n';
				counter++;
			}
			correlator << "\n \n";
	}
	cout << "DONE !" << endl;

}
