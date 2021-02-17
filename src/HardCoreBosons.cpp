//============================================================================
// Name        : HardCoreBosons_new.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>

#include "kernel.hpp"
#include "profile.h"
#include "test_runner.h"



using namespace std;

void LambdaCurve(){

	ofstream kernel; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "kernel" + to_string(2 * g.weights().size()) + ".dat";
	kernel.open(filename, mode);
	kernel.precision(15);

	SpaceTime st(X_coordinate(1.0),T_time(0.1));
		for (double Lambda = -10; Lambda <= 10; Lambda += 0.1){

			auto [detv,detw] = Determinants(Lambda, st);
			Cplx coeff = G0 (st) * 2.0/(Lambda * Lambda + 1.0) ;
			detv *= (coeff - 1.0);
			kernel << Lambda << "\t" << real(detv + detw) << "\t" << imag(detv + detw)
								 << endl;
		}
	auto f = [&](double Lambda) {
		auto [detv, detw] = Determinants(Lambda, st);
		Cplx coeff = G0(st) * 2.0 / (Lambda * Lambda + 1.0);
		detv *= (coeff - 1.0);
		return detv + detw;
	};
	double error = 100;
	cout << gauss_kronrod<double, 31>::integrate(f, -50, 50,  10, 1e-9, &error) << endl;
	cout << "error =" <<  error << endl;

//	(-2.4779,-3.15332)
//	error =0.0180585

	// 5 smaller interval
//	(-2.13641,-2.7455)
//	error =5.79921e-10
}


void CorrelatorCurve(){
	ofstream correlator; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "Correlator" + to_string(2 * g.weights().size()) + ".dat";
	correlator.open(filename, mode);
	correlator.precision(15);

	const double X_LIMITS = 1.0 * KF();
	const double T_LIMITS = 1.0  * Energy(Q_momenta(KF()));

	for (double time = 0.1*T_LIMITS; time < T_LIMITS; time += 0.01) {
		correlator << "\"t=" << time << "\"\n" ;

		for (double coordinate = X_LIMITS; coordinate <= X_LIMITS; coordinate += 0.1) {

			//SpaceTime st(X_coordinate(coordinate), T_time(time));
			Cplx result = Grep({X_coordinate(coordinate), T_time(time)});

			correlator << coordinate << "\t" << real(result) << "\t" << imag(result)
					<< "\t" << time << endl;
		}
		correlator << "\n\n" ;
	}
}


int main() {
	LOG_DURATION("DET");






	LambdaCurve();

	//CorrelatorCurve();




	return 0;
}
