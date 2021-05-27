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
//#include "kernel_lattice.hpp"
#include "profile.h"
#include "test_runner.h"


using namespace std;

void LambdaCurve(){

	ofstream kernel; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "kernel" + to_string(GAUSS_RANK) + ".dat";
	kernel.open(filename, mode);
	kernel.precision(15);

	SpaceTime st(X_coordinate(-15.5),T_time(0.01));
		for (double Lambda = -10; Lambda <= 10; Lambda += 0.01){

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
	Cplx gauss_res = gauss_kronrod<double, 61>::integrate(f, -50.0, 50.0, 15, 1e-9, &error) ;
	//Cplx trap_res  = trapezoidal(f, -50.0, 50.0);
	cout << "gauss = " << gauss_res << "\n";
	//		"rap = " << trap_res << "dif(g-t) = " << abs(gauss_res - trap_res) << endl;
	//(0.0503884,0.309736)

}

void CorrelatorCurve(){
	ofstream correlator; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "Correlator" + to_string(GAUSS_RANK) + ".dat";
	correlator.open(filename, mode);
	correlator.precision(15);

	const double X_LIMITS = 5.0 ;// * KF();
	const double T_LIMITS = 15.0 ;// * Energy(Q_momenta(KF()));

	for (double time = 0.01*T_LIMITS; time < T_LIMITS; time += 0.1) {
		correlator << "\"t=" << time << "\"\n" ;

		for (double coordinate = -X_LIMITS; coordinate <= X_LIMITS; coordinate += 0.1) {
			cout << "x = " << coordinate << " t = " << time << endl;
			//SpaceTime st(X_coordinate(coordinate), T_time(time));
			Cplx result = Grep({X_coordinate(coordinate), T_time(time)});

			correlator << coordinate << "\t" << real(result) << "\t" << imag(result)
					<< "\t" << time << endl;
		}
		correlator << "\n\n" ;
	}
}

void TsliceCurve(double time_){
	ofstream tslice; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "tslice_" + to_string((int)time_)
			+ "_" + to_string(GAUSS_RANK) + ".dat";
	tslice.open("Data/"+filename, mode);
	tslice.precision(15);

	const double X_LIMITS = 10.0 ;
	const double T_LIMITS = time_;  // Energy(Q_momenta(KF()));

	for (double coordinate = -X_LIMITS; coordinate <= X_LIMITS; coordinate += 0.1) {
		cout << coordinate << " / " << X_LIMITS  << endl;
		Cplx result = Grep( { X_coordinate(coordinate), T_time(T_LIMITS) });
		tslice << coordinate << "\t" << real(result) << "\t" << imag(result) << endl;
	}
}

void XsliceCurve(double x){
	ofstream xslice; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "xslice_" + to_string((int)x) + ".dat";
	xslice.open(filename, mode);
	xslice.precision(15);

	const double X_LIMITS = x ;
	const double T_LIMITS = 20.0 ;

	for (double time = 0.25; time < T_LIMITS; time += 0.01) {
		//xslice << "\"t=" << time << "\"\n";
		Cplx result = Grep( { X_coordinate(X_LIMITS), T_time(time) });
		xslice << time << "\t" << real(result) << "\t" << imag(result) << endl;
		cout << "time = " << time << " /" << T_LIMITS << endl;
	}
}

void PV(){
	Q_momenta q(10.5);
	SpaceTime st(X_coordinate(1.7),T_time(10.6));
	Cplx old_value = PrincipalValue_old(q, st)/M_PI;
	Cplx new_value = PrincipalValue(q, st)/M_PI;
	Cplx erf_arg = (st.x - q.value * st.t) * (-1.0 + Cplx_i) / (2.0 * sqrt(st.t));

	cout << old_value << '\n' << new_value << endl;
	cout << "erf          = " << Faddeeva::erf(erf_arg) << '\n';
	cout << "erf(wolfram) = (1.00787,-0.0223579)" << endl;

}


int main() {
	LOG_DURATION("Total");

	//LambdaCurve();
	//CorrelatorCurve();

	//XsliceCurve(0.0); // x = 1;
	//TsliceCurve(5.0);//5 second
	//PV();

	//Fourier1D();
	//Fourier2D();
	//Gpt();
	//foo();
	Q_momenta q(3.14);
	SpaceTime spacetime(X_coordinate(50.0), T_time(30.0));

	cout << EPV_l (q, spacetime) << endl;

	return 0;
}
