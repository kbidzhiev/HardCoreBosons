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



void CorrelatorCurve(){
	ofstream correlator; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)

	int Lambda = 1;
	string filename = "Data/GrepLambda/Correlator_Lambda" + to_string(Lambda)
			+ "_M" + to_string(GAUSS_RANK) + ".dat";
	correlator.open(filename, mode);
	correlator.precision(15);

	const double X_LIMITS = 5.0 ;// * KF();
	const double T_LIMITS = 15.0 ;// * Energy(Q_momenta(KF()));


	for (double time = 0.01*T_LIMITS; time < T_LIMITS; time += 0.1) {
		correlator << "\"t=" << time << "\"\n" ;

		for (double coordinate = -X_LIMITS; coordinate <= X_LIMITS; coordinate += 0.1) {
			cout << "x = " << coordinate << " t = " << time << endl;
			Cplx result = GrepLambda(Lambda,{X_coordinate(coordinate), T_time(time)});

			correlator << coordinate << "\t" << real(result) << "\t" << imag(result)
					<< "\t" << time << endl;
		}
		correlator << "\n\n" ;
	}
	correlator << endl;
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
		SpaceTime sp = {X_coordinate(coordinate), T_time(T_LIMITS)};
		//Cplx result = Asymptotics(sp);
		Cplx result = Grep(sp);
		tslice << coordinate << "\t" << real(result) << "\t" << imag(result) << endl;
	}
}

void XsliceCurve(double x){
	ofstream xslice; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "xslice_" + to_string((int)x)
		+ "_" + to_string(GAUSS_RANK) + ".dat";
	xslice.open("Data/"+filename, mode);
	xslice.precision(15);

	const double X_LIMITS = x ;
	const double T_LIMITS = 20.0 ;

	for (double time = 0.25; time < T_LIMITS; time += 0.5) {
		//xslice << "\"t=" << time << "\"\n";
		SpaceTime sp = { X_coordinate(X_LIMITS), T_time(time) };
		//Cplx result = Asymptotics(sp);
		Cplx result = Grep(sp);
		xslice << time << "\t" << real(result) << "\t" << imag(result) << endl;
		cout << "time = " << time << " /" << T_LIMITS << endl;
	}
}



void Determinant(){

	double x = 1.0;

	ofstream determinant; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	//+ to_string((int)x)
	string filename = "Determinant_"  + to_string(GAUSS_RANK) + ".dat";
	determinant.open("Data/"+filename, mode);
	determinant.precision(15);

	const double X_LIMITS = x ;
	const double T_LIMITS = 0.001 ;
	const double Lambda = 1.6;

	for (double time = 1E-6; time < T_LIMITS; time += 1E-5) {
		SpaceTime sp = { X_coordinate(X_LIMITS), T_time(time) };
		auto [V, W] = Determinants(Lambda, sp);
		determinant << time << "\t" << abs(V) << "\t" << abs(W) << endl;
		cout << "time = " << time << " /" << T_LIMITS << endl;
	}
}

void Gx_sum(){


	ofstream data; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	//+ to_string((int)x)
	string filename = "Gx_sum_"  + to_string(GAUSS_RANK) + ".dat";
	data.open("Data/"+filename, mode);
	data.precision(15);

	const double X_LIMITS = 0.05 ;
	const double T_LIMITS = 0.00001 ;

	Cplx result = 0.0;

	const double dx = 0.00001;

	double flag = 1;

	for (double coordinate = -X_LIMITS; coordinate <= X_LIMITS; coordinate += dx) {
		cout << coordinate << " / " << X_LIMITS  << '\n';
		SpaceTime sp = {X_coordinate(
				coordinate + flag*Cplx_i * coordinate/(1.0 + pow(coordinate,2))
				)
				, T_time(T_LIMITS)};
		Cplx j_cplx = - ( pow(coordinate,2) -1.0)/pow(1.0 +  pow(coordinate,2),2);
		Cplx Jacobian =  (1.0 + flag*Cplx_i * j_cplx) * dx;
		Cplx tmp = Grep(sp) * Jacobian;
		data << coordinate << "\t" << real(tmp) << "\t" << imag(tmp) << "\n";
		cout << tmp << '\n';
		result += tmp;
		//data << coordinate << "\t" << real(result) << "\t" << imag(result) << endl;
	}
	cout << "integral of Gx = " << result << ", abs = " << abs(result) << endl;

	//0.1   -> (0.743116,-0.249829)
	//0.01  -> (0.913798,-0.0685808)
	//0.001 -> (0.985748,-0.032087)
}


int main() {
	LOG_DURATION("Total");
	//UnitTests();
	//CorrelatorCurve();

	//Integrating();// verifies diff on answers btw Gauss and Gauss_Kronrod



	//LambdaCurve();
	//CorrelatorCurve();

	//XsliceCurve(0.0); //
	//TsliceCurve(10.0);//


	//Fourier1D();
	//Fourier2D();
	//Gpt();
	//Determinant();
	Gx_sum();


//	double eta = 1.0;
//	Q_momenta k(1.14 );
//	Q_momenta k_prime(1.14 - 1e-10);
//	Cplx elem = Nu_matrix_elem(eta,  k,  k_prime,  spacetime);
//	Cplx diag = Nu_diagonal_l( eta,  k,  spacetime);
//
//	cout << "elem = " << elem << "\ndiag = " << diag << endl;
//	cout << setprecision(16) << "diff = " << elem -  diag << endl;
//


//	cout << "E_PV_lattice =" << EPV_l (q, spacetime) << endl;
//	cout << "E_derivative_lattice =" << EPV_derivative_l (q, spacetime) << endl;
//	cout << "Diagonal element =" << Nu_diagonal_l(eta, q, spacetime);

	return 0;
}
