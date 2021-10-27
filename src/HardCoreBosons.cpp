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

void Gxt_sum(){

	ofstream data_profile, data; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)

	string filename_profile = "Gxt_profile_"  + to_string(GAUSS_RANK) + ".dat";
	data_profile.open("Data/" + filename_profile, mode);
	data_profile.precision(15);


	string filename = "Gxt_sum_"  + to_string(GAUSS_RANK) + ".dat";
	data.open("Data/" + filename, mode);
	data.precision(15);

	const double X_LIMITS = 15.0 ;

	const double T_min = 0.01 ;
	const double T_max = 10.0 ;
	const double dt = 0.01;
	const int n_max = (T_max - T_min)/dt;

	map<double, Cplx> m_Gp0t;

	int counter = 0;

	const double dx = 0.01;

	double deform_contour = 0;

#pragma omp parallel for num_threads(omp_get_num_procs())
	for (int n = 0; n <= n_max; ++n) {
		T_time time(T_min + n * dt);
		data_profile << "\"t = " << time.value << "\"\n";


		Cplx result = 0.0;

		for (double x = -X_LIMITS; x <= X_LIMITS; x += dx) {
			X_coordinate coord (x + deform_contour * Cplx_i * x/ (1.0 + pow(x, 2)));

			SpaceTime sp(coord, time);

			Cplx j_cplx = -(pow(coord.value, 2) - 1.0) / pow(1.0 + pow(coord.value, 2), 2);
			Cplx Jacobian = (1.0 + deform_contour * Cplx_i * j_cplx) * dx;



			Cplx tmp = Grep(sp) * Jacobian;

			data_profile << x << '\t'
					 << real(tmp) << '\t'
					 << imag(tmp) << '\t'
					 << time.value << "\n";
			//cout << tmp << '\n';
			result += tmp;
			//data_profile << coordinate << "\t" << real(result) << "\t" << imag(result) << endl;
		}
		data_profile << "\n\n";

		m_Gp0t[time.value] = result;

		cout << static_cast<double>(counter++)/n_max << " %\n";
	}

	for(auto &[t,G] : m_Gp0t){
		data << t << '\t'
				<< real(G) << '\t'
				<< imag(G) << '\n';
	}

	cout << "Done !" << endl;

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
	Gxt_sum();


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
