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


struct Parameters{
	const double mass = 1.0;

};
using namespace std;

int main() {
	LOG_DURATION("DET");

	ofstream correlator; //here I'm defining output streams, i.e. files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	string filename = "Correlator" + to_string(2 * g.weights().size()) + ".dat";
	correlator.open(filename, mode);
	correlator.precision(15);
	//correlator << "#x \t correlator \t time \n";

	//SpaceTime st(X_coordinate(20.),T_time(1.));




	for(double l = -1.0 ; l < 1.0; l += 0.1){

		SpaceTime st(X_coordinate(3.0),T_time(l));
		Cplx result = Grep(st);

//		SpaceTime st(X_coordinate(0.1),T_time(0.1));
//		auto [result,w] = Determinants (l, st);
//		Cplx coeff = G0 (st) * 2.0/(l * l + 1.0) ;
//		result *= (coeff - 1.0);
//		result += w;

		correlator << l << "\t" << real(result)<< "\t" << imag(result) << "\t" << l  << endl;
	}

	return 0;
}
