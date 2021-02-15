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




	for (double time = 0.0; time < 1.0; time += 0.1) {
		correlator << "\"time=" << time << "\"\n" ;
		for (double coordinate = -2.0; coordinate < 2.0; coordinate += 0.1) {

			//SpaceTime st(X_coordinate(coordinate), T_time(time));
			Cplx result = Grep({X_coordinate(coordinate), T_time(time)});

//			SpaceTime st(X_coordinate(0.1),T_time(0.1));
//			auto [result,w] = Determinants (l, st);
//			Cplx coeff = G0 (st) * 2.0/(l * l + 1.0) ;
//			result *= (coeff - 1.0);
//			result += w;

			correlator << coordinate << "\t" << real(result) << "\t" << imag(result)
					<< "\t" << time << endl;
		}
		correlator << "\n\n" ;
	}

	return 0;
}
