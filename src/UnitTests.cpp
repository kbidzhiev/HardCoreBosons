#include "kernel.hpp"

using namespace std;

void Integrating(){
	// Here I'm checking quality between Gauss and Gauss_Kronrod.
	// Apparently Gauss requires large matrices to reach the same answer as Gauss_Kronrod has
	// for small matrices
	// M = 11 gives diff answers,
	// whereas M = 61 is ok and agrees with Mathematica
	using namespace boost::math::quadrature;
	auto f=[](const double x){
		return cos(x);
	};
	double error ;
	double kron = gauss_kronrod<double, 11>::integrate(f, -50, 50,  10, 1e-9, &error);
	double g = gauss<double, 11>::integrate(f, -50, 50);
	cout << kron << "\n" << g << endl;
}


void PV_comparison(){
	Q_momenta q(0.5);
	const SpaceTime st(X_coordinate(1.7),T_time(0.0001));
	Cplx old_value = 0;//PrincipalValue_old(q, st);///M_PI;
	Cplx new_value = PrincipalValue(q, st);///M_PI;
	Cplx erf_arg = (st.x - q.value * st.t) * (1.0 - Cplx_i) / (2.0 * sqrt(st.t));

	cout << "PV_old = " << old_value << '\n' <<
			"PV_new = " << new_value << endl;
	cout << "erf          = " << Faddeeva::erf(erf_arg) << '\n';
	cout << "erf(wolfram) = (0.99916,0.00461777)" << endl;
}


void UnitTests(){
	//Integrating();
	PV_comparison();
}
