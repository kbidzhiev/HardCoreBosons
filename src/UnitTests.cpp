#include "kernel.hpp"

using namespace std;

void PV_comparison(){
	Q_momenta q(10.5);
	SpaceTime st(X_coordinate(1.7),T_time(10.6));
	Cplx old_value = PrincipalValue_old(q, st)/M_PI;
	Cplx new_value = PrincipalValue(q, st)/M_PI;
	Cplx erf_arg = (st.x - q.value * st.t) * (-1.0 + Cplx_i) / (2.0 * sqrt(st.t));

	cout << "PV_old = " << old_value << '\n' <<
			"PV_new = " << new_value << endl;
	cout << "erf          = " << Faddeeva::erf(erf_arg) << '\n';
	cout << "erf(wolfram) = (1.00787,-0.0223579)" << endl;

}


void UnitTests(){
	PV_comparison();
}
