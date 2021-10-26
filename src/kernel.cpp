#include "kernel.hpp"

#include <vector>
#include <omp.h>

TotalDuration elements("Elements");
TotalDuration mf("MatrixFilling");
TotalDuration det("Determinant");

//const size_t s = 2 * g.weights().size() ; //Gauss
const size_t s = 2 * g.weights().size() - 1; //GKrondrod, we always have ODD number of weights


using Eigen::MatrixXcd;
//using boost::math::quadrature::trapezoidal;

double KF(){
	return M_PI * RHO;
}

double Mu_chempot(){
	return KF()*KF() * 0.5 / MASS;

}

double Energy(Q_momenta q_momenta){
	return (q_momenta.value) * (q_momenta.value) * 0.5 / MASS ;
}

complex<double> Tau(Q_momenta q_momenta,  SpaceTime spacetime){
	return spacetime.t * Energy(q_momenta) - spacetime.x * q_momenta.value;
}

double Theta (Q_momenta q_momenta){
	double result = exp(B_BETA * (Energy(q_momenta) - Mu_chempot())) + 1.0;
	return 1./result;
}

Cplx Eminus (Q_momenta q_momenta,  SpaceTime spacetime){
	Cplx result = sqrt(Theta (q_momenta)/M_PI);
	result *= exp(Cplx_i * 0.5 * Tau(q_momenta, spacetime));
	return result;
}


Cplx G0 (SpaceTime spacetime){
	Cplx result = exp(-Cplx_i * 0.25 * M_PI );
	result *= sqrt(MASS / (2.0 * M_PI * spacetime.t));
	result *= exp( Cplx_i * 0.5 * MASS * spacetime.x * spacetime.x / spacetime.t );
	return result;
}


double Q_G (const size_t i) {
	//size_t middle_point = g.abscissa().size();
	return i < g.abscissa().size() ?
			- KF() * g.abscissa()[g.abscissa().size() - i - 1] :
			  KF() * g.abscissa()[i - g.abscissa().size()];
}

double Weight_G (const size_t i) {
	//size_t middle_point = g.weights().size();
	return i < g.weights().size() ?
			KF() * g.weights()[g.weights().size() - i - 1] :
			KF() * g.weights()[i - g.weights().size()];
}

double Q_Kr(const size_t i) {
	//size_t middle_point = g.abscissa().size() - 1;
	return i < g.abscissa().size() - 1 ?
		   - KF() * g.abscissa()[g.abscissa().size() - i - 1] :
			 KF() * g.abscissa()[i - g.abscissa().size() + 1];
}

double Weight_Kr(const size_t i) {
	//size_t middle_point = g.weights().size() - 1;
	return i < g.weights().size() - 1?
			KF() * g.weights()[g.weights().size() -1 - i ] :
			KF() * g.weights()[i - g.weights().size() + 1];
}


pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime spacetime){

	vector<Cplx> e_minus(s,0.0);
	vector<Cplx> e_infty(s,0.0);
	vector<Cplx> e_infty_derivative(s,0.0);
	vector<Cplx> e_plus_for_W(s,0.0);
	//
	vector<Cplx> e_infty_eps(s,0.0);
	vector<Cplx> e_minus_eps(s,0.0);
	//
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);




//	cout << "Weight_Kr\tQabscissa\n";
//	for(size_t i = 0; i< s; ++i){
//		cout << Q_Kr(i) <<"\t" << Weight_Kr(i)<< "\n";
//	}
//	cout << endl;
//	terminate();

	const Cplx first_term_of_PVderivative = -Cplx_i * spacetime.t * 2.0 * M_PI * G0(spacetime) / MASS;


	double epsilon = 1E-8;

//#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(elements);
		Q_momenta q_i(Q_Kr(i)); //!
		Cplx tau = Tau (q_i, spacetime);
		Cplx pv = PrincipalValue(q_i, spacetime);

		Q_momenta q_i_eps(Q_Kr(i) + epsilon); //!
		Cplx tau_eps = Tau (q_i_eps, spacetime);
		Cplx pv_eps = PrincipalValue(q_i_eps, spacetime);


		double sinsin = 1.0/(1.0 + Lambda*Lambda);
		double sincos = -Lambda/(1.0 + Lambda*Lambda);
		double sin = sqrt(sinsin);
		double cos = sincos / sin;

		e_infty[i] = sinsin * (pv/M_PI) + sincos * exp(- Cplx_i * tau);
		e_minus[i] = Eminus(q_i, spacetime) ;

		//
		e_infty_eps[i] = sinsin * (pv_eps/M_PI)	+ sincos * exp(- Cplx_i * tau_eps);
		e_minus_eps[i] = Eminus(q_i_eps, spacetime) ;
		//
		Cplx second_term_of_PVderivative = - Cplx_i * (Q_Kr(i) * spacetime.t / MASS - spacetime.x) * pv;
		Cplx e_inf_der = sinsin * (first_term_of_PVderivative + second_term_of_PVderivative)/M_PI;
		e_inf_der += sincos * exp(-Cplx_i * tau) *(-Cplx_i) * (spacetime.t * Q_Kr(i) / MASS - spacetime.x);

		e_infty_derivative[i] = e_inf_der;

		e_plus_for_W[i] = sin * (pv/M_PI) + cos*exp(- Cplx_i * tau);
		e_plus_for_W[i] /= sqrt(2.0);
		e_plus_for_W[i] *= e_minus[i];



	}


	//time for 2 core proc
	//avr 2125209 	parallel
	//avr 31592		consequential



//#pragma omp parallel for num_threads(omp_get_num_procs()) //collapse(2)
	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(mf);
		for (size_t j = 0; j < s; j++) { // j= i
			Cplx v = e_minus_eps[i] * e_minus[j];
			Cplx w =  e_plus_for_W[i] * e_plus_for_W[j];
			// In order to avoid having 0/0 in W eq (3.25), we
			// cancel sin^2 (eta/2) in denominator and enumerator

			//Cplx w = e_infty[i] * e_minus[i] //E_+ = E_infty E_-
			//					* e_infty[j] * e_minus[j];// * 0.5 * (1.0+Lambda * Lambda) ;



				//v *= e_infty_derivative[i] ;
				Cplx de =  e_infty_eps[i] - e_infty[j];
				double dq =  Q_Kr(i) + epsilon - Q_Kr(j);
				//cout << "de " << de << endl;
				v *= (de/dq);
				V(i, j) = sqrt(Weight_Kr(i)) * (v) * sqrt(Weight_Kr(j)) ; //!
				W(i, j) = sqrt(Weight_Kr(i)) * (v-w) * sqrt(Weight_Kr(j)) ; //!


			if (i == j) {
				V(i, j) += 1.0 ;
				W(i, j) += 1.0 ;
			} else {
//				v *= (e_infty[i] - e_infty[j]) ;
//				v /= Q_Kr(i) - Q_Kr(j);
//
//				V(i, j) = sqrt(Weight_Kr(i)) * v * sqrt(Weight_Kr(j));
//				W(i, j) = sqrt(Weight_Kr(i)) * (v - w) * sqrt(Weight_Kr(j)); //new part

				//V(j, i) = V(i, j);
				//W(j, i) = W(i, j);
			}
		}
	}


	//cout << W.determinant() << "\n\n";
	//cout << W << endl;
	//terminate();

	//time
	//28997668 parallel
	//142652   consequential
	//terminate();

	//W = V - W; // this part can be moved to the loop
	//cout << V - W << endl;
	//cout << W << endl;
	//terminate();


	Cplx detV = V.determinant();
	Cplx detW = W.determinant();

	double Jacobian_from_eta_to_lambda = (2.0)/(Lambda * Lambda + 1.0);
	detV *= Jacobian_from_eta_to_lambda;
	detW *= Jacobian_from_eta_to_lambda;


	return {detV, detW};
}


//for t=0
/*
pair <Cplx, Cplx> Determinants(double Lambda, SpaceTime spacetime){

	vector<Cplx> e_minus(s,0.0);
	vector<Cplx> e_infty(s,0.0);
	vector<Cplx> e_infty_derivative(s,0.0);
	vector<Cplx> e_plus_for_W(s,0.0);
	//
	vector<Cplx> e_infty_eps(s,0.0);
	vector<Cplx> e_minus_eps(s,0.0);
	//
	MatrixXcd V(s, s);
	MatrixXcd W(s, s);


	spacetime.t = 0;



	const Cplx first_term_of_PVderivative = -Cplx_i * spacetime.t * 2.0 * M_PI * G0(spacetime) / MASS;

	double eta = Lambda;

	double epsilon = 1E-8;

//#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(elements);
		Q_momenta q_i(Q_Kr(i)); //!
		Cplx tau = Tau (q_i, spacetime);
		Cplx pv = Cplx_i * M_PI * exp (-Cplx_i * Tau (q_i, spacetime));

		Q_momenta q_i_eps(Q_Kr(i) + epsilon); //!
		Cplx tau_eps = Tau(q_i_eps, spacetime);
		Cplx pv_eps = Cplx_i * M_PI * exp (-Cplx_i * Tau (q_i_eps, spacetime));




		e_infty[i] = pow(sin(0.5 * eta), 2) * (pv/M_PI)
				+ sin(0.5 * eta)*cos(0.5 * eta) * exp(- Cplx_i * tau);
		e_minus[i] = Eminus(q_i, spacetime) ;

		//
		e_infty_eps[i] = pow(sin(0.5 * eta), 2) * (pv_eps/M_PI)
				+ sin(0.5 * eta)*cos(0.5 * eta) * exp(- Cplx_i * tau_eps);
		e_minus_eps[i] = Eminus(q_i_eps, spacetime) ;
		//


		e_plus_for_W[i] = (sin(0.5 * eta) * (pv/M_PI)
						+ cos(0.5 * eta) * exp(- Cplx_i * tau)) * Eminus(q_i, spacetime);


	}


	for (size_t i = 0; i < s; i++) {
		ADD_DURATION(mf);
		for (size_t j = 0; j < s; j++) { // j= i
			Cplx v = e_minus_eps[i] * e_minus[j];
			Cplx w =  e_infty_eps[i] * e_infty[j] * e_minus_eps[i] * e_minus[j]
						/( 2.0 ) ;
			// In order to avoid having 0/0 in W eq (3.25), we
			// cancel sin^2 (eta/2) in denominator and enumerator

				//v *= e_infty_derivative[i] ;
				Cplx de =  e_infty_eps[i] - e_infty[j];
				double dq =  Q_Kr(i) + epsilon - Q_Kr(j);
				//cout << "w = " << w << endl;
				v *= (de/dq);
				V(i, j) = sqrt(Weight_Kr(i)) * (v) * sqrt(Weight_Kr(j)) ; //!
				W(i, j) = sqrt(Weight_Kr(i)) * (v-w) * sqrt(Weight_Kr(j)) ; //!


			if (i == j) {
				V(i, j) += 1.0 ;
				W(i, j) += 1.0 ;
			} else {

			}
		}
	}



	Cplx detV;
	Cplx detW;
	{
		ADD_DURATION(det);
		detV = V.determinant();
		detW = W.determinant();
	}

	return {detV, detW};
}

*/

Cplx GrepLambda(double Lambda,  SpaceTime st){
	auto [detv,detw] = Determinants(Lambda, st);
	detv *= (G0(st) - 1.0);
	return (detv + detw);
}

Cplx Grep(SpaceTime st){
	auto f= [&](double Lambda){
		return GrepLambda(Lambda, st);
	};

	double error = 100;
	/*
	 * execution time for LambdaCurve()
	 * depth = 10;  time = 4 s;
	 * depth = 5 ;  time =
	 * */
	Cplx result = gauss_kronrod<double, 61>::integrate(f, -100, 100,  10, 1e-9, &error);

	//Cplx result = gauss_kronrod<double, 61>::integrate(f, -M_PI, M_PI,  10, 1e-9, &error);


	//cout << "error value = " <<  error << endl;

	return 0.5*result/M_PI;
}



