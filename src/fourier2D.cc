#include <iostream>
#include <fstream>

#include "kernel.hpp"
#include "fftw_2d_cpp.hh"
#include "fftw_cpp.hh"

#include <stdio.h>
#include <math.h>


using namespace std;

//#include <complex.h>
//#include <fftw3.h> //fastEST Fourier transformation

//https://github.com/lschw/fftw_cpp


double EF = KF() * KF() / 2.0 ; // fermi energy

double Box (SpaceTime st){
	if ( -0.5 < st.x && st.x < 0.5
			&& -0.5 < st.t && st.t < 0.5 ){
		return 1;
	} else return 0;
}

Cplx Gauss (SpaceTime st){
	return exp(-0.5*(st.x*st.x + st.t*st.t) );
}


Cplx Asymptotics (SpaceTime st) {
	double x = st.x;
	double t = st.t;

	x *= KF();
	t *= EF;

	Cplx log_term = log(2 * t) + Cplx_i * M_PI * 0.5;
	Cplx result = -0.5 * x * x / log_term - Cplx_i * 0.5 * t;
	result = exp(result);
	result /= (sqrt(log_term) * sqrt(Cplx_i * t));
	result *= 3.5; // 3.5 \approx all the constants;

	return result;
}

void Fourier2D() {

	size_t N1 = 30;
	size_t N2 = 100;
	size_t N = N1 * N2;
	dcvector data(N);
	dcvector data_fft(N);
	dvector x1(N1);
	dvector t2(N2);
	dvector k1(N1);
	dvector w2(N2);


	const double xmax1 = 30.0;
	const double tmax2 = 30.0;
	size_t counter = 0;

	const double dx = xmax1 / N1;//Should be 1 for lattice
	const double dt = tmax2 / N2;
	const double shift_x = xmax1 / 2.0;
	const double shift_t = tmax2 / 2.0;

	const double truncation = 0.0; //0.25
	auto ExcludeVicinityOf0 = [&](double t) {
		//Here I isolate vicinity of 0 into 4 regions
		if (t <= -truncation || t >= truncation) { 			// (-T: -truncation]
			return t;
		} else { //(-truncation : 0)
			return truncation;
		}
	};


//#pragma omp parallel for num_threads(omp_get_num_procs()) collapse(2)
	for (size_t j = 0; j < N2; ++j) {
		for (size_t i = 0; i < N1; ++i) {


//#pragma omp parallel for num_threads(omp_get_num_procs())
			t2[j] = j * dt - shift_t;
			x1[i] = i * dx - shift_x;



			SpaceTime st(X_coordinate(x1[i]), T_time(abs(
					ExcludeVicinityOf0(t2[j])
					)));
			//the line replaced by abs in T_time ctr
			//if (t2[j] < 0) 	st.t = -st.t; // t2[j] is still negative, but st.t is positive;



			//const Cplx result = 2 * M_PI * GrepEta_l(0, st);
			const Cplx result = Grep_l(st); //- G0(st);
			//cout << result << endl;
			//const Cplx result = Gauss(st);
			//Cplx result = Asymptotics(st);
			//Cplx result = Box(st);
			//double lambda = 1.0;
			//const Cplx result = GrepLambda(lambda, st);


			++counter;
			cout << "x= " << i + 1 << " / " << N1 << ";\t"
				 << "t= " << j + 1 << " / " << N2 << ";\t"
				 << "Grep= " << result << ";\t"
				 << counter << " / " << N
				 << "\t" <<(double)counter/N <<endl;
			//terminate();

			//Here I isolate vicinity of 0 into 4 regions
			if (t2[j] <= 0) { 						// (-T: -truncation]
				data[i * N2 + j] = conj(result); 	// conj(result)
			} else {								// [truncation: T)
				data[i * N2 + j] = result;
			}


//			auto Positive = [&]() {
//				//Here I isolate vicinity of 0 into 4 regions
//				if (t2[j] >= 0 && t2[j] < truncation) { //[0: //truncation)
//					data[i * N2 + j] = 0;  //conj(tmp); // conj(tmp)
//				} else {									// [truncation: T)
//					data[i * N2 + j] = result;
//				}
//			};


		}

	}

	auto ValueFilter = [](double x) {
//		if (abs(x) > 0.8){
//			return x;
//		} else {
//			return 0.0;
//		}

//    	double cutoff = 0.5;
//    	if (abs(x)<1e-16) return 0.0;
//    	else if (abs(x) > cutoff){
//    		if (x > cutoff) return cutoff;
//    		else  return -cutoff;
//    	}
		return x;
	};

	// Fourier transform
	FFT2D fft(N1, N2, xmax1, tmax2);
	fft.fft(data, data_fft);
	fft.freq1(k1);
	fft.freq2(w2);
	fft.shift_freq(k1, w2, data_fft);

	// Save
	std::ofstream fh1;
	std::ofstream fh2;
	std::ofstream fourier_momentum0;
	std::ofstream function;
	fh1.open("Data/data2d_Gauss" + to_string(GAUSS_RANK)+ ".dat");
	fh2.open("Data/data2d_fft_Gauss" + to_string(GAUSS_RANK)+ ".dat");
	fourier_momentum0.open("Data/momentum0_fft_Gauss" + to_string(GAUSS_RANK)+ ".dat");
	function.open("Data/function_Gauss" + to_string(GAUSS_RANK)+ ".dat");

	fh1 << "# x \tt \tRe[f(x, t)] \tIm[f(x,t)]\n";
	fh2 << "# k \tw \tRe[f(k, w)] \tIm[f(x,t)]\n";
	for (size_t i = 0; i < N1; ++i) {
		for (size_t j = 0; j < N2; ++j) {
			fh1 << x1[i] << " \t";
			fh1 << t2[j] << " \t";
			fh1 << ValueFilter(data[i * N2 + j].real()) << " \t";
			fh1 << ValueFilter(data[i * N2 + j].imag()) << "\n";
			//if (k1[i]>=0 && w2[j]>=0)
			{
				fh2 << k1[i]/KF() << " \t";
				fh2 << w2[j]/EF << " \t";
				fh2 << pow(-1, i+j) * ValueFilter(
						data_fft[i * N2 + j].real() * xmax1 * tmax2 / (2.0 * M_PI)
						)<< " \t";
				fh2	<< pow(-1, i+j) * ValueFilter(
						data_fft[i * N2 + j].imag() * xmax1 * tmax2 / (2.0 * M_PI)
						)<< "\n";
				if(k1[i]  == 0){
					fourier_momentum0 << w2[j]  << " \t"
						<< pow(-1, i+j) * ValueFilter(
							data_fft[i * N2 + j].real() * xmax1 * tmax2 / (2.0 * M_PI)
							)<< " \t"
						<< pow(-1, i+j) * ValueFilter(
							data_fft[i * N2 + j].imag() * xmax1 * tmax2 / (2.0 * M_PI)
							) << "\n" ;
				}
				if(x1[i]  == 0){
					function << t2[j]  << " \t"
						<< ValueFilter(data[i * N2 + j].real()) << " \t"
						<< ValueFilter(data[i * N2 + j].imag()) << "\n" ;
				}
			}
		}
		fh1 << "\n";
		fh2 << "\n";
	}
	fh1.close();
	fh2.close();
	fourier_momentum0.close();
	cout << "DONE" << endl;
}

void Fourier1D() {
	size_t N = 101;
	dcvector data(N);
	dcvector data_fft(N);
	dvector t(N);
	dvector f(N);

	double xmax = 0.0;
	double tmax = 60.0;

	size_t counter = 0;
#pragma omp parallel for num_threads(omp_get_num_procs())
	for (size_t i = 0; i < N; ++i) {
		t[i] = i * tmax / N - tmax / 2.0;
		SpaceTime st(X_coordinate(xmax), T_time(
				abs(t[i]))
				);
		//data[i] = Gauss(st);//Asymptotics (st.x, st.t);
		Cplx result = GrepEta_l(2.0, st);

		if (t[i] <= 0) { 						// (-T: -truncation]
			data[i] = conj(result); 	// conj(result)
		} else {								// [truncation: T)
			data[i] = result;
		}


		//data[i] = Grep(st) ;  // Here we do Fourier for a fixed time
		// no need to introduce small time truncation
		cout << "i = " << ++counter << " / " << N << endl;
	}

	// Fourier transform
	FFT fft(N, tmax);
	fft.fft(data, data_fft);
	fft.freq(f);
	fft.shift_freq(f, data_fft);

	// Save
	std::ofstream fh1;
	std::ofstream fh2;
	fh1.open("Data/Gp/Gx_time" + to_string((int) tmax) + "Gauss" + to_string(GAUSS_RANK) +".dat");
	fh2.open("Data/Gp/Gp_time" + to_string((int) tmax) + "Gauss" + to_string(GAUSS_RANK) +".dat");
	fh1 << "# \tx \tRe[f(x)] \tIm[f(x)]\n";
	fh2 << "# \tf \tRe[f(w)] \tIm[f(w)]\n";
	for (size_t i = 0; i < N; ++i) {
		fh1 << t[i] << " \t";
		fh1 << data[i].real() << "\t" << data[i].imag() << "\n";
		fh2 << f[i]  << " \t";
		fh2 << pow(-1,i) * data_fft[i].real() * tmax/ pow(2 * M_PI, 0.5) << " \t"
			<< pow(-1,i) * data_fft[i].imag() * tmax/ pow(2 * M_PI, 0.5) << "\n";
	}
	fh1.close();
	fh2.close();

	cout << "Fourier 1D DONE" << endl;
}

void Gpt() {
	size_t N = 100;
	dcvector data(N);
	dcvector data_fft(N);
	dvector t(N);
	dvector f(N);
	// Save
	std::ofstream fh1;
	std::ofstream fh2;
	fh1.open("Data/Gp/Gxt.dat");
	fh2.open("Data/Gp/Gpt.dat");
	fh1 << "# \tx \tRe[f(x)] \tIm[f(x)]\n";
	fh2 << "# \tf \tRe[f(w)] \tIm[f(w)]\n";

	double xmax = 10.0;
	double timemax = 20.;
	for (double time = 0.0; time < timemax; time += 0.1) {
		for (size_t i = 0; i < N; ++i) {
			t[i] = i * xmax / N - xmax / 2;
			SpaceTime st(X_coordinate(t[i]), T_time(time));
			data[i] = time >= 0.1 ? Grep(st) : 0;    //sin(w1*t[i]) + sin(w2*t[i]);
			cout << "i = " << i << " / " << N
					<<" time = " << time << " / " << timemax
					<< endl;
		}
		FFT fft(N, xmax);
		fft.fft(data, data_fft);
		fft.freq(f);
		fft.shift_freq(f, data_fft);


		int i = N / 2 -1;
		fh1 << time << " \t";
		fh1 << data[i].real() << "\t" << data[i].imag() << endl;
		fh2 << time << " \t";
		fh2 << data_fft[i].real() * 2 * xmax << " \t"
				<< data_fft[i].imag() * 2 * xmax << "\t"
				<< f[i]/ (2 * M_PI)
				<< endl;
	}

//	for (size_t i = 0; i < N; ++i) {
//		fh1 << t[i] << " \t";
//		fh1 << data[i].real() << "\t" << data[i].imag() << "\n";
//		fh2 << f[i] / (2 * M_PI) << " \t";
//		fh2 << data_fft[i].real() * 2 * xmax << " \t"
//				<< data_fft[i].imag() * 2 * xmax << "\n";
//	}
	fh1.close();
	fh2.close();

	cout << "Fourier Gpt DONE" << endl;
}


void foo(){
    fftw_complex *in, *out;
    fftw_plan p;
    int w = 16;
    double shift = 0.65 * w;
    double a = 2;
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * w );
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * w );
    for (int i = 0; i < w; i++){
            in[i][0] = exp(-((i - shift)*(i - shift))/a );
            in[i][1] = 0;
            //* real(exp(-2.0 * M_PI * Cplx_i * (Cplx)(shift / w) * (Cplx)i))
    }
    p = fftw_plan_dft_1d(w, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    vector<double> fftout;
    fftout.reserve(w);
    for (int i = 0; i< w; i++){
    	if (i<w/2){
    		fftout[i+w/2] = out[i][0] ;
    	} else {
    		fftout[i-w/2] = out[i][0] ;
    	}
    }
    for (int i = 0; i < w; i++)
        printf("%4d : %+9.4f %+9.4f %+9.4f i\n", i, in[i][0], out[i][0], fftout[i]);
    fftw_destroy_plan(p); fftw_cleanup();
    fftw_free(in); fftw_free(out);
}


