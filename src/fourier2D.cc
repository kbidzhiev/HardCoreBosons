#include <iostream>
#include <fstream>

#include "kernel.hpp"
#include "fftw_2d_cpp.hh"
#include "fftw_cpp.hh"

using namespace std;

//#include <complex.h>
//#include <fftw3.h> //fastEST Fourier transformation

//https://github.com/lschw/fftw_cpp


double EF = KF() * KF() ; // fermi energy

Cplx Asymptotics (double x, double t) {
	x *= KF();
	t *= EF;

	Cplx log_term = log(2 * t) + Cplx_i * M_PI * 0.5;
	Cplx result = -0.5 * x * x / log_term - Cplx_i * 0.5 * t;
	result = exp(result);
	result /= (sqrt(log_term) * sqrt(Cplx_i * t));
	result *= 3.5; // 3.5 \approx all the constants;
	return result;
};

void Fourier2D() {

	size_t N1 = 500;
	size_t N2 = 500;
	size_t N = N1 * N2;
	dcvector data(N);
	dcvector data_fft(N);
	dvector x1(N1);
	dvector t2(N2);
	dvector k1(N1);
	dvector w2(N2);

	// Create test data
	//double w1 = 2 * M_PI;
	//double w2 = 2 * M_PI * 3;
	//double w3 = 2 * M_PI * 5;
	//double w4 = 2 * M_PI * 7;
	double xmax1 = 40.0;
	double tmax2 = 20.0;
	size_t counter = 0;
//	auto Box = [&](double x, double y){
//		if (x > 0 && y > 0){
//			return 1.0;
//		} else
//			return 0.0;
//	};



//#pragma omp parallel for num_threads(omp_get_num_procs()) collapse(2)
	for (size_t i = 0; i < N1; ++i) {
		for (size_t j = 0; j < N2; ++j) {
			x1[i] = i * xmax1 / N1 - xmax1 / 2;
			t2[j] = j * tmax2 / N2;	//- tmax2 / 2.0;
			SpaceTime st(X_coordinate(x1[i]), T_time(t2[j]));
			cout << "x= " << i + 1 << " / " << N1 << " ; t= " << j + 1 << " / "
					<< N2 << " ;\t" << ++counter << " / " << N << endl;

			data[i * N2 + j] = t2[j] >= 0.5 ? Grep(st) : 0;
			//data[i * N2 + j] = t2[j] >= 0.5 ? Asymptotics(x1[i],t2[j]) : 0;


		}
	}

	auto ValueFilter = [](double x) {

//    	double cutoff = 10.0;
//    	if (abs(x)<1e-14) return 0.0;
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
	fh1.open("Data/data2d.dat");
	fh2.open("Data/data2d_fft.dat");
	fh1 << "# x \tt \tRe[f(x, t)] \tIm[f(x,t)]\n";
	fh2 << "# k \tw \tRe[f(k, w)] \tIm[f(x,t)]\n";
	for (size_t i = 0; i < N1; ++i) {
		for (size_t j = 0; j < N2; ++j) {
			fh1 << x1[i] << " \t";
			fh1 << t2[j] << " \t";
			fh1 << data[i * N2 + j].real() << " \t";
			fh1 << data[i * N2 + j].imag() << "\n";
			//if (k1[i]>=0 && w2[j]>=0)
			{
				fh2 << k1[i] / (2 * M_PI) << " \t";
				fh2 << w2[j] / (2 * M_PI) << " \t";
				fh2 << ValueFilter(data_fft[i * N2 + j].real())
						* (4 * xmax1 * tmax2) << " \t";
				fh2	<< ValueFilter(data_fft[i * N2 + j].imag())
						* (4 * xmax1 * tmax2) << "\n";

			}
		}
		fh1 << "\n";
		fh2 << "\n";
	}
	fh1.close();
	fh2.close();
	cout << "DONE" << endl;
}

void Fourier1D() {
	size_t N = 1000;
	dcvector data(N);
	dcvector data_fft(N);
	dvector t(N);
	dvector f(N);

	double xmax = 50.0;
	double time = 1.0;


	for (size_t i = 0; i < N; ++i) {
		t[i] = i * xmax / N - xmax / 2;
		SpaceTime st(X_coordinate(t[i]), T_time(time));
		//data[i] = Asymptotics (st.x, st.t);
		data[i] = Grep(st) ;  // Here we do Fourier for a fixed time
		// no need to introduce small time truncation
		cout << "i = " << i << " / " << N << endl;
	}

	// Fourier transform
	FFT fft(N, xmax);
	fft.fft(data, data_fft);
	fft.freq(f);
	fft.shift_freq(f, data_fft);

	// Save
	std::ofstream fh1;
	std::ofstream fh2;
	fh1.open("Data/Gp/Gx_time" + to_string((int) time) + "Gauss" + to_string(GAUSS_RANK) +".dat");
	fh2.open("Data/Gp/Gp_time" + to_string((int) time) + "Gauss" + to_string(GAUSS_RANK) +".dat");
	fh1 << "# \tx \tRe[f(x)] \tIm[f(x)]\n";
	fh2 << "# \tf \tRe[f(w)] \tIm[f(w)]\n";
	for (size_t i = 0; i < N; ++i) {
		fh1 << t[i] << " \t";
		fh1 << data[i].real() << "\t" << data[i].imag() << "\n";
		fh2 << f[i] / (2 * M_PI) << " \t";
		fh2 << data_fft[i].real() * 2 * xmax << " \t"
				<< data_fft[i].imag() * 2 * xmax << "\n";
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

