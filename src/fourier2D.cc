#include <iostream>
#include <fstream>

#include "kernel.hpp"
#include "fftw_2d_cpp.hh"
#include "fftw_cpp.hh"

using namespace std;

//#include <complex.h>
//#include <fftw3.h> //fastEST Fourier transformation




void Fourier2D(){

	size_t N1 = 100;
	size_t N2 = 100;
	size_t N = N1 * N2;
	dcvector data(N);
	dcvector data_fft(N);
	dvector t1(N1);
	dvector t2(N2);
	dvector f1(N1);
	dvector f2(N2);

	// Create test data
	//double w1 = 2 * M_PI;
	//double w2 = 2 * M_PI * 3;
	//double w3 = 2 * M_PI * 5;
	//double w4 = 2 * M_PI * 7;
	double tmax1 = 10;
	double tmax2 = 10;
	size_t counter = 0;
	for (size_t i = 0; i < N1; ++i) {
		t1[i] = i * tmax1 / N1;
		for (size_t j = 0; j < N2; ++j) {
			t2[j] = j * tmax2 / N2;
			SpaceTime st(X_coordinate(t1[i]),T_time(t2[j] + 0.01));
			cout << "x= " << i << " / " << N1
			<< " ; t= " << j << " / " << N2 << ";\t"
			<< ++counter << " / " << N - N1  << endl;

			data[i * N2 + j] = Grep(st);
					//(sin(w1 * t1[i]) + sin(w2 * t1[i]))
					//* (sin(w3 * t2[j]) + sin(w4 * t2[j]));
		}
	}

	// Fourier transform
	FFT2D fft(N1, N2, tmax1, tmax2);
	fft.fft(data, data_fft);
	fft.freq1(f1);
	fft.freq2(f2);
	fft.shift_freq(f1, f2, data_fft);

	// Save
	std::ofstream fh1;
	std::ofstream fh2;
	fh1.open("Data/data2d.dat");
	fh2.open("Data/data2d_fft.dat");
	fh1 << "# x \tt \tRe[f(x, t)] \tIm[f(x,t)]\n";
	fh2 << "# k \tw \tRe[f(k, w)] \tIm[f(x,t)]\n";
	for (size_t i = 0; i < N1; ++i) {
		for (size_t j = 0; j < N2; ++j) {
			fh1 << t1[i] << " \t";
			fh1 << t2[j] << " \t";
			fh1 << data[i * N2 + j].real() << " \t";
			fh1 << data[i * N2 + j].imag() << "\n";
			fh2 << f1[i] << " \t";
			fh2 << f2[j] << " \t";
			fh2 << data_fft[i * N2 + j].real() << " \t";
			fh2 << data_fft[i * N2 + j].imag() << "\n";
		}
		fh1 << "\n";
		fh2 << "\n";
	}
	fh1.close();
	fh2.close();
	cout << "DONE" << endl;
}



void Fourier1D(){
    size_t N = 10'000;
    dcvector data(N);
    dcvector data_fft(N);
    dvector t(N);
    dvector f(N);

    // Create some test data
    //double w1 = 2*M_PI;
    //double w2 = 2*M_PI*3;
    double tmax = 5.0;
    for(size_t i = 0; i < N; ++i) {
        t[i] = i*tmax/N;
        SpaceTime st(X_coordinate(t[i]),T_time(10.0));
        data[i] = Grep(st);//sin(w1*t[i]) + sin(w2*t[i]);
        cout << "i = " << i << " / " << N << endl;
    }

    // Fourier transform
    FFT fft(N, tmax);
    fft.fft(data, data_fft);
    fft.freq(f);
    fft.shift_freq(f, data_fft);

    // Save
    std::ofstream fh1;
    std::ofstream fh2;
    fh1.open("Data/data.dat");
    fh2.open("Data/data_fft.dat");
    fh1 << "# \tx \tRe[f(x)] \tIm[f(x)]\n";
    fh2 << "# \tf \tRe[f(w)] \tIm[f(w)]\n";
    for(size_t i = 0; i < N; ++i) {
        fh1 << t[i] << " \t";
        fh1 << data[i].real() << " \t" << data[i].imag()<< "\n";
        fh2 << f[i] << " \t";
        fh2 << data_fft[i].real() << " \t" << data_fft[i].imag()<< "\n";
    }
    fh1.close();
    fh2.close();

    cout << "Fourier 1D DONE" << endl;
}
