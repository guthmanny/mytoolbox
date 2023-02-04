//
// Created by wcy on 2023/1/4.
//

#ifndef MYFFT_FFT_H
#define MYFFT_FFT_H

#include <complex>
#include <cstring>
#include <iostream>

namespace fft {

    class fft {

    };

    void fft_test();

    int decompose(int N, int *p);

    void fft_2d(const int *in, int *out, int N);

    void my_fft(std::complex<double> *in, std::complex<double> *out, int N, int *p, int n);

    void my_fft1(std::complex<double> *in, std::complex<double> *out, int N);

} // fft

#endif //MYFFT_FFT_H
