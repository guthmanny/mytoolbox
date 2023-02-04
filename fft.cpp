//
// Created by wcy on 2023/1/4.
//


#include "fft.h"

namespace fft {

    int decompose(int N, int *p) {
        int n = 0;
        while (N % 2 == 0) {
            N /= 2;
            p[n++] = 2;
        }
        for (int i = 3; i * i <= N; i += 2) {
            while (N % i == 0) {
                N /= i;
                p[n++] = i;
            }
        }
        if (N > 1) {
            p[n++] = N;
        }
        return n;
    }

    void fft_test() {


        const int N = 4521;
        int p[20] = {0};

        std::complex<double> in[N], out[N];
        for (int i = 0; i < N; ++i) {
            in[i] = std::complex<double>(i, 0);
        }
        my_fft1(in, out, N);
        for (int i = 0; i < N; ++i) {
            std::cout<<"in[" << i << "]\t= " << in[i]<<std::endl;
            std::cout << "out[" << i << "]\t= " << out[i] << std::endl;
        }
    }

    void fft_2d(const int *in, int *out, int N) {
        int p[10] = {0};
        int n = decompose(N, p);

        int pos[N];
        for (int i = 0; i < N; ++i) {
            int temp = i;
            pos[i] = 0;

            for (int j = 0; j < n; ++j) {
                pos[i] = pos[i] * p[j] + temp % p[j];
                temp /= p[j];
            }
        }

        int temp[N];
        for (int i = 0; i < N; ++i) {
            temp[i] = in[pos[i]];
        }

        int step = 1;
        for (int i = 0; i < n; ++i) {
            step *= p[n - 1 - i];

            for (int j = 0; j < N / step; ++j) {
                for (int k = 0; k < step; ++k) {
                    out[j * step + k] += temp[j * step + k];
                }
            }
        }
    }

    void my_fft(std::complex<double> *in, std::complex<double> *out, int N, int *p, int n) {
        int P = p[n - 1];
        int Q = N / P;
        for (int a = 0; a < P; ++a) {
            for (int b = 0; b < Q; ++b) {
                out[a * Q + b] = std::complex<double>(0, 0);
                for (int alpha = 0; alpha < P; ++alpha) {
                    out[a * Q + b] += in[alpha * Q + b] * std::polar(1.0, -2 * M_PI * a * alpha / P);
                }
                out[a * Q + b] *= std::polar(1.0, -2 * M_PI * a * b / N);
            }
        }
        if (n > 1) {
            memcpy(in, out, N * sizeof(std::complex<double>));
            for (int i = 0; i < P; ++i) {
                my_fft(in + i * Q, out + i * Q, Q, p, n - 1);
            }
        }
    }


    void my_fft1(std::complex<double> *in, std::complex<double> *out, int N) {
        int p[20] = {0};
        int n = decompose(N, p);

        auto *new_in = new std::complex<double>[N];
        auto *new_out = new std::complex<double>[N];
        memcpy(new_in, in, N * sizeof(std::complex<double>));
        my_fft(new_in, new_out, N, p, n);

        int serial[N];
        for (int i = 0; i < N; ++i) {
            int temp = i;
            serial[i] = 0;

            for (int j = 0; j < n; ++j) {
                serial[i] = serial[i] * p[j] + temp % p[j];
                temp /= p[j];
            }
        }
        for (int i = 0; i < N; ++i) {
            out[serial[i]] = new_out[i];
        }
    }

} // fft