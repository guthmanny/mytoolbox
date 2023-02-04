#ifndef __MYTOOLBOX__H
#define __MYTOOLBOX__H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc/malloc.h>
//#include <fftw3.h>
#include <complex.h>
//#include "fft.h"

typedef struct filter
{
	double b0;
	double b1;
	double b2;
	double a1;
	double a2;
}filter;

extern double* linspace(double d1, double d2, int nlength);
extern double* logspace(double d1, double d2, int nlength);
extern double* expSineSweep(double fstart, double fstop, double nlength, double fs);
extern double* freqz(int freqPoint, double* dBMag, filter* coe, double fstart, double fstop, double fs);
//extern void freqLogResp(char* filename, double fstart, double fstop, double* in, int inLens, double fs);

#endif /* __MYTOOLBOX__H */