//
// genSine.c
// Created by meng on 1/31/23.

#include "genSine.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double* genSine(double freq,double fs,double lensInSeconds)
{
	int sampleNums = (int)(lensInSeconds*fs);
	double* result =(double*) malloc(sizeof(double)*sampleNums);
	for (int n = 0; n < sampleNums; n++) {
		result[n] = sin(2*M_PI*(freq/fs)*n);
	}
	return result;
}

