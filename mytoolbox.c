#include "mytoolbox.h"

double* linspace(double d1, double d2, int nlength)
{
	int n1 = nlength - 1;
	double* y = (double*)calloc(nlength, sizeof(double));
	if (y==NULL) {
		printf("Could not allocate memory for wd\n");
		exit(0);
	}
	for (int i = 0; i < nlength; i++) {
		y[i] = d1 + i*(d2 - d1)/n1;
	}
	return y;
}

/* logspace generates a row vector of logrithmically equally
   spaced point between d1 and d2.
*/
double* logspace(double d1, double d2, int nlength)
{
	double start = log10(d1);
	double end   = log10(d2);
	double* y = linspace(start,end,nlength);
	for (int i = 0; i < nlength; i++) {
		y[i] = pow(10, y[i]);
	}
	return y;
}

double* expSineSweep(double fstart, double fstop, double lensInSeconds, double fs)
{
	int n = 0;
	double tmp;
	double nlength = lensInSeconds*fs;
	double n1 = nlength - 1;
	double w1 = 2 * M_PI * fstart / fs;
	double w2 = 2 * M_PI * fstop / fs;
	double lgw2w1 = log(w2/w1);
	double* y = (double*)calloc(nlength,sizeof(double));
	if (y==NULL) {
		printf("Could not allocate memory for wd\n");
		exit(0);
	}
	for (int i = 0; i < nlength; i++) {
		n = n+1;
		y[i] = sin(w1*n1/lgw2w1*(exp((n/n1)*lgw2w1)-1));
	}
	return y;
}

 double* freqz(int freqPoint, double* dBMag, filter* coe, double fstart, double fstop, double fs)
{
	double* freq = logspace(fstart, fstop, freqPoint);
	
	double* wd = (double *)calloc(freqPoint, sizeof(double));
	if (wd==NULL) {
		printf("Could not allocate memory for wd\n");
		exit(0);
	}
	
	for (int i = 0; i < freqPoint; i++) {
		wd[i] = 2 * M_PI * freq[i]/fs;
	}
	
	double complex * e = (double complex *)calloc(freqPoint,sizeof(double complex));
	if(e==NULL){
		printf("Could not allocate memory for e\n");
		exit(0);
	}
	for (int i = 0; i < freqPoint; i++){
		e[i] = cexp(-I * wd [i]);
	}
	double complex * Be    = (double complex *)calloc(freqPoint,sizeof(double complex));
	double complex * Ae    = (double complex *)calloc(freqPoint,sizeof(double complex));
	double complex * H     = (double complex *)calloc(freqPoint,sizeof(double complex));
	double* mag   = (double*)calloc(freqPoint,sizeof(double));
	
	for (int i=0; i < freqPoint; i++) {
		Be[i]    = coe->b0 + e[i] * (coe->b1 + coe->b2*e[i]);
		Ae[i]    = 1 + e[i] * (coe->a1 + coe->a2*e[i]);
		H[i]     = Be[i]/Ae[i];
		mag[i]   = cabs(H[i]);
		dBMag[i] = 20*log10(mag[i]); 
	}
	return freq;
	
	free(wd);
	free(e);
	free(Be);
	free(Ae);
	free(H);
	free(mag);
}

//void freqLogResp(char* filename, double fstart, double fstop, double* input, int inLens, double fs)
//{
//	
//	double* in;
//	fftw_complex* out;
//	fftw_complex* out1;
//	fftw_plan plan;
//	int result;
//	int i;
//	FILE *fp;
//	
//	double tlens = ((double)inLens)/fs; //time length in Seconds
//	
//	// generate the reference signal
//	double * refsig = expSineSweep(fstart, fstop, tlens, fs);
//	
//	//the nearest power of two sequence length for FFT operations.
//	int fftSize = pow(2, ceil(log2(inLens)));
//	
//	// Calculate size of result data
//	int resultSize = (fftSize / 2) + 1;
//	
//	// Allocate memory to hold input and output data
//	in = (double *) fftw_malloc(fftSize * sizeof(double));
//	out = (fftw_complex *) fftw_malloc(resultSize * sizeof(fftw_complex));
//	out1 = (fftw_complex *) fftw_malloc(resultSize * sizeof(fftw_complex));
//	if (in == NULL || out == NULL || out1 == NULL) {
//		result = 1;
//		fprintf(stderr, "outputFFT: Could not allocate input/output data\n");
//		goto finalise;
//	}
//	
//	// Create the plan and check for success
//	plan = fftw_plan_dft_r2c_1d(fftSize, in, out, FFTW_MEASURE); 
//	if (plan == NULL) {
//		result = 1;
//		fprintf(stderr, "outputFFT: Could not create plan\n");
//		goto finalise;
//	}
//	
//	// Copy window and add zero padding (if required)
//	for (i=0 ; i < inLens ; i++) in[i] = input[i];
//	for ( ; i<fftSize ; i++) in[i] = 0;
//	
//	// Perform fft
//	fftw_execute(plan);
//
//
//	// Create the plan and check for success
//	plan = fftw_plan_dft_r2c_1d(fftSize, in, out1, FFTW_MEASURE); 
//	if (plan == NULL) {
//		result = 1;
//		fprintf(stderr, "outputFFT: Could not create plan\n");
//		goto finalise;
//	}
//
//	// reference signal fft
//	for (i=0 ; i < inLens ; i++) in[i] = refsig[i];
//	for ( ; i<fftSize ; i++) in[i] = 0;
//	
//	// Perform fft
//	fftw_execute(plan);
//
//	
//	// Open file for writing
//	fp = fopen(filename, "w");
//	if (fp == NULL) {
//		result = 1;
//		fprintf(stderr, "outputFFT: Could open output file for writing\n");
//		goto finalise;
//	}
//	
//	// Output result
//	for (i=0 ; i<resultSize ; i++)
//		{
//			double freq = fs * i / fftSize;
//			double mag = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
//			double magdB = 20 * log10(mag);
//			double mag1 = sqrt(out1[i][0] * out1[i][0] + out1[i][1] * out1[i][1]);
//			double magdB1 = 20 * log10(mag1);
//			fprintf(fp, "%f %f\n",freq, magdB-magdB1);
//		}
//
//	// Perform any cleaning up
//	finalise:
//	if (plan != NULL) fftw_destroy_plan(plan);
//	if (in != NULL) fftw_free(in);
//	if (out != NULL) fftw_free(out);
//	if (fp != NULL) fclose(fp);
//}

