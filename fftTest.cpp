#include <iostream>

#include "/Users/meng/Documents/mycode/mytoolbox/genSine.h"
#include "/Users/meng/Documents/mycode/mytoolbox/fft.h"

using namespace std;

void myplot(const char* filename);

typedef struct{
	double sineFreq;
	double fs;
	double lensInSeconds;
} sineCoe;

int main(int argc, char *argv[]) {
	
	double* sinRes;
	sineCoe sine1;
	sine1.sineFreq = 441.0;
	sine1.fs = 48000.0;
	sine1.lensInSeconds = 0.1;
	
	int samplesNums = (int)(sine1.fs * sine1.lensInSeconds);
	
	auto* in  = new std::complex<double>[samplesNums];
	auto* out = new std::complex<double>[samplesNums];
	double* magdB = new double[samplesNums];
	double* freq  = new double[samplesNums];
	

	// generating sine wave
	sinRes = genSine(sine1.sineFreq	, sine1.fs, sine1.lensInSeconds);
	for (int i = 0; i<samplesNums; i++) {
		in[i].real(sinRes[i]);
	}
	// Open file for writing
	FILE* fp;
	int result = 0;
	fp = fopen("test.dat", "w");
	if (fp == NULL) {
		result = 1;
		fprintf(stderr, "outputFFT: Could open output file for writing\n");
		exit(1);
	}

	//fft;
	fft::my_fft1(in, out, samplesNums);
	for (int i=0 ; i < samplesNums; i++)
		{
			freq[i] = sine1.fs * i / samplesNums;
			magdB[i] = 20 * log10(abs(out[i])/(samplesNums/2));
			//std::cout<<freq[i]<<" "<<magdB[i]<<endl;
			fprintf(fp, "%lf %lf\n",freq[i],magdB[i]);
		}
	if (fp != NULL) fclose(fp);	
	myplot("test.dat");

	delete[] in;
	delete[] out;
	delete[] magdB;
	delete[] freq;
	free(sinRes);

	return 0;
}

void myplot(const char* filename)
{
	FILE* command = popen("gnuplot","w");
	fprintf(command,"set boxwidth 0.5 relative\n");
	fprintf(command,"set logscale x\n");
	fprintf(command,"plot '%s' using 1:2 w l\n",filename);
	fflush(command);
	getchar();
	pclose(command);
}