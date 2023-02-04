/******************************************************************************
* This file contains several functions that calculate a Discrete Fourier
* transform (DFT) for a 2^3 = 8 digit sequence. The dittt() and diftt()
* functions perform the FFT, while the others are just mainly sad
* experimentations to speed up the function dftt2(), which does the basic DFT
* according to the mathematical definition of the DFT.
* By using these functions it is possible to verify the validity of
* any other DFT/FFT algorithm, since the results should be the same for
* all of the functions below. The FFT functions are restricted
* to work on only data sequences of the order 2^x, while the definition
* of DFT can be used with any amount of data.
* This file is ready to run as is without input parameters,
* the input arrays to the functions are hardcoded to every function separately. 
* Composed by: J.L. @ sometime in 2012. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRSIZ 8
#define DOUPING 6.28318530717959
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr


/******************************************************************************
* FUNCTION: dittt()
*
* PURPOSE:
*     - Calculates Fast Fourier Transform of given data series using bit
*       reversal prior to FFT. This function is directly taken from the book
*       "Numerical Recipes in C". The algorithm follows the butterfly diagram
*       of the radix-2 Decimation In Time (DIT) implementation of FFT, where
*       the bit reversal is mandatory. The function replaces the (input) array
*       data[1..2*nn] by its discrete Fourier transfor, if the parameter isign
*       is given as -1; or replaces data[1..2*nn] by nn times its inverse
*       discrete Fourier transform, if isign is given as 1. The array of data[]
*       is a complex array of length nn or, equivalently, a real array of
*       length 2*nn. nn MUST be an integer power of 2, but this is not checked!
*         
*       The needed increments of trigonometric functions are handled via 
*       a clever identity using the fact that:
*
*       exp(ia + ib) = exp(ia)*exp(ib) = exp(ia) + Z*exp(ia), 
*
*       where a is the starting angle and b is the increment. Solving for Z
*       gives:
*         
*       Z = exp(ib) - 1 = -(1 - cos(b)) + isin(b) = -2*sin^2(b/2) + isin(b).
*
*       Then the increments can be calculated with respect to the zero-angle
*       set by wr (omega-real) and wi (omega-imaginary) by relation:
*
*       (wr + wi) + Z*(wr + wi) = (wr + wi) + [-2*sin^2(b/2) + isin(b)](wr + wi).
*
*       By setting wpr (omega-phase-real) = -2*sin^2(b/2) and wpi = sin(b),
*       one has
*        
*       (wr + wi) + (wpr + wpi)*(wr + wi) = exp(ia) + Z*exp(ia),
*        
*       where the real part is: wr + wpr*wr - wpi*wi
*       and the imaginary part is: wi + wpr*wi + wpi*wr.
*       The actual incremental parts here are:
*       wpr*wr - wpi*wi for real part and wpr*wi + wpi*wr for imaginary part.  
*        
* 
* INPUT: 
*       - data[] = Array containing the data to be transformed
*                  The transformed data is stored back to this array
*                  so that the real and imaginary parts are following
*                  each other --> size of array = 2*size of data
*       - nn     = size of data = size of array/2, has to be a power of 2
*       - isign  = if -1, calculates the FFT; 
*                  if 1, calculates the IFFT i.e. the inverse.
*
* OUTPUT: - (void)
*/
void dittt()
{
double wtemp, wr, wpr, wpi, wi, theta;
double tempr, tempi;
int N = TRSIZ;
int i = 0, j = 0, n = 0, k = 0, m = 0, isign = -1,istep,mmax;
double data1[2*TRSIZ] = {0,0,1,0,4,0,9,0,2,0,3,0,4,0,5,0};
double *data;

data = &data1[0] - 1;

n = N*2;
j = 1;
// do the bit-reversal
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      SWAP(data[j], data[i]);
      SWAP(data[j+1], data[i+1]);
    }

    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j +=m;
  }

// calculate the FFT
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr*data[j] - wi*data[j+1];
        tempi = wr*data[j+1] + wi*data[j];
        data[j] = data[i] - tempr;
        data[j+1] = data[i+1] - tempi;
        data[i] =  data[i] + tempr;
        data[i+1] = data[i+1] + tempi;
        printf("\ni = %d ,j = %d, m = %d, wr = %f , wi = %f",(i-1)/2,(j-1)/2,m,wr,wi);
      }
      printf("\nm = %d ,istep = %d, mmax = %d, wr = %f , wi = %f, Z = %f",m,istep,mmax,wr,wi,atan(wi/wr)/(6.28318530717959/(1.0*n/2)));
      wtemp = wr;
      wr += wtemp*wpr - wi*wpi;
      wi += wtemp*wpi + wi*wpr;
    }
    mmax = istep;      
  }

// print the results
printf("\nFourier components from the DIT algorithm:");
for (k = 0; k < 2*N; k +=2 ) 
  printf("\n%f  %f", data[k+1], data[k+2]);
}


/******************************************************************************
* FUNCTION: diftt()
*
* PURPOSE:
*     - Calculates Fast Fourier Transform of given data series using bit
*       reversal after the FFT. This algorithm follows the butterfly diagram
*       of the radix-2 Decimation In Frequency (DIF) implementation of FFT,
*       which is modified from the previous Decimation In Time (DIT) algorithm.
*       The function replaces the (input) array
*       data[1..2*nn] by its discrete Fourier transfor, if the parameter isign
*       is given as -1; or replaces data[1..2*nn] by nn times its inverse
*       discrete Fourier transform, if isign is given as 1. The array of data[]
*       is a complex array of length nn or, equivalently, a real array of
*       length 2*nn. nn MUST be an integer power of 2, but this is not checked!
*         
*       The needed increments of trigonometric functions are handled via 
*       a clever identity using the fact that:
*
*       exp(ia + ib) = exp(ia)*exp(ib) = exp(ia) + Z*exp(ia), 
*
*       where a is the starting angle and b is the increment. Solving for Z
*       gives:
*         
*       Z = exp(ib) - 1 = -(1 - cos(b)) + isin(b) = -2*sin^2(b/2) + isin(b).
*
*       Then the increments can be calculated with respect to the zero-angle
*       set by wr (omega-real) and wi (omega-imaginary) by relation:
*
*       (wr + wi) + Z*(wr + wi) = (wr + wi) + [-2*sin^2(b/2) + isin(b)](wr + wi).
*
*       By setting wpr (omega-phase-real) = -2*sin^2(b/2) and wpi = sin(b),
*       one has
*        
*       (wr + wi) + (wpr + wpi)*(wr + wi) = exp(ia) + Z*exp(ia),
*        
*       where the real part is: wr + wpr*wr - wpi*wi
*       and the imaginary part is: wi + wpr*wi + wpi*wr.
*       The actual incremental parts here are:
*       wpr*wr - wpi*wi for real part and wpr*wi + wpi*wr for imaginary part.  
*        
* 
* INPUT: 
*       - data[] = Array containing the data to be transformed
*                  The transformed data is stored back to this array
*                  so that the real and imaginary parts are following
*                  each other --> size of array = 2*size of data
*       - nn     = size of data = size of array/2, has to be a power of 2
*       - isign  = if -1, calculates the FFT; 
*                  if 1, calculates the IFFT i.e. the inverse.
*
* OUTPUT: - (void)
*/
void diftt()
{
double wtemp, wr, wpr, wpi, wi, theta;
double tempr, tempi;
int N = TRSIZ;
int i = 0, j = 0, n = 0, k = 0, m = 0, isign = -1,istep,mmax;
double data1[2*TRSIZ] = {0,0,1,0,4,0,9,0,2,0,3,0,4,0,5,0};
double *data;            

data = &data1[0] - 1;

n = N*2;
mmax = n/2;
// calculate the FFT
  while (mmax >= 2) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = data[i];
        tempi = data[i+1];
        data[i] = data[i] + data[j];
        data[i+1] = data[i+1] + data[j+1];
        data[j] = tempr - data[j];
        data[j+1] = tempi - data[j+1];
        tempr = wr*data[j] - wi*data[j+1];
        tempi = wr*data[j+1] + wi*data[j];
        data[j] = tempr;
        data[j+1] = tempi;
        printf("\ni = %d ,j = %d, m = %d, wr = %f , wi = %f",(i-1)/2,(j-1)/2,m,wr,wi);
      }
      printf("\nm = %d ,istep = %d, mmax = %d, wr = %f , wi = %f, Z = %f",m,istep,mmax,wr,wi,atan(wi/wr)/(6.28318530717959/(1.0*n/2)));
      wtemp = wr;
      wr += wtemp*wpr - wi*wpi;
      wi += wtemp*wpi + wi*wpr;
    }
    mmax = mmax/2;      
  }

// do the bit-reversal  
j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      SWAP(data[j], data[i]);
      SWAP(data[j+1], data[i+1]);
    }

    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j +=m;
  }

// print the results
printf("\nFourier components from the DIF algorithm:");
for (k = 0; k < 2*N; k +=2 ) 
  printf("\n%f  %f", data[k+1], data[k+2]);
}

void fftt()
{
int N = TRSIZ, N_2 = TRSIZ/2;
int n = 0, k = 0;
int data[TRSIZ] = {0,1,4,9,2,3,4,5};
double ER, OR, EI, OI, KR, KI, KXR, KXI;
double reall[TRSIZ];
double imagg[TRSIZ];

// Do the transform
  for (n = 0; n < N_2 ; n++) {
    ER = OR = EI = OI = 0.0;
    KXR = cos( 6.28318530717959*n/N );
    KXI = -sin( 6.28318530717959*n/N );
    for (k = 0; k < N_2; k++) {
      KR = cos( 6.28318530717959*k*n/N_2 );
      KI = -sin( 6.28318530717959*k*n/N_2 );
      ER += KR*data[2*k];  
      OR += (KXR*KR - KXI*KI)*data[2*k + 1];
      EI += KI*data[2*k];  
      OI += (KXR*KI + KXI*KR)*data[2*k + 1];
    }
    reall[n] = ER + OR;
    imagg[n] = EI + OI;
    reall[n + N_2] = ER - OR;
    imagg[n + N_2] = EI - OI;
}

// print the results
printf("\nFourier components from another experiment");
for (k = 0; k < N; k++) 
  printf("\n%f  %f", reall[k], imagg[k]);
}

void dftt()
{
int N = TRSIZ, n = 0, k = 0;
int data[TRSIZ] = {0,1,4,9,2,3,4,5};
double reall[TRSIZ];
double imagg[TRSIZ];
double Xa, Xb, K;
double wtemp, wr, wpr, wpi, wi;

// Do the transform
  for (k = 0; k < N; k++) {
    Xa = Xb = 0.0;
    K = 6.28318530717959*k/N;
    wtemp = sin(0.5*K);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(K);
    wr = 1.0;
    wi = 0.0;
    for (n = 0; n < N; n++) {
      Xa += data[n] * wr;  
      Xb -= data[n] * wi;
      wtemp = wr;
      wr += wtemp*wpr - wi*wpi;
      wi += wtemp*wpi + wi*wpr;
    }
    reall[k] = Xa;
    imagg[k] = Xb;
  }

// print the results
printf("\nFourier components from an experiment");
for (k = 0; k < N; k++) 
  printf("\n%f  %f", reall[k], imagg[k]);
}

void dftt2()
{
int N = TRSIZ, n = 0, k = 0;
int data[TRSIZ] = {0,1,4,9,2,3,4,5};
double reall[TRSIZ];
double imagg[TRSIZ];
double Xa, Xb, K;

// Do the transform
  for (k = 0; k < N; k++) {
    Xa = Xb = 0.0;
    for (n = 0; n < N; n++) {
      K = 6.28318530717959*k*n/N;
      Xa += data[n] * cos(K);  
      Xb -= data[n] * sin(K);
    }
    reall[k] = Xa;
    imagg[k] = Xb;
  }

// print the results
printf("\nFourier components from the definition");
for (k = 0; k < N; k++) 
  printf("\n%f  %f", reall[k], imagg[k]);
}

// launch a test for comparing the results
// from each transform algorithm
int main(void)
{
// experiment
dftt();
printf("\n\n");
// experiment2
fftt();
printf("\n\n");
// definition
dftt2();
printf("\n\n");
// Decimation In Time (DIT)
dittt();
printf("\n\n");
// Decimation In Frequency (DIF)
diftt();
printf("\n\n");
return 0;
}
