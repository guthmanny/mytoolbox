//
// Biquad.cpp
// Created by meng on 9/9/22.

#include "Biquad.h"
#include <cmath>
namespace SndLab {
	
Biquad::Biquad(){
}
std::vector<double> Biquad::getCoeffs(){
	return coeffArray;
}
bool Biquad::calculateFilterCoeffs(){
	
	double alpha, a0,a1,a2,b0,b1,b2;
	const double temp_pi=3.1415926535897932384626433832795;
	double omega,tsin,tcos;
	
	switch (params.algorithm) {
		
		case BiquadAlgorithm::LPF:
			
			omega	=	2.0*temp_pi*params.fc/sampleRate;
			tsin	=	sin(omega);
			tcos	=	cos(omega);			
			if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*params.Q*omega/tsin);
			else
				alpha=tsin/(2.0*params.Q);
			
			b0=(1.0-tcos)/2.0;
			b1=1.0-tcos;
			b2=(1.0-tcos)/2.0;
			a0=1.0+alpha;
			a1=-2.0*tcos;
			a2=1.0-alpha;			
			
			break;
		
		case BiquadAlgorithm::HPF:
			omega	=	2.0*temp_pi*params.fc/sampleRate;
			tsin	=	sin(omega);
			tcos	=	cos(omega);			
			if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*params.Q*omega/tsin);
			else
				alpha=tsin/(2.0*params.Q);
			
			b0=(1.0+tcos)/2.0;
			b1=-(1.0+tcos);
			b2=(1.0+tcos)/2.0;
			a0=1.0+ alpha;
			a1=-2.0*tcos;
			a2=1.0-alpha;	
			
			break;
		
		case BiquadAlgorithm::BPF:
			omega	=	2.0*temp_pi*params.fc/sampleRate;
			tsin	=	sin(omega);
			tcos	=	cos(omega);			
			if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*params.Q*omega/tsin);
			else
				alpha=tsin/(2.0*params.Q);
			
			b0=tsin/2.0;
			b1=0.0;
			b2=-tsin/2;
			a0=1.0+alpha;
			a1=-2.0*tcos;
			a2=1.0-alpha;
			
			break;
		
		case BiquadAlgorithm::Notch:
			
			omega	=	2.0*temp_pi*params.fc/sampleRate;
			tsin	=	sin(omega);
			tcos	=	cos(omega);			
			if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*params.Q*omega/tsin);
			else
				alpha=tsin/(2.0*params.Q);
			
			b0=1.0;
			b1=-2.0*tcos;
			b2=1.0;
			a0=1.0+alpha;
			a1=-2.0*tcos;
			a2=1.0-alpha;

			break;
		
		case BiquadAlgorithm::Allpass:
			
			omega	=	2.0*temp_pi*params.fc/sampleRate;
			tsin	=	sin(omega);
			tcos	=	cos(omega);			
			if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*params.Q*omega/tsin);
			else
				alpha=tsin/(2.0*params.Q);

			b0=1.0-alpha;
			b1=-2.0*tcos;
			b2=1.0+alpha;
			a0=1.0+alpha;
			a1=-2.0*tcos;
			a2=1.0-alpha;
			
			break;
		
	}
	coeffArray[0] = b0/a0;
	coeffArray[1] = b1/a0;
	coeffArray[2] = b2/a0;
	coeffArray[3] = 1.0;
	coeffArray[4] = a1/a0;
	coeffArray[5] = a2/a0;
	
	return true;
}
	
}// SndLab namespace

