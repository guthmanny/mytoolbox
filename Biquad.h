//
// Biquad.h
// Created by meng on 9/9/22.

#ifndef _Biquad_h_
#define _Biquad_h_
#include <vector>

namespace SndLab {
	
enum class BiquadAlgorithm {
	LPF,
	HPF,
	BPF,
	Notch,
	Allpass
};

struct BiquadParameters{
	BiquadParameters(){};
	BiquadAlgorithm algorithm = BiquadAlgorithm::LPF;
	double fc = 1000.0;
	double Q = 0.707;
	double dBgain = 0.0;
};
	
class Biquad
{
public:
	Biquad();
	BiquadParameters getParameters();
	void setParameters(BiquadParameters _parameters);
	bool reset(double sampleRate);
	bool calculateFilterCoeffs();
	std::vector<double> getCoeffs();
	bool processAudioSample();
private:
	BiquadParameters params;
	double sampleRate = 48000.0;
	std::vector <double> coeffArray{0.0,0.0,0.0,0.0,0.0,0.0};
	bool q_is_bandwidth = 0;
};
	
}// Sndlab namespace
#endif