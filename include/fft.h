#ifndef FFT
#define FFT

#include <vector>
#include <fftw3.h>
#include <complex>
#include <map>
#include "params.h"

using namespace std;

class FFT_engine
{
    int fft_dim;
    int fft_dim2;

    fftw_plan plan_to_fft;
    fftw_plan plan_from_fft;

    double* in_array;
    fftw_complex* out_array;

public:
    vector<FFTPoly> pos_powers;
    vector<FFTPoly> neg_powers;

    FFT_engine() = delete;
    FFT_engine(const int dim);

    void to_fft(FFTPoly& out, const ModQPoly& in) const;
    void from_fft(vector<long>& out, const FFTPoly& in) const;

    ~FFT_engine();
};

FFTPoly operator *(const FFTPoly& a, const FFTPoly& b);
void operator *=(FFTPoly& a, const FFTPoly& b);
FFTPoly operator *(const FFTPoly& a, const int b);
FFTPoly operator +(const FFTPoly& a, const FFTPoly& b);
void operator +=(FFTPoly& a, const FFTPoly& b);
void operator +=(FFTPoly& a, const complex<double> b);
FFTPoly operator -(const FFTPoly& a, const FFTPoly& b);
void operator -=(FFTPoly& a, const FFTPoly& b);

// global FFT engine of dimension N
const FFT_engine fftN(Param::N);


#endif