#include "fft.h"
#include <cassert>
#include <iostream>
#include <algorithm>
#include <iterator>

FFT_engine::FFT_engine(const int dim): fft_dim(dim)
{
    assert(dim%2 == 0);

    fft_dim2 = (dim >> 1) + 1;

    in_array = (double*) fftw_malloc(sizeof(double) * 2*dim);
    out_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (dim + 2));
    plan_to_fft = fftw_plan_dft_r2c_1d(2*dim, in_array, out_array,  FFTW_PATIENT);
    plan_from_fft = fftw_plan_dft_c2r_1d(2*dim, out_array, in_array,  FFTW_PATIENT);

    pos_powers = vector<FFTPoly>(dim,FFTPoly(fft_dim2));
    neg_powers = vector<FFTPoly>(dim,FFTPoly(fft_dim2));
    for(int i = 0; i < dim; i++)
    {
        ModQPoly x_power(dim,0);
        x_power[i] += 1;
        FFTPoly x_power_fft(fft_dim2);
        to_fft(x_power_fft, x_power);
        pos_powers[i] = x_power_fft;

        x_power[i] -= 2;
        to_fft(x_power_fft, x_power);
        neg_powers[i] = x_power_fft;
    }
}

void FFT_engine::to_fft(FFTPoly& out, const ModQPoly& in) const
{
    assert(out.size() == fft_dim2);

    double* in_arr = in_array;
    fftw_complex* out_arr = out_array;
    int N = fft_dim;

    for (int i = 0; i < N; ++i)
    {
        in_arr[i] = double(in[i]);
        in_arr[i+N] = 0.0;
    }
    fftw_execute(plan_to_fft);
    int tmp = 1;
    for (auto it = out.begin(); it < out.end(); ++it)
    {
        fftw_complex& out_z = out_arr[tmp];
        complex<double>& outi = *it;
        outi.real(out_z[0]);
        outi.imag(out_z[1]);
        tmp += 2;
    }
}

void FFT_engine::from_fft(vector<long>& out, const FFTPoly& in) const
{
    int tmp = 0;
    double* in_arr = in_array;
    fftw_complex* out_arr = out_array;
    int N = fft_dim;
    int Nd = double(N);

    for (auto it = in.begin(); it < in.end(); ++it) 
    {
        //std::cout << "i: " << i << ", number: " << in[i] << std::endl;
        out_arr[tmp+1][0] = real(*it)/Nd;
        out_arr[tmp+1][1] = imag(*it)/Nd;
        out_arr[tmp][0] = 0.0;
        out_arr[tmp][1] = 0.0;
        tmp += 2;
    }
    fftw_execute(plan_from_fft);
    out.resize(fft_dim); 
    for (int i = 0; i < N; ++i)
    {	
        out[i] = long(round(in_arr[i]));
        //std::cout << "i: " << i << ", number: " << out[i] << std::endl;
    }
}

FFT_engine::~FFT_engine()
{
    fftw_destroy_plan(plan_to_fft);
    fftw_destroy_plan(plan_from_fft);
    fftw_free(in_array);
    fftw_free(out_array);
}

FFTPoly operator +(const FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    FFTPoly res(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i]+b[i];

    return res;
}

void operator +=(FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); ++i)
        a[i]+=b[i];
}

void operator +=(FFTPoly& a, const complex<double> b)
{
    for (size_t i = 0; i < a.size(); ++i)
        a[i]+=b;
}

FFTPoly operator -(const FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    FFTPoly res(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i]-b[i];

    return res;
}

void operator -=(FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); ++i)
        a[i]-=b[i];
}

FFTPoly operator *(const FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    FFTPoly res(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i]*b[i];

    return res;
}

void operator *=(FFTPoly& a, const FFTPoly& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); ++i)
        a[i]*=b[i];
}

FFTPoly operator *(const FFTPoly& a, const int b)
{
    FFTPoly res(a.size());
    double bd = double(b);
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i] * bd;

    return res;
}