#include "lwehe.h"
#include "sampler.h"
#include "fft.h"
#include "ntruhe.h"

#include <chrono>
#include <vector>
#include <cassert>

using namespace std;

Ctxt_LWE::Ctxt_LWE(const Ctxt_LWE& ct)
{
    a = ct.a;
    b = ct.b;
}

Ctxt_LWE& Ctxt_LWE::operator=(const Ctxt_LWE& ct)
{
    a = ct.a;
    b = ct.b;
    return *this;   
}

Ctxt_LWE Ctxt_LWE::operator +(const Ctxt_LWE& ct) const
{
    Ctxt_LWE res;
    for (size_t i = 0; i < parLWE.n; i++)
        res.a[i] = parLWE.mod_q_base(a[i] + ct.a[i]);

    res.b = parLWE.mod_q_base(b + ct.b);

    return res;
}

Ctxt_LWE Ctxt_LWE::operator -(const Ctxt_LWE& ct) const
{
    Ctxt_LWE res;
    for (size_t i = 0; i < parLWE.n; i++)
        res.a[i] = parLWE.mod_q_base(a[i] - ct.a[i]);

    res.b = parLWE.mod_q_base(b + ct.b);

    return res;
}

Ctxt_LWE operator -(const int c, const Ctxt_LWE& ct)
{
    Ctxt_LWE res;
    res.a = vector<int>(parLWE.n);
    const vector<int>& a = ct.a;
    for (size_t i = 0; i < a.size(); i++)
        res.a[i] = parLWE.mod_q_base(-a[i]);

    res.b = parLWE.mod_q_base(c-ct.b);
    return res;
}

void SchemeLWE::encrypt(Ctxt_LWE& ct, int m) const
{
    clock_t start = clock();

    int n = parLWE.n;

    vector<int> a(n,0L);
    Sampler s(parLWE);
    s.get_uniform_vector(a);
    ct.a = a;
    normal_distribution<double> gaussian_sampler(0.0, Param::e_st_dev);
    int b = parLWE.delta_base*m + static_cast<int>(round(gaussian_sampler(rand_engine)));
    for (int i = 0; i < n; i++)
    {
        b -= sk_base[i] * a[i];
    }
    parLWE.mod_q_base(b);
    ct.b = b;

    //cout << "Encryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

int SchemeLWE::decrypt(const Ctxt_LWE& ct) const
{
    clock_t start = clock();

    int output = ct.b;
    for (int i = 0; i < parLWE.n; i++)
    {
        output += ct.a[i] * sk_base[i];
    }
    output = parLWE.mod_q_base(output);
    output = int(round(double(output*Param::t)/double(parLWE.q_base)));
    //cout << "Decryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    return output;
}

void external_product(vector<long>& res, const vector<int>& poly, const vector<FFTPoly>& poly_vector, int b, int shift, int l)
{ 
    int N = Param::N;
    int N2p1 = Param::N2p1;

    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    vector<int> poly_decomp(N);

    for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    FFTPoly res_fft(N2p1);
    FFTPoly tmp_fft(N2p1);
    int mask = b-1;
    int bound = b >> 1;
    int digit, sgn;
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int& abs_val = poly_abs[i];
            digit = abs_val & mask;
            if (digit > bound)
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? (digit - b): (b - digit);
                abs_val >>= shift;
                ++abs_val;
            }
            else
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? digit: -digit;
                abs_val >>= shift;
            }
        }
        fftN.to_fft(tmp_fft, poly_decomp);
        tmp_fft *= poly_vector[j];
        res_fft += tmp_fft;
    }
    fftN.from_fft(res, res_fft);
}

void SchemeLWE::key_switch(Ctxt_LWE& ct, const ModQPoly& poly) const
{
    int N = Param::N;
    int B_ksk = Param::B_ksk;
    int l_ksk = parLWE.l_ksk;
    int Nl = parLWE.Nl;
    int n = parLWE.n;

    vector<int> poly_decomp(Nl); 
    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    int digit;
    int il = 0;
    int tmp;
    int sgn;
    int bound = B_ksk >> 1;
    for (int i = 0; i < N; ++i)
    {
        tmp = poly_abs[i];
        sgn = poly_sign[i];
        for (int j = 0; j < l_ksk; ++j)
        {
            digit = tmp % B_ksk;
            if (digit > bound)
            {
                poly_decomp[il+j] = (sgn == 1) ? (digit - B_ksk): (B_ksk - digit);
                tmp /= B_ksk;
                ++tmp;
            }
            else
            {
                poly_decomp[il+j] = (sgn == 1) ? digit:  - digit;
                tmp /= B_ksk;
            }
        }
        il += l_ksk;
    }
    vector<long> a(n);
    for (int i = 0; i < Nl; ++i)
    {
        long tmp_int = long(poly_decomp[i]);
        const vector<int>& ksk_row = ksk.A[i];
        for (int j = 0; j < n; ++j)
        {
            a[j] += long(ksk_row[j]) * tmp_int;
        }
    }
    parLWE.mod_q_base(ct.a, a);
    long b = 0L;
    const vector<int>& ksk_b = ksk.b;
    for (int i = 0; i < Nl; ++i)
        b += ksk_b[i] * long(poly_decomp[i]);
    ct.b = parLWE.mod_q_base(b);
}

void SchemeLWE::bootstrap(Ctxt_LWE& ct) const
{
    //clock_t start = clock();
    int N = Param::N;
    int N2 = Param::N2;
    int N2p1 = Param::N2p1;
    int B_bsk_size = Param::B_bsk_size;
    int half_delta_boot = Param::half_delta_boot;

    // switch to modulus 2*N
    modulo_switch_to_boot(ct);
    // initialize accumulator and rotate accumulator by X^ct.b
    vector<int> acc(N, half_delta_boot);
    int b_pow = (N/2 + ct.b)%N2;
    if (b_pow < 0)
        b_pow += N2;
    int b_sign = 1;
    if (b_pow >= N)
    {
        b_pow -= N;
        b_sign = -1;
    }
    for (int i = 0; i < b_pow; ++i)
        acc[i] = (b_sign == 1) ? -acc[i]: acc[i];
    for (int i = b_pow; i < N; ++i)
        acc[i] = (b_sign == 1) ? acc[i]: -acc[i];

    //accumulator loop
    int coef_counter = 0;
    vector<int>& a = ct.a;
    int coef, coef_sign, B, shift, l;
    double Bd;
    //auto start = clock();
    //float cmux_time = 0.0;
    //float extprod_time = 0.0;
    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);

    const BSKey_LWE& boot_key = bk;
    for (int iBase = 0; iBase < B_bsk_size; ++iBase)
    {
        B = parLWE.B_bsk[iBase];
        Bd = double(B);
        shift = parLWE.shift_bsk[iBase];
        l = parLWE.l_bsk[iBase];
        const vector<NGSFFTctxt>& bk_coef_row = boot_key[iBase];
        for (int iCoef = 0; iCoef < parLWE.bsk_partition[iBase]; ++iCoef)
        { 
            //auto start = clock();
            coef = a[iCoef+coef_counter];
            if (coef == 0) continue;
            coef_sign = 1;
            if (coef < 0) coef += N2;
            if (coef >= N)
            {
                coef -= N;
                coef_sign = -1;
            }

            // acc * (X^coef - 1)
            if (coef_sign == 1)
            {
                for (int i = 0; i<coef; ++i)
                    tmp_poly[i] = mod_q_boot(-acc[i-coef+N] - acc[i]);
                for (int i = coef; i < N; ++i)
                    tmp_poly[i] = mod_q_boot(acc[i-coef] - acc[i]);
            }
            else
            {
                for (int i = 0; i<coef; ++i)
                    tmp_poly[i] = mod_q_boot(acc[i-coef+N] - acc[i]);
                for (int i = coef; i < N; ++i)
                    tmp_poly[i] = mod_q_boot(-acc[i-coef] - acc[i]);
            }
            //cmux_time += float(clock()-start)/CLOCKS_PER_SEC;    
  
            //start = clock();
            // acc * (X^coef - 1) x bk[i]
            external_product(tmp_poly_long, tmp_poly, bk_coef_row[iCoef], B, shift, l);
            mod_q_boot(tmp_poly, tmp_poly_long);
            // acc * (X^coef - 1) x bk[i] + acc
            for (int i = 0; i<N; ++i)
                acc[i] += tmp_poly[i];
            //extprod_time += float(clock()-start)/CLOCKS_PER_SEC;
        }
        coef_counter += parLWE.bsk_partition[iBase];
    }
    //cout << "Cmux: " << cmux_time << endl;
    //cout << "Ext. prod: " << extprod_time << endl;

    // add floor(q_boot/(2*t)) to all coefficients of the accumulator
    for (auto it = acc.begin(); it < acc.end(); ++it)
        *it += half_delta_boot;

    //mod q_boot of the accumulator
    mod_q_boot(acc);
    
    //mod switch to q_base
    modulo_switch_to_base_lwe(acc);
    
    //key switch
    //auto start = clock();
    key_switch(ct, acc);
    //cout << "Key-switching: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    //cout << "Bootstrapping: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeLWE::bootstrap2(Ctxt_LWE& ct) const
{
    //clock_t start = clock();
    int N = Param::N;
    int N2 = Param::N2;
    int N2p1 = Param::N2p1;
    int B_bsk_size = Param::B_bsk_size;
    int half_delta_boot = Param::half_delta_boot;

    // switch to modulus 2*N
    modulo_switch_to_boot(ct);
    // initialize accumulator and rotate accumulator by X^ct.b
    vector<int> acc(N, half_delta_boot);
    int b_pow = (N/2 + ct.b)%N2;
    if (b_pow < 0)
        b_pow += N2;
    int b_sign = 1;
    if (b_pow >= N)
    {
        b_pow -= N;
        b_sign = -1;
    }
    for (int i = 0; i < b_pow; ++i)
        acc[i] = (b_sign == 1) ? -acc[i]: acc[i];
    for (int i = b_pow; i < N; ++i)
        acc[i] = (b_sign == 1) ? acc[i]: -acc[i];

    //accumulator loop
    int coef_counter = 0;
    vector<int>& a = ct.a;
    int coef1, coef2, coef_sum, coef1_sign, coef2_sign, coef_sum_sign, B, shift, l;
    double Bd;
    //auto start = clock();
    //float cmux_time = 0.0;
    //float extprod_time = 0.0;
    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);

    const BSKey_LWE& boot_key = bk;
    for (int iBase = 0; iBase < B_bsk_size; ++iBase)
    {
        B = parLWE.B_bsk[iBase];
        Bd = double(B);
        shift = parLWE.shift_bsk[iBase];
        l = parLWE.l_bsk[iBase];
        const vector<NGSFFTctxt>& bk_coef_row = boot_key[iBase];
        vector<FFTPoly> mux_fft(l,FFTPoly(N2p1,complex<double>(0.0,0.0)));
        for (int iCoef = 0; iCoef < parLWE.bsk_partition[iBase]; iCoef+=2)
        { 
            //auto start = clock();
            // normalize coef1
            coef1 = a[iCoef+coef_counter];
            coef1_sign = 1;
            if (coef1 < 0) coef1 += N2;
            if (coef1 >= N)
            {
                coef1 -= N;
                coef1_sign = -1;
            }
            // normalize coef2
            coef2 = a[iCoef+coef_counter+1];
            coef2_sign = 1;
            if (coef2 < 0) coef2 += N2;
            if (coef2 >= N)
            {
                coef2 -= N;
                coef2_sign = -1;
            }
            // normalize coef_sum
            coef_sum = (a[iCoef+coef_counter] + a[iCoef+coef_counter+1]) % N2;
            coef_sum_sign = 1;
            if (coef_sum < 0) coef_sum += N2;
            if (coef_sum >= N)
            {
                coef_sum -= N;
                coef_sum_sign = -1;
            }

            // bk_0 * X^(c_0+c_1) + bk_1 * X^(c_0) + bk_2 * X^(c_1) + bk_3
            const FFTPoly& x_sum = (coef_sum_sign == 1) ? fftN.pos_powers[coef_sum]: fftN.neg_powers[coef_sum];
            const FFTPoly& x_c1 = (coef1_sign == 1) ? fftN.pos_powers[coef1]: fftN.neg_powers[coef1];
            const FFTPoly& x_c2 = (coef2_sign == 1) ? fftN.pos_powers[coef2]: fftN.neg_powers[coef2];
            const NGSFFTctxt& bk_part_row0 = bk_coef_row[4*(iCoef >> 1)];
            const NGSFFTctxt& bk_part_row1 = bk_coef_row[4*(iCoef >> 1)+1];
            const NGSFFTctxt& bk_part_row2 = bk_coef_row[4*(iCoef >> 1)+2];
            const NGSFFTctxt& bk_part_row3 = bk_coef_row[4*(iCoef >> 1)+3];
            for (int iPart = 0; iPart < l; ++iPart)
            {
                mux_fft[iPart] = bk_part_row0[iPart];
                mux_fft[iPart] *= x_sum;

                if (coef1 == coef2 && coef1_sign == coef2_sign)
                {
                    FFTPoly tmp_fft = bk_part_row1[iPart];
                    tmp_fft += bk_part_row2[iPart];
                    tmp_fft *= x_c1;
                    mux_fft[iPart] += tmp_fft;
                }
                else
                {
                    FFTPoly tmp_fft = bk_part_row1[iPart];
                    tmp_fft *= x_c1;
                    mux_fft[iPart] += tmp_fft;

                    tmp_fft = bk_part_row2[iPart];
                    tmp_fft *= x_c2;
                    mux_fft[iPart] += tmp_fft;
                }

                mux_fft[iPart] += bk_part_row3[iPart];
            }
            //cmux_time += float(clock()-start)/CLOCKS_PER_SEC;    
  
            //start = clock();
            external_product(tmp_poly_long, acc, mux_fft, B, shift, l);
            mod_q_boot(acc, tmp_poly_long);
            //extprod_time += float(clock()-start)/CLOCKS_PER_SEC;
        }
        coef_counter += parLWE.bsk_partition[iBase];
    }
    //cout << "Cmux: " << cmux_time << endl;
    //cout << "Ext. prod: " << extprod_time << endl;

    // add floor(q_boot/(2*t)) to all coefficients of the accumulator
    for (auto it = acc.begin(); it < acc.end(); ++it)
        *it += half_delta_boot;

    //mod q_boot of the accumulator
    mod_q_boot(acc);
    
    //mod switch to q_base
    modulo_switch_to_base_lwe(acc);
    
    //key switch
    //auto start = clock();
    key_switch(ct, acc);
    //cout << "Key-switching: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    //cout << "Bootstrapping: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeLWE::nand_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const
{
    //clock_t start = clock();
    ct_res = parLWE.nand_const - (ct1 + ct2);
    bootstrap(ct_res);
    //cout << "NAND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeLWE::and_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const
{
    //clock_t start = clock();
    ct_res = parLWE.and_const - (ct1 + ct2);
    bootstrap(ct_res);
    //cout << "AND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeLWE::or_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const
{
    //clock_t start = clock();
    ct_res = parLWE.or_const - (ct1 + ct2);
    bootstrap(ct_res);
    //cout << "OR: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}


void SchemeLWE::xor_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const
{
    if (ct_res.a.size() != parLWE.n)
        ct_res.a = vector<int>(parLWE.n);

    ct_res = ct1 + ct2;
    ct_res = ct_res + ct_res; // ct = 2*(ct1 + ct2)

    bootstrap(ct_res);
}

void SchemeLWE::not_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct) const
{
    if (ct_res.a.size() != parLWE.n)
        ct_res.a = vector<int>(parLWE.n);

    ct_res = parLWE.delta_base - ct;
}

