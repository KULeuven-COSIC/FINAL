#include "ntruhe.h"
#include "sampler.h"
#include "fft.h"
#include "lwehe.h"

#include "time.h"
#include <cassert>

Ctxt_NTRU::Ctxt_NTRU(const Ctxt_NTRU& ct)
{
    data = ct.data;
}

Ctxt_NTRU& Ctxt_NTRU::operator=(const Ctxt_NTRU& ct)
{
    data = ct.data;
    return *this;   
}

Ctxt_NTRU Ctxt_NTRU::operator +(const Ctxt_NTRU& ct) const
{
    Ctxt_NTRU res;
    for (size_t i = 0; i < parNTRU.n; ++i)
        res.data[i] = parNTRU.mod_q_base(data[i] + ct.data[i]);

    return res;
}

Ctxt_NTRU Ctxt_NTRU::operator -(const Ctxt_NTRU& ct) const
{
    Ctxt_NTRU res;
    for (size_t i = 0; i < parNTRU.n; ++i)
        res.data[i] = parNTRU.mod_q_base(data[i] - ct.data[i]);
    return res;
}

void Ctxt_NTRU::operator -=(const Ctxt_NTRU& ct)
{
    for (size_t i = 0; i < parNTRU.n; ++i)
    {
        data[i] -= ct.data[i];
        data[i] = parNTRU.mod_q_base(data[i]);
    }     
}

void SchemeNTRU::encrypt(Ctxt_NTRU& ct, const int b) const
{
    clock_t start = clock();

    int n = parNTRU.n;

    vector<int> g(n,0L);
    Sampler::get_ternary_vector(g);
    g[0] += b * parNTRU.delta_base;
    ct.data = vector<int>(n,0);
    vector<long> ct_long(n,0L);
    for (int i = 0; i < n; i++)
    {
        long g_coef = long(g[i]);
        const vector<int>& sk_row = sk_base.sk_inv[i];
        for (int j = 0; j < n; j++)
            ct_long[j] += long(sk_row[j]) * g_coef;
    }
    parNTRU.mod_q_base(ct.data, ct_long);

    //cout << "Encryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

int SchemeNTRU::decrypt(const Ctxt_NTRU& ct) const
{
    clock_t start = clock();

    int output = 0;
    for (int i = 0; i < parNTRU.n; i++)
    {
        output += ct.data[i] * sk_base.sk[i][0];
    }
    output = parNTRU.mod_q_base(output);
    output = int(round(double(output)/double(parNTRU.q_base)*Param::t));
    //cout << "Decryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    return output;
}

void SchemeNTRU::key_switch(Ctxt_NTRU& ct, const ModQPoly& poly) const
{
    int N = Param::N;
    int B_ksk = Param::B_ksk;
    int Nl = parNTRU.Nl;
    int l_ksk = parNTRU.l_ksk;
    int n = parNTRU.n;

    vector<int> poly_decomp(Nl); 
    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    int bound = B_ksk >> 1;
    int il = 0;
    int digit, tmp, sgn;
    for (int i = 0; i < N; ++i)
    {
        tmp = poly_abs[i];
        sgn = poly_sign[i];
        for (int j = 0; j < l_ksk; j++)
        {
            int digit = tmp % B_ksk;
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
    vector<long> ct_long(n);
    for (int i = 0; i < Nl; ++i)
    {
        long tmp_int = long(poly_decomp[i]);
        const vector<int>& ksk_row = ksk[i];
        for (int j = 0; j < n; ++j)
        {
            ct_long[j] += long(ksk_row[j]) * tmp_int;
        }
    }
    parNTRU.mod_q_base(ct.data, ct_long);
}

void SchemeNTRU::mask_constant(Ctxt_NTRU& ct, int constant)
{
    int n = parNTRU.n;

    vector<int> g(n);
    Sampler::get_ternary_vector(g);
    g[0] += constant;
    ct.data = vector<int>(n);
    vector<long> ct_long(n);
    for (int i = 0; i < n; i++)
    {
        long g_coef = long(g[i]);
        const vector<int>& sk_row = sk_base.sk_inv[i];
        for (int j = 0; j < n; j++)
            ct_long[j] += sk_row[j] * g_coef;
    }
    parNTRU.mod_q_base(ct.data, ct_long);
}

void SchemeNTRU::bootstrap(Ctxt_NTRU& ct) const
{
    //clock_t start = clock();
    int N = Param::N;
    int N2 = Param::N2;
    int N2p1 = Param::N2p1;
    int B_bsk_size = Param::B_bsk_size;
    int half_delta_boot = Param::half_delta_boot;

    // switch to modulus 2*N
    modulo_switch_to_boot(ct);
    // initialize accumulator
    ModQPoly acc(N, half_delta_boot);
    //vector<long> acc_long(N);
    for (size_t i = 0; i < N/2; i++)
        acc[i] = -acc[i];

    //accumulator loop
    int coef_counter = 0;
    vector<int>& data = ct.data;
    int coef, neg_coef, coef_sign, neg_coef_sign, B, shift, l;
    double Bd;
    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);
    
    const BSKey_NTRU& boot_key = bk;
    for (int iBase = 0; iBase < B_bsk_size; ++iBase)
    { 
        B = parNTRU.B_bsk[iBase];
        Bd = double(B);
        shift = parNTRU.shift_bsk[iBase];
        l = parNTRU.l_bsk[iBase];
        const vector<vector<NGSFFTctxt>>& bk_coef_row = boot_key[iBase];
        vector<FFTPoly> mux_fft(l, FFTPoly(N2p1));
        for (int iCoef = 0; iCoef < parNTRU.bsk_partition[iBase]; ++iCoef)
        { 
            coef = data[iCoef+coef_counter];
            if (coef == 0) continue;
            coef_sign = 1;
            if (coef < 0) coef += N2;
            if (coef >= N)
            {
                coef -= N;
                coef_sign = -1;
            }
            neg_coef = N - coef;
            neg_coef_sign = -coef_sign;
            if(neg_coef == N)
            {
                neg_coef = 0;
                neg_coef_sign = -neg_coef_sign;
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
            
            // bk_0[i] - bk_1[i]*X^-coef
            const FFTPoly& x_power_neg = (neg_coef_sign == 1) ? fftN.pos_powers[neg_coef]: fftN.neg_powers[neg_coef];
            const NGSFFTctxt& bk_part_row0 = bk_coef_row[iCoef][0];
            const NGSFFTctxt& bk_part_row1 = bk_coef_row[iCoef][1];
            for (int iPart = 0; iPart < l; ++iPart)
            {
                FFTPoly tmp_fft = bk_part_row1[iPart];
                tmp_fft *= x_power_neg;
                mux_fft[iPart] = bk_part_row0[iPart];
                mux_fft[iPart] -= tmp_fft;
            }
            // acc * (X^coef - 1) x (bk_0[i] - bk_1[i]*X^-coef)
            external_product(tmp_poly_long, tmp_poly, mux_fft, B, shift, l);
            mod_q_boot(tmp_poly, tmp_poly_long);
            // acc * (X^coef - 1) x (bk_0[i] - bk_1[i]*X^-coef) + acc
            for (int i = 0; i<N; ++i)
                acc[i] += tmp_poly[i];
        }
        coef_counter += parNTRU.bsk_partition[iBase];
    }
    
    // add floor(q_boot/(2*t)) to all coefficients of the accumulator
    for (auto it = acc.begin(); it < acc.end(); ++it)
        *it += half_delta_boot;

    //mod q_boot of the accumulator
    mod_q_boot(acc);

    //mod switch to q_base
    modulo_switch_to_base_ntru(acc);

    //key switch
    key_switch(ct, acc);

    //cout << "Bootstrapping: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeNTRU::nand_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const
{
    //clock_t start = clock();
    ct_res = ct_nand_const - ct1 - ct2;
    bootstrap(ct_res);
    //cout << "NAND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeNTRU::and_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const
{
    //clock_t start = clock();
    ct_res = ct_and_const - ct1 - ct2;
    bootstrap(ct_res);
    //cout << "AND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeNTRU::or_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const
{
    //clock_t start = clock();
    ct_res = ct_or_const - ct1 - ct2;
    bootstrap(ct_res);
    //cout << "OR: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void SchemeNTRU::xor_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const
{
    if (ct_res.data.size() != parNTRU.n)
        ct_res.data = vector<int>(parNTRU.n);

    ct_res = ct1 + ct2;
    ct_res = ct_res + ct_res; // ct_res = 2*(ct1 + ct2) % q

    bootstrap(ct_res);
}

void SchemeNTRU::not_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct) const
{
    ct_res = ct_not_const - ct;
}

