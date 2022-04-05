#include "keygen.h"
#include "params.h"
#include "fft.h"

#include<iostream>
#include<time.h>
#include<algorithm>

using namespace NTL;
using namespace std;

void KeyGen::get_sk_boot(SKey_boot& sk_boot)
{
    cout << "Started generating the secret key of the bootstrapping scheme" << endl;
    clock_t start = clock();
    sk_boot.sk = ModQPoly(Param::N,0);
    sk_boot.sk_inv = ModQPoly(Param::N,0);

    sampler.get_invertible_vector(sk_boot.sk, sk_boot.sk_inv, Param::t, 1L);
    cout << "Generation time of the secret key of the bootstrapping scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_sk_base(SKey_base_NTRU& sk_base)
{
    cout << "Started generating the secret key of the base scheme" << endl;
    clock_t start = clock();
    sk_base.sk = ModQMatrix(param.n, vector<int>(param.n,0L));
    sk_base.sk_inv = ModQMatrix(param.n, vector<int>(param.n,0L));

    sampler.get_invertible_matrix(sk_base.sk, sk_base.sk_inv, 1L, 0L);
    cout << "Generation time of the secret key of the base scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_sk_base(SKey_base_LWE& sk_base)
{
    cout << "Started generating the secret key of the base scheme" << endl;
    clock_t start = clock();
    sk_base.clear();
    sk_base = vector<int>(param.n,0L);

    sampler.get_binary_vector(sk_base);
    cout << "Generation time of the secret key of the base scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_ksk(KSKey_NTRU& ksk, const SKey_base_NTRU& sk_base, const SKey_boot& sk_boot)
{
    //cout << "Started key-switching key generation" << endl;
    clock_t start = clock();
    // reset key-switching key
    ksk.clear();
    ksk = ModQMatrix(param.Nl, vector<int>(param.n,0));
    vector<vector<long>> ksk_long(param.Nl, vector<long>(param.n,0L));
    //cout << "Reset time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    // noise matrix G as in the paper
    ModQMatrix G(param.Nl, vector<int>(param.n,0L));
    sampler.get_ternary_matrix(G);
    //cout << "G gen time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    
    // matrix G + P * Phi(f) * E as in the paper
    int coef_w_pwr = sk_boot.sk[0];
    for (int i = 0; i < param.l_ksk; i++)
    {
        G[i][0] += coef_w_pwr;
        coef_w_pwr *= Param::B_ksk;
    }
    for (int i = 1; i < Param::N; i++)
    {
        coef_w_pwr = -sk_boot.sk[Param::N-i];
        for (int j = 0; j < param.l_ksk; j++)
        {
            G[i*param.l_ksk+j][0] += coef_w_pwr;
            coef_w_pwr *= Param::B_ksk;
        }
    }
    //cout << "G+P time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    // parameters of the block optimization of matrix multiplication
    int block = 4;
    int blocks = (param.n/block)*block;
    int rem_block = param.n%block;
    // (G + P * Phi(f) * E) * F^(-1) as in the paper
    for (int i = 0; i < param.Nl; i++)
    {
        //cout << "i: " << i << endl;
        vector<long>& k_row = ksk_long[i];
        vector<int>& g_row = G[i];
        for (int k = 0; k < param.n; k++)
        {
            const vector<int>& f_row = sk_base.sk_inv[k];
            //cout << "j: " << j << endl;
            long coef = long(g_row[k]);
            for (int j = 0; j < blocks; j+=block)
            {
                k_row[j] += (coef * f_row[j]);
                k_row[j+1] += (coef * f_row[j+1]);
                k_row[j+2] += (coef * f_row[j+2]);
                k_row[j+3] += (coef * f_row[j+3]);
            }
            for (int j = 0; j < rem_block; j++)
                k_row[blocks+j] += (coef * f_row[blocks+j]);
        }
    }
    //cout << "After K time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    // reduce modulo q_base
    for (int i = 0; i < param.Nl; i++)
        param.mod_q_base(ksk[i], ksk_long[i]);
    cout << "KSKey-gen time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_ksk(KSKey_LWE& ksk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot)
{
    //cout << "Started key-switching key generation" << endl;
    clock_t start = clock();
    // reset key-switching key
    ksk.A.clear();
    ksk.b.clear();
    for (int i = 0; i < param.Nl; i++)
    {
        vector<int> row(param.n,0L);
        ksk.A.push_back(row);
    }
    ksk.b = vector<int>(param.Nl, 0L);
    //cout << "Reset time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    // noise matrix G as in the paper
    sampler.get_uniform_matrix(ksk.A);
    //cout << "A gen time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    
    // matrix P * f_0 as in the paper
    vector<int> Pf0(param.Nl, 0L);
    int coef_w_pwr = sk_boot.sk[0];
    for (int i = 0; i < param.l_ksk; i++)
    {
        Pf0[i] += coef_w_pwr;
        coef_w_pwr *= Param::B_ksk;
    }
    for (int i = 1; i < Param::N; i++)
    {
        coef_w_pwr = -sk_boot.sk[Param::N-i];
        for (int j = 0; j < param.l_ksk; j++)
        {
            Pf0[i*param.l_ksk+j] += coef_w_pwr;
            coef_w_pwr *= Param::B_ksk;
        }
    }
    //cout << "Pf0 time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

    // A*s_base + e + Pf0 as in the paper
    normal_distribution<double> gaussian_sampler(0.0, Param::e_st_dev);
    for (int i = 0; i < param.Nl; i++)
    {
        //cout << "i: " << i << endl;
        vector<int>& k_row = ksk.A[i];
        for (int k = 0; k < param.n; k++)
            ksk.b[i] -= k_row[k] * sk_base[k];
        ksk.b[i] += (Pf0[i] + static_cast<int>(round(gaussian_sampler(rand_engine))));
        param.mod_q_base(ksk.b[i]);
    }   
    cout << "KSKey-gen time: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_bsk(BSKey_NTRU& bsk, const SKey_base_NTRU& sk_base, const SKey_boot& sk_boot)
{
    clock_t start = clock();
    
    // index of a secret key coefficient of the base scheme
    int coef_counter = 0;

    // reset the input
    bsk.clear();

    // loop over different decomposition bases
    for (int iBase = 0; iBase < Param::B_bsk_size; iBase++)
    {
        vector<vector<NGSFFTctxt>> base_row;
        for (int iCoef = coef_counter; iCoef < coef_counter+param.bsk_partition[iBase]; iCoef++)
        {
            vector<NGSFFTctxt> coef_row;
            int sk_base_coef = sk_base.sk[iCoef][0];
            /** 
             * represent coefficient of the secret key 
             * of the base scheme using 2 bits.
             * The representation rule is as follows:
             * -1 => [0,1]
             * 0 => [0,0]
             * 1 => [1,0]
             * */
            int coef_bits[2] = {0,0};
            if (sk_base_coef == -1)
                coef_bits[1] = 1;
            else if (sk_base_coef == 1)
                coef_bits[0] = 1;
            // encrypt each bit using the NGS scheme
            for (int iBit = 0; iBit < 2; iBit++)
            {
                NGSFFTctxt bit_row;
                enc_ngs(bit_row, coef_bits[iBit], param.l_bsk[iBase], param.B_bsk[iBase], sk_boot);
                coef_row.push_back(bit_row);
            }
            base_row.push_back(coef_row);
        }
        bsk.push_back(base_row);
        coef_counter += param.bsk_partition[iBase];
    }    

    cout << "Bootstrapping key generation: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_bsk(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot)
{
    clock_t start = clock();
    
    // index of a secret key coefficient of the base scheme
    int coef_counter = 0;

    // reset the input
    bsk.clear();
    bsk = vector<vector<NGSFFTctxt>>(Param::B_bsk_size);
    for (int i = 0; i < Param::B_bsk_size; i++)
        bsk[i] = vector<NGSFFTctxt>(param.bsk_partition[i], NGSFFTctxt(param.l_bsk[i], FFTPoly(Param::N2p1)));
    
    // loop over different decomposition bases
    for (int iBase = 0; iBase < Param::B_bsk_size; iBase++)
    {
        vector<NGSFFTctxt> base_row(param.bsk_partition[iBase], NGSFFTctxt(param.l_bsk[iBase], FFTPoly(Param::N2p1)));
        for (int iCoef = coef_counter; iCoef < coef_counter+param.bsk_partition[iBase]; iCoef++)
        {
            NGSFFTctxt coef_row(param.l_bsk[iBase], FFTPoly(Param::N2p1));
            int sk_base_coef = sk_base[iCoef];
            
            // encrypt each bit using the NGS scheme
            enc_ngs(coef_row, sk_base_coef, param.l_bsk[iBase], param.B_bsk[iBase], sk_boot);
            base_row[iCoef-coef_counter] = coef_row;
        }
        bsk[iBase] = base_row;
        coef_counter += param.bsk_partition[iBase];
    }    

    cout << "Bootstrapping key generation: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_bsk2(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot)
{
    clock_t start = clock();
    
    // index of a secret key coefficient of the base scheme
    int coef_counter = 0;

    // reset the input
    bsk.clear();
    bsk = vector<vector<NGSFFTctxt>>(Param::B_bsk_size);
    for (int i = 0; i < Param::B_bsk_size; i++)
        bsk[i] = vector<NGSFFTctxt>(4 * (param.bsk_partition[i] >> 1), NGSFFTctxt(param.l_bsk[i], FFTPoly(Param::N2p1)));
    
    // loop over different decomposition bases
    int bits[4];
    for (int iBase = 0; iBase < Param::B_bsk_size; iBase++)
    {
        vector<NGSFFTctxt> base_row(4 * (param.bsk_partition[iBase] >> 1), NGSFFTctxt(param.l_bsk[iBase], FFTPoly(Param::N2p1)));
        for (int iCoef = coef_counter; iCoef < coef_counter+param.bsk_partition[iBase]; iCoef+=2)
        {
            NGSFFTctxt coef_row(param.l_bsk[iBase], FFTPoly(Param::N2p1));
            // bits to encrypt: s[coef]*s[coef+1], s[coef]*(1-s[coef+1]), (1-s[coef])*s[coef+1] 
            bits[0] = sk_base[iCoef]*sk_base[iCoef+1];
            bits[1] = sk_base[iCoef]*(1-sk_base[iCoef+1]);
            bits[2] = (1-sk_base[iCoef])*sk_base[iCoef+1];
            bits[3] = (1-sk_base[iCoef])*(1-sk_base[iCoef+1]);
            // encrypt each bit using the NGS scheme
            for (int iBit = 0; iBit < 4; iBit++)
            {
                enc_ngs(coef_row, bits[iBit], param.l_bsk[iBase], param.B_bsk[iBase], sk_boot);
                base_row[4*((iCoef-coef_counter) >> 1)+iBit] = coef_row;
            } 
        }
        bsk[iBase] = base_row;
        coef_counter += param.bsk_partition[iBase];
    }    

    cout << "Bootstrapping2 generation: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void enc_ngs(NGSFFTctxt& ct, int m, int l, int B, const SKey_boot& sk_boot)
{
    ModQPoly msg(Param::N,0L);
    msg[0] = m; // msg = m (degree-0 polynomial)
    enc_ngs(ct, msg, l, B, sk_boot);
}

void mult_fft_poly_by_int(FFTPoly& a, const int b){
    for(int i = 0; i < a.size(); i++)
        a[i] *= b;
}

void enc_ngs(NGSFFTctxt& ct, const ModQPoly& m, int l, int B, const SKey_boot& sk_boot)
{
    if(ct.size() != l)
        ct = NGSFFTctxt(l);

    FFTPoly sk_boot_inv_fft(Param::N2p1); // f^-1 in FFT form
    fftN.to_fft(sk_boot_inv_fft, sk_boot.sk_inv);
    FFTPoly g_fft(Param::N2p1);
    ModQPoly msg(m); // at each iteration i, msg will be equal to m * B^i
    FFTPoly msg_powB(Param::N2p1);
    fftN.to_fft(msg_powB, msg); // FFT of m * B^i
    FFTPoly tmp_ct(Param::N2p1);
    vector<long> tmp_ct_long(Param::N);
    vector<int> tmp_ct_int(Param::N);

    for (int i = 0; i < l; i++)
    {
        // sample random ternary vector
        ModQPoly g(Param::N,0L);
        Sampler::get_ternary_vector(g);
        // FFT transform it
        fftN.to_fft(g_fft, g);
        // compute g * sk_boot^(-1)
        tmp_ct = g_fft * sk_boot_inv_fft;
        // compute g * sk_boot^(-1) + B^i * m
        tmp_ct += msg_powB;
        // inverse FFT of the above result
        fftN.from_fft(tmp_ct_long, tmp_ct);
        // reduction modulo q_boot
        mod_q_boot(tmp_ct_int, tmp_ct_long);
        // FFT transform for further use
        fftN.to_fft(tmp_ct, tmp_ct_int);

        ct[i] = tmp_ct;

        mult_fft_poly_by_int(msg_powB, B); // msg_powB = msg * B^i
    }
}

