#include <iostream>
#include <cassert>
#include "params.h"
#include "sampler.h"
#include "keygen.h"
#include "fft.h"
#include "ntruhe.h"
#include "lwehe.h"

#include <time.h>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <limits.h>

#include <NTL/ZZX.h>

using namespace std;
using namespace NTL;

void test_params()
{
    {
        Param param(LWE);
        cout << "Ciphertext modulus of the base scheme (LWE): " << param.q_base << endl;
        cout << "Dimension of the base scheme (LWE): " << param.n << endl;
        cout << "Ciphertext modulus for bootstrapping (LWE): " << q_boot << endl;
        cout << "Polynomial modulus (LWE): " << Param::get_def_poly() << endl;
        assert(param.l_ksk == int(ceil(log(double(param.q_base))/log(double(Param::B_ksk)))));
        cout << "Decomposition length for key-switching (LWE): " << param.l_ksk << endl;
        cout << "Decomposition bases for key-switching (LWE): " << Param::B_ksk << endl;
        cout << "Dimension for bootstrapping (LWE): " << Param::N << endl;
        cout << "Decomposition bases for bootstrapping (LWE): ";
        for (const auto &v: param.B_bsk) cout << v << ' ';
        cout << endl;
        cout << "Delta (LWE): " << param.delta_base << endl;
        cout << "Half Delta (LWE): " << param.half_delta_base << endl;
    }
    {
        Param param(NTRU);
        cout << "Ciphertext modulus of the base scheme (MNTRU): " << param.q_base << endl;
        cout << "Dimension of the base scheme (NTRU): " << param.n << endl;
        cout << "Ciphertext modulus for bootstrapping (NTRU): " << q_boot << endl;
        cout << "Polynomial modulus (NTRU): " << Param::get_def_poly() << endl;
        assert(param.l_ksk == int(ceil(log(double(param.q_base))/log(double(Param::B_ksk)))));
        cout << "Decomposition length for key-switching (MNTRU): " << param.l_ksk << endl;
        cout << "Decomposition bases for key-switching (MNTRU): " << Param::B_ksk << endl;
        cout << "Dimension for bootstrapping (MNTRU): " << Param::N << endl;
        cout << "Decomposition bases for bootstrapping (MNTRU): ";
        for (const auto &v: param.B_bsk) cout << v << ' ';
        cout << endl;
        cout << "Decomposition lengths for bootstrapping (MNTRU): ";
        for (int i = 0; i < Param::B_bsk_size; i++) 
        {
            assert(param.l_bsk[i] == int(ceil(log(double(q_boot))/log(double(param.B_bsk[i])))));
            cout << param.l_bsk[i] << ' ';
        }
        cout << endl;
        cout << "Decomposition lengths for bootstrapping (MNTRU): ";
        for (int i = 0; i < Param::B_bsk_size; i++) 
        {
            assert(param.l_bsk[i] == int(ceil(log(double(q_boot))/log(double(param.B_bsk[i])))));
            cout << param.l_bsk[i] << ' ';
        }
        cout << endl;
        cout << "Delta (MNTRU): " << param.delta_base << endl;
        cout << "Half Delta (MNTRU): " << param.half_delta_base << endl;

        {
            assert(0L == mod_q_boot(0L));
            assert(1L == mod_q_boot(1L));
            assert(0L == mod_q_boot(q_boot));
            assert(half_q_boot == mod_q_boot(half_q_boot));
            assert(-half_q_boot == mod_q_boot(-half_q_boot));
            cout << "MODULO REDUCTION IS OK" << endl;
        }
    }
    
    cout << "Plaintext modulus: " << Param::t << endl;
    cout << endl;
    cout << "PARAMS ARE OK" << endl;

    {
        vector<int> res;
        decompose(res, 0, 2, 3);
        assert(res.size() == 3);
        for (auto iter=res.begin(); iter < res.end(); iter++)
            assert(0L == *iter);
    }
    {
        vector<int> res;
        decompose(res, 1, 2, 3);
        assert(res.size() == 3);
        assert(res[0] == 1);
        for (auto iter=res.begin()+1; iter < res.end(); iter++)
            assert(0L == *iter);
    }
    {
        vector<int> res;
        decompose(res, 2, 3, 3);
        assert(res.size() == 3);
        assert(res[0] == -1 && res[1] == 1 && res[2] == 0);
    }
    {
        vector<int> res;
        decompose(res, 2, 4, 3);
        assert(res.size() == 3);
        assert(res[0] == 2 && res[1] == 0 && res[2] == 0);
        decompose(res, 3, 4, 3);
        assert(res.size() == 3);
        assert(res[0] == -1 && res[1] == 1 && res[2] == 0);
    }
    {
        vector<int> res;
        try
        {
            decompose(res, 14, 3, 3);
            assert(false);
        }
        catch (overflow_error)
        {
            assert(true);
        }
    }
    {
        vector<int> res;
        try
        {
            decompose(res, -14, 3, 3);
            assert(false);
        }
        catch (overflow_error)
        {
            assert(true);
        }
    }
    {
        vector<int> res;
        decompose(res, 13, 3, 3);
        assert(res.size() == 3);
        assert(res[0] == 1 && res[1] == 1 && res[2] == 1);
        decompose(res, -13, 3, 3);
        assert(res.size() == 3);
        assert(res[0] == -1 && res[1] == -1 && res[2] == -1);
    }


    cout << "DECOMPOSITION IS OK" << endl;

    
}

void test_sampler()
{
    int N = Param::N;

    Param pLWE(LWE);
    Param pNTRU(NTRU);
    for (int run = 0; run < 1; run++)
    {
        //cout << "Run: " << run+1 << endl;
        {
            vector<int> vec(pNTRU.n, 0L);
            Sampler::get_ternary_vector(vec);
            
            assert(vec.size() == pNTRU.n);
            for (int i = 0; i < pNTRU.n; i++)
            {
                assert((vec[i]==0) || (vec[i]==-1) || (vec[i]==1) );
            }
        }

        {
            vector<int> vec(N,0L);
            Sampler::get_ternary_vector(vec);
            
            assert(vec.size() == N);
            for (int i = 0; i < N; i++)
            {
                assert((vec[i]==0) || (vec[i]==-1) || (vec[i]==1) );
            }
        }

        {
            vector<int> vec(N,0L);
            Sampler::get_binary_vector(vec);
            
            assert(vec.size() == N);
            for (int i = 0; i < N; i++)
            {
                assert((vec[i]==0) || (vec[i]==1) );
            }
        }

        {
            int n = pLWE.n;
            vector<vector<int>> mat(n, vector<int>(N,0L));
            Sampler::get_ternary_matrix(mat);
            
            assert(mat.size() == n && mat[0].size() == N);
            for (int i = 0; i < n; i++)
            {
                vector<int>& row = mat[i];
                for (int j = 0; j < N; j++)
                    assert((row[j]==0) || (row[j]==-1) || (row[j]==1) );
            }
        }

        {
            int n = pLWE.n;
            vector<int> vec(n, 0L);
            double st_dev = 4.0;
            Sampler::get_gaussian_vector(vec, st_dev);
            
            assert(vec.size() == n);
            for (int i = 0; i < n; i++)
            {
                assert(conv<double>(abs(vec[i])) < 6*st_dev);
            }
        }

        {
            int n = pNTRU.n;
            vector<vector<int>> mat(n, vector<int>(N,0L));
            double st_dev = 4.0;
            Sampler::get_gaussian_matrix(mat, st_dev);
            
            assert(mat.size() == n && mat[0].size() == N);
            for (int i = 0; i < n; i++)
            {
                vector<int>& row = mat[i];
                for (int j = 0; j < N; j++)
                    assert(conv<double>(abs(row[j])) < 6*st_dev);
            }
        }

        {
            vector<int> vec_inv(N,0L);
            vector<int> vec(N,0L);
            Sampler s(pNTRU);
            s.get_invertible_vector(vec, vec_inv, 4, 1);
            
            assert(vec.size() == N && vec_inv.size() == N);
            assert((vec[0]==1) || (vec[0]==-3) || (vec[0]==5) );
            for (int i = 1; i < N; i++)
            {
                assert((vec[i]==0) || (vec[i]==-4) || (vec[i]==4));
            }
        }

        {
            int n = pLWE.n;
            vector<vector<int>> mat_inv(n, vector<int>(n,0L));
            vector<vector<int>> mat(n, vector<int>(n,0L));
            Sampler s(pLWE);
            s.get_invertible_matrix(mat, mat_inv, 5, 1);
            
            assert(mat.size() == n && mat[0].size() == n 
                && mat_inv.size() == n && mat_inv[0].size() == n);
            for (int i = 0; i < n; i++)
                assert((mat[i][i]==1) || (mat[i][i]==-4) || (mat[i][i]==6) );
            
            for (int i = 0; i < n; i++)
                for (int j = 0; (j < n) && (j != i); j++)
                {
                    assert((mat[i][j]==0) || (mat[i][j]==-5) || (mat[i][j]==5) );
                }
        }
    }
    cout << "SAMPLER IS OK" << endl;
}

void test_ntru_key_gen()
{
    Param param(NTRU);
    int n = param.n;
    int Nl = param.Nl;
    int half_q_base = param.half_q_base;
    int q_base = param.q_base;
    int l_ksk = param.l_ksk;
    int N = Param::N;
    int t = Param::t;
    int B_ksk = Param::B_ksk;
    int B_bsk_size = Param::B_bsk_size;
    int N2p1 = Param::N2p1;

    SKey_base_NTRU sk_base;
    KeyGen k(param);
    k.get_sk_base(sk_base);
    cout << "Secret key of the base scheme is generated" << endl;
    assert(sk_base.sk.size() == n && sk_base.sk[0].size() == n 
        && sk_base.sk_inv.size() == n && sk_base.sk_inv[0].size() == n);  
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            assert((sk_base.sk[i][j]==0) || (sk_base.sk[i][j]==-1) || (sk_base.sk[i][j]==1) );
        }

    SKey_boot sk_boot;
    k.get_sk_boot(sk_boot);
    cout << "Secret key of the bootstrapping scheme is generated" << endl;
    assert(sk_boot.sk.size() == N && sk_boot.sk_inv.size() == N);
    assert((sk_boot.sk[0]==1) || (sk_boot.sk[0]==(-t+1)) || (sk_boot.sk[0]==(t+1))); 
    for (int i = 1; i < N; i++)
    {
        assert((sk_boot.sk[i]==0) || (sk_boot.sk[i]==-t) || (sk_boot.sk[i]==t) );
    }

    KSKey_NTRU ksk;
    k.get_ksk(ksk, sk_base, sk_boot);
    cout << "Key-switching key is generated" << endl;
    assert(ksk.size() == Nl && ksk[0].size() == n);
    for (int i = 0; i < Nl; i++)
        for (int j = 0; j < n; j++)
        {
            //cout << ksk[i][j] << endl;
            assert(ksk[i][j] <= half_q_base && ksk[i][j] >= -half_q_base);
        }

    vector<int> q4_decomp;
    decompose(q4_decomp, q_base/4, B_ksk, l_ksk);
    vector<int> ks_res(n,0L);
    for (int i = 0; i < l_ksk; i++)
    {
        int tmp_int = q4_decomp[i];
        vector<int>& ksk_row = ksk[i];
        for (int j = 0; j < n; j++)
        {
            ks_res[j] += ksk_row[j] * tmp_int;
        }
    }
    param.mod_q_base(ks_res);
    int ks_int = 0;
    for (int i = 0; i < n; i++)
    {
        ks_int += ks_res[i] * sk_base.sk[i][0];
    }
    ks_int = param.mod_q_base(ks_int);
    ks_int = int(round(double(ks_int*4)/double(q_base)));
    assert(ks_int == 1L);

    // bootstrapping key test
    BSKey_NTRU bsk;
    k.get_bsk(bsk, sk_base, sk_boot);
    cout << "Bootstrapping key is generated" << endl;

    // check dimensions
    assert(bsk.size() == B_bsk_size);
    for (int i = 0;  i < bsk.size(); i++)
    {
        assert(bsk[i].size() == param.bsk_partition[i]);
        for (int j = 0; j < bsk[i].size(); j++)
        {
            assert(bsk[i][j].size() == 2);
            assert(bsk[i][j][0].size() == param.l_bsk[i]);
            assert(bsk[i][j][1].size() == param.l_bsk[i]);
        }
    }
    // convert sk_boot to FFT
    vector<complex<double>> sk_boot_fft(N2p1);
    fftN.to_fft(sk_boot_fft, sk_boot.sk);

    int coef_counter = 0;
    for (int iBase = 0; iBase < B_bsk_size; iBase++)
    {
        decompose(q4_decomp, q_boot/4, param.B_bsk[iBase], param.l_bsk[iBase]);
        for (size_t iCoef = 0; iCoef < bsk[iBase].size(); iCoef++)
        {
            int sk_coef = 0;
            int sk_base_coef_bits[2];
            for (int iBit = 0; iBit < 2; iBit++)
            {
                vector<complex<double>> tmp_fft(N2p1, complex<double>(0.0,0.0));
                for (int iPart = 0; iPart < param.l_bsk[iBase]; iPart++)
                {
                    tmp_fft = tmp_fft + bsk[iBase][iCoef][iBit][iPart] * q4_decomp[iPart];
                }
                tmp_fft = tmp_fft * sk_boot_fft;
                vector<int> tmp_int;
                vector<long> tmp_long;
                fftN.from_fft(tmp_long, tmp_fft);
                mod_q_boot(tmp_int, tmp_long);
                sk_base_coef_bits[iBit] = int(round(double(tmp_int[0]*4)/double(q_boot)));
            }
            if (sk_base_coef_bits[1] == 1)
                sk_coef = -1;
            else if (sk_base_coef_bits[0] == 1)
                sk_coef = 1;

            assert(sk_coef == sk_base.sk[coef_counter + iCoef][0]);
        }
        coef_counter += param.bsk_partition[iBase];
    }

    cout << "KEYGEN IS OK" << endl;
}

void test_lwe_key_gen()
{
    Param param(LWE);
    int n = param.n;
    int N = Param::N;
    int t = Param::t;

    SKey_base_LWE sk_base;
    KeyGen k(param);
    k.get_sk_base(sk_base);
    cout << "Secret key of the base scheme is generated" << endl;
    assert(sk_base.size() == n);  
    for (int j = 0; j < n; j++)
    {
        assert((sk_base[j]==0) || (sk_base[j]==1));
    }

    SKey_boot sk_boot;
    k.get_sk_boot(sk_boot);
    cout << "Secret key of the bootstrapping scheme is generated" << endl;
    assert(sk_boot.sk.size() == N && sk_boot.sk_inv.size() == N);
    assert((sk_boot.sk[0]==1) || (sk_boot.sk[0]==(-t+1)) || (sk_boot.sk[0]==(t+1))); 
    for (int i = 1; i < N; i++)
    {
        assert((sk_boot.sk[i]==0) || (sk_boot.sk[i]==-t) || (sk_boot.sk[i]==t) );
    }

    cout << "KEYGEN IS OK" << endl;
}

void test_fft()
{
    int N = Param::N;
    int N2p1 = Param::N2p1;

    FFT_engine fft_engine(N);
    {
        vector<int> in(N,0L);
        vector<complex<double>> out(N2p1);
        clock_t start = clock();
        fft_engine.to_fft(out, in);
        cout << "Forward FFT (zero): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        for (size_t i = 0; i < N/2; i++)
        {
            if (int(round(real(out[i])))!=0 || int(round(imag(out[i])))!=0)
            {
                cout << i << " " << out[i] << endl;
                assert(false);
            }
        }
    }

    {
        vector<long> out;
        vector<complex<double>> in(N2p1, complex<double>(0.0,0.0));
        clock_t start = clock();
        fft_engine.from_fft(out, in);
        cout << "Backward FFT (zero): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        for (size_t i = 0; i < N/2; i++)
        {
            assert(out[i] == 0L);
        }
    }

    {
        vector<int> in(N,0L);
        in[0] = 1L;
        vector<complex<double>> out(N2p1);
        clock_t start = clock();
        fft_engine.to_fft(out, in);
        cout << "Forward FFT (1,0,...0): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        for (size_t i = 0; i < N/2; i++)
        {
            assert(int(round(real(out[i])))==1 && int(round(imag(out[i])))==0);
        }
    }

    {
        vector<long> out;
        vector<complex<double>> in(N2p1, complex<double>(1.0,0.0));
        clock_t start = clock();
        fft_engine.from_fft(out, in);
        cout << "Backward FFT (1,1,...1): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        assert(out[0] == 1L);
        for (size_t i = 1; i < N; i++)
        {
            assert(out[i] == 0L);
        }
    }

    {
        uniform_int_distribution<int> sampler(INT_MIN, INT_MAX);
        int coef = sampler(rand_engine);
        vector<long> out;
        vector<complex<double>> in(N2p1, complex<double>(double(coef),0.0));
        clock_t start = clock();
        fft_engine.from_fft(out, in);
        cout << "Backward FFT (a,a,...a): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        assert(out[0] == coef);
        for (size_t i = 1; i < N; i++)
        {
            assert(out[i] == 0L);
        }
    }

    {
        uniform_int_distribution<int> sampler(INT_MIN, INT_MAX);
        vector<int> in;
        for (int i = 0; i < N; i++)
            in.push_back(sampler(rand_engine));
        vector<complex<double>> interm(N2p1);
        vector<long> out;
        clock_t start = clock();
        fft_engine.to_fft(interm, in);
        cout << "Forward FFT (random): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        start = clock();
        fft_engine.from_fft(out, interm);
        cout << "Backward FFT (random): " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        for (size_t i = 0; i < N; i++)
        {
            //cout << "i: " << i << "in[i]: " << in[i] << " out[i]: " << out[i] << endl; 
            assert(in[i] == out[i]);
        }
    }

    {
        uniform_int_distribution<int> sampler(-100, 100);
        vector<int> in1, in2, res;
        for (int i = 0; i < N; i++)
        {
            in1.push_back(sampler(rand_engine));
            in2.push_back(sampler(rand_engine));
            res.push_back(in1[i]+in2[i]);
        }
        vector<complex<double>> interm1(N2p1);
        vector<complex<double>> interm2(N2p1);
        vector<complex<double>> intermres(N2p1);
        vector<long> out;
        
        fft_engine.to_fft(interm1, in1);
        fft_engine.to_fft(interm2, in2);

        clock_t start = clock();
        intermres = interm1 + interm2;
        cout << "FFT addition: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        
        fft_engine.from_fft(out, intermres);
        
        for (size_t i = 0; i < N; i++)
        {
            //cout << "i: " << i << " in1[i]: " << in1[i] << " in2[i]: " << in2[i]  << " res[i]: " << res[i] << " out[i]: " << out[i] << endl; 
            assert(res[i] == out[i]);
        }
    }

    {
        uniform_int_distribution<int> sampler(-100, 100);
        vector<int> in1, in2, res;
        for (int i = 0; i < N; i++)
        {
            in1.push_back(sampler(rand_engine));
            in2.push_back(sampler(rand_engine));
        }
        ZZX poly1, poly2, poly_res;
        for (int i = 0; i < N; i++)
        {
            SetCoeff(poly1, i, in1[i]);
            SetCoeff(poly2, i, in2[i]);
        }
        ZZX poly_mod;
        SetCoeff(poly_mod, 0, 1);
        SetCoeff(poly_mod, N, 1);

        clock_t start = clock();
        MulMod(poly_res, poly1, poly2, poly_mod);
        cout << "NTL multiplication: " << float(clock()-start)/CLOCKS_PER_SEC << endl;

        for (int i = 0; i < N; i++)
        {
            res.push_back(conv<long>(poly_res[i]));
        }

        vector<complex<double>> interm1(N2p1);
        vector<complex<double>> interm2(N2p1);
        vector<complex<double>> intermres(N2p1);
        vector<long> out;
        
        fft_engine.to_fft(interm1, in1);
        fft_engine.to_fft(interm2, in2);

        start = clock();
        intermres = interm1 * interm2;
        cout << "FFT multiplication: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
        
        fft_engine.from_fft(out, intermres);
        
        for (size_t i = 0; i < N; i++)
        {
            //cout << "i: " << i << " in1[i]: " << in1[i] << " in2[i]: " << in2[i]  << " res[i]: " << res[i] << " out[i]: " << out[i] << endl; 
            assert(res[i] == out[i]);
        }
    }

    cout << "FFT is OK" << endl;
}

void test_ntruhe_encrypt()
{
    SchemeNTRU s;

    {
        int input = 0;
        Ctxt_NTRU ct;
        s.encrypt(ct, input);
        int output = s.decrypt(ct);
        assert(output == input);
    }
    {
        int input = 1;
        Ctxt_NTRU ct;
        s.encrypt(ct, input);
        int output = s.decrypt(ct);
        assert(output == input);
    }
    cout << "NTRU ENCRYPTION IS OK" << endl;
}

void test_lwehe_encrypt()
{
    SchemeLWE s;

    {
        int input = 0;
        Ctxt_LWE ct;
        s.encrypt(ct, input);
        int output = s.decrypt(ct);
        assert(output == input);
    }
    {
        int input = 1;
        Ctxt_LWE ct;
        s.encrypt(ct, input);
        int output = s.decrypt(ct);
        assert(output == input);
    }
    cout << "LWE ENCRYPTION IS OK" << endl;
}

/*
void test_mod_switch()
{
    SchemeNTRU s;
    {
        int input = 0;
        Ctxt_NTRU ct;
        s.encrypt(ct, input);
        s.modulo_switch_to_base(ct.data);
        int output = 0;
        for (int i = 0; i < ntru_he::n; i++)
        {
            output += ct.data[i] * sk_base.sk[i][0];
        }
        output = output%Param::N2;
        if (output > Param::N)
            output -= Param::N2;
        else if (output <= -Param::N)
            output += Param::N2;
        output = int(round(double(output*t)/double(Param::N2)));
        assert(output == input);
    }
    {
        int input = 1;
        ntru_he::Ctxt ct;
        ntru_he::encrypt(ct, input, sk_base);
        ntru_he::modulo_switch(ct, ntru_he::q_base, Param::N2);
        int output = 0;
        for (int i = 0; i < ntru_he::n; i++)
        {
            output += ct[i] * sk_base.sk[i][0];
        }
        output = output%Param::N2;
        if (output > Param::N)
            output -= Param::N2;
        else if (output <= -Param::N)
            output += Param::N2;
        output = int(round(double(output*t)/double(Param::N2)));
        assert(output == input);
    }
    {
        int input = 0;
        ntru_he::Ctxt ct;
        ntru_he::encrypt(ct, input, sk_base);
        ntru_he::modulo_switch_to_boot(ct);
        int output = 0;
        for (int i = 0; i < ntru_he::n; i++)
        {
            output += ct[i] * sk_base.sk[i][0];
        }
        output = output%Param::N2;
        if (output > Param::N)
            output -= Param::N2;
        else if (output <= -Param::N)
            output += Param::N2;
        output = int(round(double(output*t)/double(Param::N2)));
        assert(output == input);
    }
    {
        int input = 1;
        ntru_he::Ctxt ct;
        ntru_he::encrypt(ct, input, sk_base);
        ntru_he::modulo_switch_to_boot(ct);
        int output = 0;
        for (int i = 0; i < ntru_he::n; i++)
        {
            output += ct[i] * sk_base.sk[i][0];
        }
        output = output%Param::N2;
        if (output > Param::N)
            output -= Param::N2;
        else if (output <= -Param::N)
            output += Param::N2;
        output = int(round(double(output*t)/double(Param::N2)));
        assert(output == input);
    }
    cout << "MODULO SWITCHING IS OK" << endl;
}
*/

void test_bootstrap()
{
    SchemeNTRU s;
    {
        int input = 2;
        Ctxt_NTRU ct;
        s.encrypt(ct, input);

        s.bootstrap(ct);

        int output = s.decrypt(ct);
        cout << "Bootstrapping output: " << output << endl;
        assert(output == 1L);
    }
    {
        int input = 0;
        Ctxt_NTRU ct;
        s.encrypt(ct, input);

        s.bootstrap(ct);

        int output = s.decrypt(ct);
        cout << "Bootstrapping output: " << output << endl;
        assert(output == 0L);
    }

    cout << "BOOTSTRAPPING IS OK" << endl;
}

/*
void test_nand_aux()
{
    ntru_he::SKey_base sk_base;
    ntru_he::get_sk_base(sk_base);

    ntru_he::Ctxt ct;
    ntru_he::get_nand_aux(ct, sk_base);
    int output = 0;
    for (int i = 0; i < ntru_he::n; i++)
    {
        output += ct[i] * sk_base.sk[i][0];
    }
    output = ntru_he::mod_q_base(output);
    assert(
        output == (ntru_he::nand_const-ntru_he::q_base) 
        || output == (ntru_he::nand_const-ntru_he::q_base+1) 
        || output == (ntru_he::nand_const-ntru_he::q_base-1)
        );
    cout << "NAND ENCRYPTION IS OK" << endl;  
}*/

enum GateType {NAND, AND, OR};

void test_ntruhe_gate_helper(int in1, int in2, const SchemeNTRU& s, GateType g)
{
    float avg_time = 0.0;
    for (int i = 0; i < 100; i++)
    {
        Ctxt_NTRU ct_res, ct1, ct2, ct_nand;
        s.encrypt(ct1, in1);
        s.encrypt(ct2, in2);
        
        if (g == NAND)
        {
            auto start = clock();
            s.nand_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "NAND output: " << output << endl;
            assert(output == !(in1 & in2));
        }
        else if (g == AND) {
            auto start = clock();
            s.and_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "AND output: " << output << endl;
            assert(output == (in1 & in2));
        }
        else if (g == OR) {
            auto start = clock();
            s.or_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "OR output: " << output << endl;
            assert(output == (in1 | in2));
        }
    }
    cout << "Avg. time" << avg_time/100.0 << endl;
}

void test_ntru_gate(GateType g)
{
    SchemeNTRU s;

    test_ntruhe_gate_helper(0, 0, s, g);
    test_ntruhe_gate_helper(0, 1, s, g);
    test_ntruhe_gate_helper(1, 0, s, g);
    test_ntruhe_gate_helper(1, 1, s, g);
}

void test_ntruhe_nand()
{
    GateType g = NAND;

    test_ntru_gate(g);

    cout << "NAND IS OK" << endl;
}

void test_ntruhe_and()
{
    GateType g = AND;

    test_ntru_gate(g);

    cout << "AND IS OK" << endl;
}

void test_ntruhe_or()
{
    GateType g = OR;

    test_ntru_gate(g);

    cout << "OR IS OK" << endl;
}

void test_lwehe_gate_helper(int in1, int in2, SchemeLWE& s, GateType g)
{
    float avg_time = 0.0;
    for (int i = 0; i < 100; i++)
    {
        Ctxt_LWE ct_res, ct1, ct2, ct_nand;
        s.encrypt(ct1, in1);
        s.encrypt(ct2, in2);
        
        if (g == NAND)
        {
            auto start = clock();
            s.nand_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "NAND output: " << output << endl;
            assert(output == !(in1 & in2));
        }
        else if (g == AND) {
            auto start = clock();
            s.and_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "AND output: " << output << endl;
            assert(output == (in1 & in2));
        }
        else if (g == OR) {
            auto start = clock();
            s.or_gate(ct_res, ct1, ct2);
            avg_time += float(clock()-start)/CLOCKS_PER_SEC;

            int output = s.decrypt(ct_res);

            //cout << "OR output: " << output << endl;
            assert(output == (in1 | in2));
        }
    }
    cout << "Avg. time" << avg_time/100.0 << endl;
}

void test_lwe_gate(GateType g)
{
    SchemeLWE s;

    test_lwehe_gate_helper(0, 0, s, g);
    test_lwehe_gate_helper(0, 1, s, g);
    test_lwehe_gate_helper(1, 0, s, g);
    test_lwehe_gate_helper(1, 1, s, g);
}

void test_lwehe_nand()
{
    GateType g = NAND;

    test_lwe_gate(g);

    cout << "NAND IS OK" << endl;
}

void test_lwehe_and()
{
    GateType g = AND;

    test_lwe_gate(g);

    cout << "AND IS OK" << endl;
}

void test_lwehe_or()
{
    GateType g = OR;

    test_lwe_gate(g);

    cout << "OR IS OK" << endl;
}

int main()
{
    //test_params();
    //test_sampler();
    //test_ntru_key_gen();
    //test_lwe_key_gen();
    //test_fft();
    //test_ntruhe_encrypt();
    //test_lwehe_encrypt();
    //test_mod_switch();
    //test_bootstrap();
    //test_nand_aux();
    //test_ntruhe_nand();
    test_lwehe_nand();
    //test_ntruhe_and();
    test_lwehe_and();
    //test_ntruhe_or();
    test_lwehe_or();
    return 0;
}