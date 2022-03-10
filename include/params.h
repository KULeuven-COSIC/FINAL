#ifndef PARAMS
#define PARAMS

#include <NTL/ZZ_pX.h>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <complex>
#include <cassert>

using namespace NTL;

enum SchemeType {NTRU, LWE};

// representation of a polynomial modulo some integer
typedef std::vector<int> ModQPoly;
// matrix modulo some integer
typedef std::vector<std::vector<int>> ModQMatrix;
// representation of an FFT transformation of some poly
typedef std::vector<std::complex<double>> FFTPoly;
// NGS ciphertest in NTT form
typedef std::vector<FFTPoly> NGSFFTctxt;

/**
 * Reduction modulo q in the symmetric interval (-q/2, q/2]
 * Correct only for odd q.
 * @param[in] input integer to reduce.
 * @param[in] q modulus
 * @param[in] half_q (q-1)/2
 * @returns reduced integer in the symmetric interval (-q/2, q/2]
 */
inline long lazy_mod_q(const long input, const long q, const long half_q)
{
    int coef = input%q;
    if (coef > half_q)
        return coef - q;
    if (coef < -half_q)
        return coef + q;
    return coef;
}

inline int lazy_mod_q(const int input, const int q, const int half_q)
{
    int coef = input%q;
    if (coef > half_q)
        return coef - q;
    if (coef < -half_q)
        return coef + q;
    return coef;
}

inline void lazy_mod_q(std::vector<int>& output, const std::vector<long>& input, const long q, const long half_q)
{
    assert(output.size() == input.size());

    std::vector<int>::iterator oit = output.begin();
    for (auto iit = input.begin(); iit < input.end(); iit++, oit++)
        *oit = static_cast<int>(lazy_mod_q(*iit, q, half_q));
}

/** ciphertext modulus of the ring-based scheme used 
 * for bootstrapping keys and test vectors/lookup-table encodings
 */
const int q_boot = 912829; // ~2^19.8, prime
const long q_boot_long = long(q_boot);

//half of the above modulus
const int half_q_boot = q_boot/2;
const long half_q_boot_long = long(half_q_boot);

/**
 * Reduction modulo q_boot in the symmetric interval (-q_boot/2, q_boot/2]
 * Correct only for odd q_boot.
 * @param[in] input integer to reduce.
 * @returns reduced integer in the symmetric interval (-q_boot/2, q_boot/2]
 */
inline int mod_q_boot(const int input)
{
    return lazy_mod_q(input, q_boot, half_q_boot);
}

inline int mod_q_boot(const long input)
{
    return lazy_mod_q(input, long(q_boot), long(half_q_boot));
}

/**
 * Reduction modulo q_boot in the symmetric interval (-q_boot/2, q_boot/2]
 * Correct only for odd q_boot.
 * @param[in,out] input vector to reduce.
 * @param[in] q integer modulus.
 */
inline void mod_q_boot(std::vector<int>& input)
{
    for (auto it = input.begin(); it < input.end(); it++)
        *it = mod_q_boot(*it);
}

/**
 * Reduction modulo q_boot in the symmetric interval (-q_boot/2, q_boot/2]
 * Correct only for odd q_boot.
 * @param[in,out] input vector to reduce.
 * @param[in] q integer modulus.
 */
inline void mod_q_boot(std::vector<int>& output, std::vector<long>& input)
{
    lazy_mod_q(output, input, q_boot_long, half_q_boot_long);
}

class Param
{
    public:
        // decomposition base of key-switching keys (DO NOT CHANGE!)
        const static int B_ksk = 3;

        // dimension of bootstrapping keys and test vectors/lookup-table encodings
        const static int N = 1024;
        // N/2 + 1, needed for FFT
        const static int N2p1 = N/2+1;
        // order of the cyclotomic ring, which is used for the bootstrapping scheme
        const static int N2 = 2*N;

        static ZZ_pX get_def_poly()
        {
            ZZ_pX poly;
            // polynomial modulus of the ring-based scheme
            // element of Z_(q_boot)
            ZZ_p coef;
            coef.init(ZZ(q_boot));
            // polynomial modulus X^N+1
            coef = 1;
            SetCoeff(poly, 0, coef);
            SetCoeff(poly, N, coef);

            return poly;
        }

        ZZ_pX xToNplus1;

        const static int B_bsk_size = 2;

        // plaintext modulus
        const static int t = 4;

        // Delta scalar used in bootstrapping
        const static int half_delta_boot = q_boot/(2*t);

        // standard deviation for discrete Gaussian distribution
        constexpr static double e_st_dev = 4.39;

        // Type of the base scheme
        SchemeType scheme_type;

        // ciphertext modulus of the base scheme used for encryption
        int q_base;
        // half of the above modulus
        int half_q_base;
        // dimension of the ciphertext space
        int n;

        // decomposition base of key-switching keys 
        int l_ksk;

        // number of rows of key-switching key matrices
        int Nl;

        // decomposition bases of bootstrapping keys
        int B_bsk[B_bsk_size]; 
        // binary logarithms of decomposition bases
        int shift_bsk[B_bsk_size];
        // partition of bootstrapping keys per decomposition base
        int bsk_partition[B_bsk_size];

        // decomposition lengths of bootstrapping keys
        int l_bsk[2];

        // Delta scalars used in encryption
        int half_delta_base;
        int delta_base;

        // NAND constant
        int nand_const;
        // AND constant
        int and_const;
        // OR constant
        int or_const;

        void init()
        {
            if (scheme_type == SchemeType::NTRU)
            {
                q_base = 131071; // ~2^17, prime
                n = 800;

                B_bsk[0] = 8;
                B_bsk[1] = 16;
                shift_bsk[0] = 3;
                shift_bsk[1] = 4;
                bsk_partition[0] = 750;
                bsk_partition[1] = 50;
            }
            else if (scheme_type == SchemeType::LWE)
            {
                q_base = 92683; // ~2^16.5, prime
                n = 610;

                B_bsk[0] = 8;
                B_bsk[1] = 16;
                shift_bsk[0] = 3;
                shift_bsk[1] = 4;
                bsk_partition[0] = 140;
                bsk_partition[1] = 470;
            }
            
            half_q_base = q_base/2;
            l_ksk = int(ceil(log(double(q_base))/log(double(B_ksk))));
            Nl = N * l_ksk;

            for (size_t i = 0; i < 2; i++)
                l_bsk[i] = int(ceil(log(double(q_boot))/log(double(B_bsk[i]))));

            half_delta_base = q_base/(2*t);
            delta_base = 2*half_delta_base;

            nand_const = 5*half_delta_base;
            and_const = half_delta_base;
            or_const = 7*half_delta_base;

            xToNplus1 = get_def_poly();
        }
        Param(){};
        Param(SchemeType _scheme_type): scheme_type(_scheme_type)
        {
            init();
        }
        Param(const Param &param): scheme_type(param.scheme_type)
        {
            init();
        }
        Param& operator=(const Param& param)
        {
            scheme_type = param.scheme_type;
            init();
            return *this;
        }

        bool operator == (const Param& param) const
        {
            return this->scheme_type == param.scheme_type;
        }

        /**
         * Reduction modulo q_base in the symmetric interval (-q_base/2, q_base/2]
         * Correct only for odd q_base.
         * @param[in] input integer to reduce.
         * @returns reduced integer in the symmetric interval (-q_base/2, q_base/2]
         */
        inline int mod_q_base(const int input) const
        {
            return lazy_mod_q(input, q_base, half_q_base);
        }
        inline int mod_q_base(const long input) const
        {
            return lazy_mod_q(input, long(q_base), long(half_q_base));
        }

        /**
         * Reduction modulo q_base in the symmetric interval (-q_base/2, q_base/2]
         * Correct only for odd q_base.
         * @param[in,out] input vector to reduce.
         */
        inline void mod_q_base(std::vector<int>& input) const
        {
            for (auto it = input.begin(); it < input.end(); it++)
                *it = mod_q_base(*it);
        }
        inline void mod_q_base(std::vector<int>& output, std::vector<long>& input) const
        {
            output.resize(input.size());
            for (size_t i = 0; i < input.size(); i++)
                output[i] = mod_q_base(input[i]);
        }
};

const Param parLWE(LWE);
const Param parNTRU(NTRU);

/**
 * Switches a given polynomial to a given modulus.
 * @param[in,out] poly polynomial
 * @param[in] old_q old modulus
 * @param[in] new_q old modulus
 */ 
inline void modulo_switch(ModQPoly& poly, int old_q, int new_q)
{
    double ratio = double(new_q)/double(old_q);
    for (auto it = poly.begin(); it < poly.end(); it++)
        *it = int(round(double(*it)*ratio));
}

/**
 * Balanced decomposition of an integer in base b.
 * The length of decomposition must be l. 
 * @param[out] res decomposition vector
 * @param[in] input integer to decompose
 * @param[in] b decomposition base
 * @param[in] l decomposition length
 */
inline void decompose(std::vector<int>& res, const int input, const int b, const int l)
{
    res.clear();

    int input_sign = (input < 0) ? -1: 1;
    int input_rem = abs(input);
    for (int i=0; i<l; i++)
    {
        int digit = input_rem % b;
        int digit2 = 2*digit;
        if (digit2 > b)
        {
            res.push_back(input_sign * (digit - b));
            input_rem = (input_rem - digit)/b + 1;
        }
        else
        {
            res.push_back(input_sign * digit);
            input_rem = (input_rem - digit)/b;
        }
    }
    if (input_rem != 0)
        throw std::overflow_error("Input is too big for given length\n");
}

#endif
