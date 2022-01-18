#ifndef KEYGEN
#define KEYGEN

#include <NTL/mat_ZZ.h>
#include <vector>
#include <fftw3.h>
#include <complex>
#include "params.h"
#include "sampler.h"

// secret key of the bootstrapping scheme
typedef struct {
    ModQPoly sk;
    ModQPoly sk_inv;
} SKey_boot;

// secret key of the NTRU base scheme
typedef struct {
    ModQMatrix sk;
    ModQMatrix sk_inv;
} SKey_base_NTRU;

// secret key of the LWE base scheme
typedef std::vector<int> SKey_base_LWE;

/** 
 * Bootstrapping key.
 * It consists of several sets of keys corresponding to different
 * decomposition bases of the bootstrapping key B_bsk.
 * The i-th set contains vectors with l_bsk[i] complex vectors.
 * These complex vectors are an encryption of some bit of the secret key
 * of the base scheme in the NGS form.
 */
typedef std::vector<std::vector<std::vector<NGSFFTctxt>>> BSKey_NTRU;

/** 
 * Bootstrapping key.
 * It consists of several sets of keys corresponding to different
 * decomposition bases of the bootstrapping key B_bsk.
 * The i-th set contains vectors with l_bsk[i] complex vectors.
 * These complex vectors are an encryption of some bit of the secret key
 * of the base scheme in the NGS form.
 */
typedef std::vector<std::vector<NGSFFTctxt>> BSKey_LWE;

// key-switching key from NTRU to NTRU
typedef ModQMatrix KSKey_NTRU;

// key-switching key from NTRU to LWE
typedef struct{
    ModQMatrix A;
    std::vector<int> b;
}  KSKey_LWE;


class KeyGen
{
    Param param;
    Sampler sampler;

    public:

        KeyGen(Param _param): param(_param), sampler(_param)
        {}

        /**
         * Generate a secret key of the bootstrapping scheme. 
         * @param[out] sk_boot secret key of the bootstrapping scheme.
         */
        void get_sk_boot(SKey_boot& sk_boot);

        /**
         * Generate a secret key of the base scheme. 
         * @param[out] sk_base secret key of the base scheme.
         */
        void get_sk_base(SKey_base_NTRU& sk_base);

        /**
         * Generate a secret key of the base scheme. 
         * @param[out] sk_base secret key of the base scheme.
         */
        void get_sk_base(SKey_base_LWE& sk_base);

        /**
         * Generate a key-switching key from the bootstrapping scheme to the base scheme. 
         * @param[out] ksk key-switching key.
         * @param[in] sk_base secret key of the base scheme.
         * @param[in] sk_boot secret key of the bootstrapping scheme.
         */
        void get_ksk(KSKey_NTRU& ksk, const SKey_base_NTRU& sk_base, const SKey_boot& sk_boot);

        /**
         * Generate a key-switching key from the bootstrapping scheme to the base scheme. 
         * @param[out] ksk key-switching key.
         * @param[in] sk_base secret key of the base scheme.
         * @param[in] sk_boot secret key of the bootstrapping scheme.
         */
        void get_ksk(KSKey_LWE& ksk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot);

        /**
         * Generate a bootstrapping key
         * @param[out] bsk bootstrapping key.
         * @param[in] sk_base secret key of the base scheme.
         * @param[in] sk_boot secret key of the bootstrapping scheme.
         */
        void get_bsk(BSKey_NTRU& bsk, const SKey_base_NTRU& sk_base, const SKey_boot& sk_boot);

        /**
         * Generate a bootstrapping key
         * @param[out] bsk bootstrapping key.
         * @param[in] sk_base secret key of the base scheme.
         * @param[in] sk_boot secret key of the bootstrapping scheme.
         */
        void get_bsk(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot);

        /**
         * Generate a bootstrapping key (EXPERIMENTAL)
         * @param[out] bsk bootstrapping key.
         * @param[in] sk_base secret key of the base scheme.
         * @param[in] sk_boot secret key of the bootstrapping scheme.
         */
        void get_bsk2(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot);
};

#endif