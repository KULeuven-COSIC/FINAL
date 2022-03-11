#ifndef LWEHE
#define LWEHE

#include "params.h"
#include "keygen.h"

class Ctxt_LWE
{
    public:
    std::vector<int> a;
    int b;

    Ctxt_LWE()
    {
        a.clear();
        a.resize(parLWE.n);
    }
    Ctxt_LWE(const Ctxt_LWE& ct);
    Ctxt_LWE& operator=(const Ctxt_LWE& ct);

    Ctxt_LWE operator +(const Ctxt_LWE& ct) const;
    Ctxt_LWE operator -(const Ctxt_LWE& ct) const;
};

Ctxt_LWE operator -(const int c, const Ctxt_LWE& ct);

/**
 * Switches a given ciphertext to a given modulus.
 * @param[in,out] ct ciphertext
 * @param[in] old_q old modulus
 * @param[in] new_q old modulus
 */ 
inline void modulo_switch_lwe(Ctxt_LWE& ct, int old_q, int new_q)
{
    std::vector<int>& a = ct.a;
    for (size_t i = 0; i < a.size(); i++)
        a[i] = int((a[i]*new_q)/old_q);

    ct.b = int((ct.b*new_q)/old_q);
}

/**
 * Switches a given polynomial from q_base to modulus 2*N.
 * @param[in,out] poly polynomial
 */ 
inline void modulo_switch_to_boot(Ctxt_LWE& poly)
{
    modulo_switch_lwe(poly, parLWE.q_base, Param::N2);
}

/**
 * Switches a given polynomial from q_boot to q_base.
 * @param[in,out] poly polynomial
 */ 
inline void modulo_switch_to_base_lwe(ModQPoly& poly)
{
    modulo_switch(poly, q_boot, parLWE.q_base);
}

/**
 * Computes the external product of a given polynomial ciphertext
 * with an NGS ciphertext in the FFT form
 * @param[in,out] poly polynomial ciphertext
 * @param[in] poly_vector NGS ciphertext
 * @param[in] b decomposition base, power of 2
 * @param[in] shift bit shift to divide by b
 * @param[in] l decomposition length
 */ 
void external_product(std::vector<long>& res, const std::vector<int>& poly, const std::vector<FFTPoly>& poly_vector, int b, int shift, int l);

class SchemeLWE
{
    SKey_base_LWE sk_base;
    SKey_boot sk_boot;
    KSKey_LWE ksk;
    BSKey_LWE bk;

    public:

        SchemeLWE()
        {
            KeyGen keygen(parLWE);

            keygen.get_sk_base(sk_base);
            keygen.get_sk_boot(sk_boot);
            keygen.get_ksk(ksk,sk_base,sk_boot);
            keygen.get_bsk(bk,sk_base,sk_boot);
        }
        /**
         * Encrypts a bit using LWE.
         * @param[out] ct ciphertext encrypting the input bit 
         * @param[in] m bit to encrypt
         */ 
        void encrypt(Ctxt_LWE& ct, int m) const;

        /**
         * Decrypts a bit using LWE.
         * @param[out] ct ciphertext encrypting a bit
         * @return b bit
         */ 
        int decrypt(const Ctxt_LWE& ct) const;
        
        /**
         * Performs key switching of a given ciphertext from a polynomial NTRU
         * to LWE
         * @param[out] ct LWE ciphertext (vector of dimension n)
         * @param[in] poly polynomial ciphertext (vector of dimension N)
         */
        void key_switch(Ctxt_LWE& ct, const ModQPoly& poly) const;

        /**
         * Bootstrapps a given ciphertext
         * @param[in,out] ct ciphertext to bootstrap
         */
        void bootstrap(Ctxt_LWE& ct) const;

        /**
         * Bootstrapps a given ciphertext
         * @param[in,out] ct ciphertext to bootstrap
         */
        void bootstrap2(Ctxt_LWE& ct) const;

        /**
         * Computes the NAND gate of two given ciphertexts ct1 and ct2
         * @param[out] ct_res encryptions of the outuput of the NAND gate
         * @param[in] ct_1 encryption of the first input bit
         * @param[in] ct_2 encryption of the second input bit
         */
        void nand_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const;

        /**
         * Computes the AND gate of two given ciphertexts ct1 and ct2
         * @param[out] ct_res encryptions of the outuput of the AND gate
         * @param[in] ct_1 encryption of the first input bit
         * @param[in] ct_2 encryption of the second input bit
         */
        void and_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const;

        /**
         * Computes the OR gate of two given ciphertexts ct1 and ct2
         * @param[out] ct_res encryptions of the outuput of the OR gate
         * @param[in] ct_1 encryption of the first input bit
         * @param[in] ct_2 encryption of the second input bit
         */
        void or_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const;

        /**
         * Computes the XOR gate of two given ciphertexts ct1 and ct2
         * @param[out] ct_res encryptions of the outuput of the XOR gate
         * @param[in] ct_1 encryption of the first input bit
         * @param[in] ct_2 encryption of the second input bit
         */
        void xor_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const;

        /**
         * Computes the NOT gate of a given ciphertext ct
         * @param[out] ct_res encryption of the outuput of the NOT gate
         * @param[in] ct encryption of the input bit
         */
        void not_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct) const;
};

#endif
