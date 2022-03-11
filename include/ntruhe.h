#ifndef NTRUHE
#define NTRUHE

#include "params.h"
#include "keygen.h"

class Ctxt_NTRU
{
    public:
    std::vector<int> data;

    Ctxt_NTRU()
    {
        data.clear();
        data.resize(parNTRU.n);
    }
    Ctxt_NTRU(const Ctxt_NTRU& ct);
    Ctxt_NTRU& operator=(const Ctxt_NTRU& ct);

    Ctxt_NTRU operator +(const Ctxt_NTRU& ct) const;
    Ctxt_NTRU operator -(const Ctxt_NTRU& ct) const;
    void operator -=(const Ctxt_NTRU& ct);
};

/**
 * Switches a given ciphertext to a given modulus.
 * @param[in,out] ct ciphertext
 * @param[in] old_q old modulus
 * @param[in] new_q old modulus
 */ 
inline void modulo_switch_ntru(Ctxt_NTRU& ct, int old_q, int new_q)
{
    std::vector<int>& a = ct.data;
    for (size_t i = 0; i < a.size(); i++)
        a[i] = int((a[i]*new_q)/old_q);
}

/**
 * Switches a given polynomial from q_base to modulus 2*N.
 * @param[in,out] poly polynomial
 */ 
inline void modulo_switch_to_boot(Ctxt_NTRU& poly)
{
    modulo_switch_ntru(poly, parNTRU.q_base, Param::N2);
}

/**
 * Switches a given polynomial from q_boot to q_base.
 * @param[in,out] poly polynomial
 */ 
inline void modulo_switch_to_base_ntru(ModQPoly& poly)
{
    modulo_switch(poly, q_boot, parNTRU.q_base);
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
//void external_product(std::vector<long>& res, const std::vector<int>& poly, const std::vector<FFTPoly>& poly_vector, const int b, const int shift, const int l);

class SchemeNTRU
{
    SKey_base_NTRU sk_base;
    SKey_boot sk_boot;
    KSKey_NTRU ksk;
    BSKey_NTRU bk;

    Ctxt_NTRU ct_nand_const;
    Ctxt_NTRU ct_and_const;
    Ctxt_NTRU ct_or_const;
    Ctxt_NTRU ct_not_const;

    void mask_constant(Ctxt_NTRU& ct, int constant);

    inline void set_nand_const()
    {
        //clock_t start = clock();
        mask_constant(ct_nand_const, parNTRU.nand_const);
        //cout << "Encryption of NAND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    }

    inline void set_and_const()
    {
        //clock_t start = clock();
        mask_constant(ct_and_const, parNTRU.and_const);
        //cout << "Encryption of AND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    }

    inline void set_or_const()
    {
        //clock_t start = clock();
        mask_constant(ct_or_const, parNTRU.or_const);
        //cout << "Encryption of OR: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    }

    inline void set_not_const()
    {
        encrypt(ct_not_const, 1);
    }



    public:
    
    SchemeNTRU()
    {
        KeyGen keygen(parNTRU);

        keygen.get_sk_base(sk_base);
        keygen.get_sk_boot(sk_boot);
        keygen.get_ksk(ksk,sk_base,sk_boot);
        keygen.get_bsk(bk,sk_base,sk_boot);

        set_nand_const();
        set_and_const();
        set_or_const();
        set_not_const();
    }
    /**
     * Encrypts a bit using matrix NTRU.
     * @param[out] ct ciphertext encrypting the input bit 
     * @param[in] b bit to encrypt
     */ 
    void encrypt(Ctxt_NTRU& ct, const int b) const;

    /**
     * Decrypts a bit using matrix NTRU.
     * @param[out] ct ciphertext encrypting a bit
     * @return b bit
     */ 
    int decrypt(const Ctxt_NTRU& ct) const;
    
    /**
     * Performs key switching of a given ciphertext from a polynomial NTRU
     * to a matrix NTRU
     * @param[out] ct matrix NTRU ciphertext (vector of dimension n)
     * @param[in] poly polynomial ciphertext (vector of dimension N)
     */
    void key_switch(Ctxt_NTRU& ct, const ModQPoly& poly) const;

    /**
     * Bootstrapps a given ciphertext
     * @param[in,out] ct ciphertext to bootstrap
     */
    void bootstrap(Ctxt_NTRU& ct) const;

    /**
     * Computes the NAND gate of two given ciphertexts ct1 and ct2
     * @param[out] ct_res encryptions of the outuput of the NAND gate
     * @param[in] ct_1 encryption of the first input bit
     * @param[in] ct_2 encryption of the second input bit
     */
    void nand_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const;

    /**
     * Computes the AND gate of two given ciphertexts ct1 and ct2
     * @param[out] ct_res encryptions of the outuput of the AND gate
     * @param[in] ct_1 encryption of the first input bit
     * @param[in] ct_2 encryption of the second input bit
     */
    void and_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const;

    /**
     * Computes the OR gate of two given ciphertexts ct1 and ct2
     * @param[out] ct_res encryptions of the outuput of the OR gate
     * @param[in] ct_1 encryption of the first input bit
     * @param[in] ct_2 encryption of the second input bit
     */
    void or_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const;

    /**
     * Computes the XOR gate of two given ciphertexts ct1 and ct2
     * @param[out] ct_res encryptions of the outuput of the XOR gate
     * @param[in] ct_1 encryption of the first input bit
     * @param[in] ct_2 encryption of the second input bit
     */
    void xor_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct1, const Ctxt_NTRU& ct2) const;

    /**
    * Computes the NOT gate of a given ciphertext ct
    * @param[out] ct_res encryption of the outuput of the NOT gate
    * @param[in] ct encryption of the input bit
    */
    void not_gate(Ctxt_NTRU& ct_res, const Ctxt_NTRU& ct) const;

};

#endif
