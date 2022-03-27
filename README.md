# FINAL 

## Faster FHE instantiated with NTRU and LWE

The FINAL library contains the implementation of the fully homorphic encryption schemes
presented in the paper ["FINAL: Faster FHE instantiated with NTRU and LWE"](http://eprint.iacr.org/2022/074),
by Charlotte Bonte (<charlotte.bonte@intel.com>), Ilia Iliashenko (<ilia@esat.kuleuven.be>), Jeongeun Park (<Jeongeun.Park@esat.kuleuven.be>), Hilder V. L. Pereira (<HilderVitor.LimaPereira@esat.kuleuven.be>), and Nigel P. Smart (<nigel.smart@kuleuven.be>).

It is distributed under the MIT license. Please, check the LICENSE file for more details.

### Requirements 

A C++ compiler, the [NTL](https://libntl.org) and [FFTW 3](http://www.fftw.org) libraries.

## Run the code

1. Run `make` in the main repository folder.
2. Run the `test` program and check that all the homomorphic gates are computed correctly.

## Usage

Use `test.cpp` and `Makefile` as reference points to create and compile your own program with FINAL.

### Example
```c++
// Input bits
int b1 = 0;
int b2 = 1;

// LWE encryption base scheme
SchemeLWE s;
// LWE ciphertexts
Ctxt_LWE ct1, ct2, ct_or, ct_nand, ct_xor, ct_and, ct_not;
// Encryption of bits
s.encrypt(ct1, b1);
s.encrypt(ct2, b2);

// Computes AND
s.and_gate(ct_and, ct1, ct2);
assert(s.decrypt(ct_and) == 0);
// Computes NAND
s.nand_gate(ct_nand, ct1, ct2);
assert(s.decrypt(ct_nand) == 1);
// Computes OR
s.or_gate(ct_or, ct1, ct2);
assert(s.decrypt(ct_or) == 1);
// Computes XOR
s.xor_gate(ct_xor, ct1, ct2);
assert(s.decrypt(ct_xor) == 1);
// Computes NOT
s.not_gate(ct_not, ct1);
assert(s.decrypt(ct_not) == 1);
```
