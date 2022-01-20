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
