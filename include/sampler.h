#ifndef SAMPLER
#define SAMPLER

#include <NTL/mat_ZZ.h>
#include <random>
#include <vector>
#include <chrono>
#include "params.h"

using namespace NTL;
using namespace std;


// random engine
static default_random_engine rand_engine(std::chrono::system_clock::now().time_since_epoch().count());
// uniform distribution on the ternary set
static uniform_int_distribution<int> ternary_sampler(-1,1);
// uniform distribution on the binary set
static uniform_int_distribution<int> binary_sampler(0,1);

class Sampler
{
    Param param;
    uniform_int_distribution<int> mod_q_base_sampler;

    public:
        Sampler(Param _param): param(_param)
        {
            mod_q_base_sampler = uniform_int_distribution<int>(-param.half_q_base, param.half_q_base);
        }

        /**
         * Generate a uniformly random vector modulo q_base.
         *
         * @param[out] vec vector with uniformly random coefficients.
         */
        void get_uniform_vector(vector<int>& vec);

        /**
         * Generate a uniformly random matrix modulo q_base.
         *
         * @param[out] mat matrix with uniformly random coefficients.
         */
        void get_uniform_matrix(vector<vector<int>>& mat);

        /**
         * Generate a matrix of the form scale*M+shift*I
         * where M has uniformly random ternary coefficients
         * and I is an identity matrix.
         * This matrix must be invertible 
         * @param[out] mat random matrix of the above form.
         * @param[out] mat_inv inverse matrix.
         * @param[in] scale scale in the above form.
         * @param[in] shift shift in the above form
         */
        void get_invertible_matrix(vector<vector<int>>& mat, vector<vector<int>>& mat_inv, int scale, int shift);

        /**
         * Generate a uniformly random matrix with ternary coefficients.
         *
         * @param[out] mat matrix with ternary coefficients.
         */
        static void get_ternary_matrix(vector<vector<int>>& mat);

        /**
         * Generate a uniformly random vector with ternary coefficients.
         *
         * @param[out] vec vector with ternary coefficients.
         */
        static void get_ternary_vector(vector<int>& vec);

        /**
         * Generate a uniformly random vector with binary coefficients.
         *
         * @param[out] vec vector with binary coefficients.
         */
        static void get_binary_vector(vector<int>& vec);

        /**
         * Generate a random matrix with coefficients distributed 
         * according to the discrete Gaussian distribution 
         * with zero mean and standard deviation st_dev.
         * @param[out] mat random matrix.
         * @param[in] st_dev standard deviation.
         */
        static void get_gaussian_matrix(vector<vector<int>>& mat, double st_dev);

        /**
         * Generate a random vector with coefficients distributed 
         * according to the discrete Gaussian distribution 
         * with zero mean and standard deviation st_dev.
         * @param[out] vec random vector.
         * @param[in] st_dev standard deviation.
         */
        static void get_gaussian_vector(vector<int>& vec, double st_dev);

        /**
         * Generate a vector of the form scale*v+shift
         * where v has uniformly random ternary coefficients.
         * The polynomial with the above vector of coefficients
         * must be invertible modulo X^N+1 and q_boot. 
         * @param[out] vec random vector of the above form.
         * @param[out] vec_inv coefficient vector of the polynomial inverse modulo X^N+1.
         * @param[in] scale scale in the above form.
         * @param[in] shift shift in the above form
         */
        void get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale, int shift);
};


#endif