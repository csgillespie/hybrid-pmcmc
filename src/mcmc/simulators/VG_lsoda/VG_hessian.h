#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_deriv.h>

#ifndef HESSIAN_H
#define HESSIAN_H

/**
 * Estimate the Hessian of a function f. The vector x and matrix H should satisfy
 * x->size == H->size1 == H->size2 and h should be positive. The function checks
 * these values and returns GSL_EINVAL or GSL_EDOM if they are not satisfied.
 *
 * The Hessian is estimated with central difference quotient estimates and uses
 * \f$4n^2-2n+4\f$ function calls.
 * 
 * @param x A vector.
 * @param f A function that takes as its first argument a vector of the same size
 * as x. The void* is used to pass any additional function parameters.
 * @param params Additional parameters to be given to each function call.
 * @param h A small positive constant used in the calculations. See standard texts
 * on numerical analysis for typical values.
 * @param H Stores the resulting Hessian.
 */
int hessian( const gsl_vector* x, double (*f)( const gsl_vector*, void* ),
	     void* params, double h, gsl_matrix* H );

#endif
