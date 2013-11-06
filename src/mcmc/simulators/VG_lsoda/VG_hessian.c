#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_errno.h>
/*This code was generously provided  by John Lamb from GSL discussion forum*/

/** Estimate the Hessian of a function f. The vector x and matrix H should satisfy
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
	     void* params, double h, gsl_matrix* H ){
  /* variables */
  size_t n;              /* size of x (and H) */
  size_t row, col;       /* for indexing H    */
  long double pp, pm, mp, mm; /* temporaries for calculations */
  long double cc, twh2, fh2, th;
  long double xrow, xcol;     /* local storage */
  gsl_vector_view view;  /* last row of H */
  
  /* carry out some sanity checks */
  if( x == 0 || f == 0 || H == 0 )
    return GSL_EINVAL;
  n = x->size;
  if( H->size1 != n || H->size2 != n )
    return GSL_EINVAL;
  if( h <= 0 )
    return GSL_EDOM;
  
  /* set up */
  view = gsl_matrix_row( H, n - 1 );
  gsl_vector_memcpy( &view.vector, x );  /* use a copy of x */
  th = 2 * h;
  fh2 = th * th;
  twh2 = 3 * fh2;
  cc = f( &view.vector, params );
  
  /* calculate upper triangle */
  for( row = 0; row < n; ++row ){
    for( col = row; col < n; ++col ){
      /* copy element to be modified */
      xrow = gsl_vector_get( &view.vector, row );
      if( row == col ){
	gsl_vector_set( &view.vector, row, xrow + th );
	pp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow + h );
	pm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - th );
	mm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - h );
	mp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow );
	/* set matrix entry */
	gsl_matrix_set( H, row, col, (16 * (pm + mp) - (pp + 30 * cc + mm)) / twh2 );
	/* very last call will overwrite view and view should not then be changed */
      } else {
	/* copy element to be modified */
	xcol = gsl_vector_get( &view.vector, col );
	gsl_vector_set( &view.vector, row, xrow + h );
	gsl_vector_set( &view.vector, col, xcol + h );
	pp = f( &view.vector, params );
	gsl_vector_set( &view.vector, col, xcol - h );
	pm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - h );
	mm = f( &view.vector, params );
	gsl_vector_set( &view.vector, col, xcol + h );
	mp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow );
	gsl_vector_set( &view.vector, col, xcol );
	/* set matrix entry */
	gsl_matrix_set( H, row, col, (pp + mm - (pm + mp)) / fh2 );
      }
    }
  }

  /* calculate lower triangle */
  for( row = 1; row < n; ++row ){ /* will work even for 1 x 1 matrix */
    for( col = 0; col < row; ++col ){
      gsl_matrix_set( H, row, col, gsl_matrix_get( H, col, row ) );
    }
  }

  return 0;
}
