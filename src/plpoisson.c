#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

void R_init_plpoisson(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void plBinom(int *xo, double *p, double *a,
             int *up, int *lw) {
  /* Prediction limits with Binomial distribution
     INPUT:
       *xo - Number of observed occurrences (integer)
        *p - probability (real number between 0 and 1)
        *a - coverage probability (type I error)
       *up - upper bound (integer; used for output)
       *lw - lower bound (integer; used for output)
     OUTPUT:
       No output is returned.  The values of `up` and `lw` are replaced by the
       resulting upper and lower bounds respectively.
   */

  int s;
  double tmp;

  /* Compute the lower bound */
  *lw = 0;
  s = 1;
  while(pbinom(*xo, *xo + s, *p, 0, 0) < *a) {
      s <<= 1;
  }
  do {
    do {
      R_CheckUserInterrupt();
      s >>= 1;
      tmp = pbinom(*xo, *xo + *lw + s, *p, 0, 0);
/*      Rprintf("pbinom(%d, %d + %d + %d, p) = %g\n", *xo, *xo, *lw, s, tmp); */
    } while(tmp > *a);
    *lw += s;
  } while (s > 1);

  /* Compute the upper bound */
  *up = 0;
  s = 1;
  while(pbinom(*xo, *xo + s, *p, 1, 0) > *a) {
      s <<= 1;
  }
  do {
    do {
      R_CheckUserInterrupt();
      s >>= 1;
      tmp = pbinom(*xo, *xo + *up + s, *p, 1, 0);
/*      Rprintf("1 - pbinom(%d, %d + %d + %d, p) = %g\n", *xo, *xo, *up, s, tmp); */
    } while(tmp < *a);
    *up += s;
  } while (s > 1);
  *up += 1;
}

