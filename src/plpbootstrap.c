#include <R.h>
#include <Rmath.h>

#define TOL_H 1e-9
#define NBOOT 1030
#define MAX_COUNT 20

/**
 * @brief Chi-square quantile ratio
 * @param a Pointer to a real scalar value
 * @return The ratio between the 95% and 5 % Chi-square quantiles
 */
static double chirat(double *a) {
    return Rf_qchisq(0.95, exp(*a) * 2.0, 1, 0) / Rf_qchisq(0.05, exp(*a) * 2.0, 1, 0);
}
/**
 * @brief Approximate first derivative of Chi-square quantile ratio
 * @param a Pointer to a real scalar value
 * @param tol Double precision floating-point value for accuracy tollerance
 * @return The approximation of the first derivative of a Chi-square quantile ratio
 */
static double derivata(double *a, double tol) {
    double d, p;
    d = *a - 0.5 * tol;
    p = *a + 0.5 * tol;
    d = chirat(&p) - chirat(&d);
    return d / tol;
}
/**
 * @brief Computes the Gamma hyper-parameters `a`, `b` from a bootstrap sample `x` of length `n` using the method of moments
 * @param x Pointer to an array of Poisson distributed data
 * @param n Pointer to the length of the array of data `x`
 * @param a Pointer to the Gamma hyper-parameters `a`
 * @param b Pointer to the Gamma hyper-parameters `b`
 */
void ABmoments(int *x, int *n, double *a, double *b) {
    int i;
    double s = 0.0, s2 = 0.0;
    double in = 1.0 / (double) *n;
    double m, r;

    /* Compute the momemnts of `x`... */
    for (i = 0; i < *n; i++) {
        s += (double) x[i];
        s2 += (double) x[i] * x[i];
    } /* ... to compute the moments of lambda */
    m = s * in;
    s = s2 * in - m * m;
    r = m / (s * in);

    /* Compute the estimates */
    *a = m * r;
    *b = 1.0 / r;
}

/**
 * @brief Computes the Gamma hyper-parameters `a`, `b` from a bootstrap sample `x` of length `n` using the maximum likelihood
 * @param x Pointer to an array of Poisson distributed data
 * @param n Pointer to the length of the array of data `x`
 * @param a Pointer to the Gamma hyper-parameters `a`
 * @param b Pointer to the Gamma hyper-parameters `b`
 */
void ABmle(int *x, int *n, double *a, double *b) {
    int i, bs;
    double lambda, in = 1.0 / (double) *n;
    double step, sl = 0.0, sll = 0.0;
    size_t count = 0;

    /* Bootstrapping with replacement */
    GetRNGstate();
    for (bs = 0; bs < NBOOT; bs++) {
        lambda = 0.0;
        for (i = 0; i < *n; i++)
            lambda += (double) x[(int) runif(0.0, (double) *n)];
        lambda *= in;
        sl += lambda;
        sll += log(lambda);
    }
    PutRNGstate();

    /* Compu/te the estimates */
    in = 1.0 / ((double) NBOOT);
    sl *= in;
    sll *= in;
    in = log(sl) - sll;
    *a = 0.5 / in;
    do {
        /* Find the zero of the equation via Gauss-Newton's step */
        step = digamma(*a) - log(*a) + in;
        step /= trigamma(*a) - 1.0 / *a;
        *a -= step;
        if (++count > MAX_COUNT) break;
    } while (fabs(step) > TOL_H);
    *b = sl / *a;
    if (isnan(*a) || *a <= 0.0) {
        *a = NAN;
        *b = NAN;
    }
}

/**
 * @brief Compute Chi-square approximation of `a` and `b` from the sample `x` of length `n`
 * @param x Pointer to an array of Poisson distributed data
 * @param n Pointer to the length of the array of data `x`
 * @param a Pointer to the Gamma hyper-parameters `a`
 * @param b Pointer to the Gamma hyper-parameters `b`
 */
void ABChisq(int *x, int *n, double *a, double *b) {
    int i, bs;
    double *lambda, in, ratio;
    double step;

    in = 1.0 / (double) *n;
    lambda = Calloc(NBOOT, double);

    /* Bootstrapping with replacement */
    GetRNGstate();
    for (bs = 0; bs < NBOOT; bs++) {
        lambda[bs] = 0.0;
        for (i = 0; i < *n; i++)
            lambda[bs] += (double) x[(int) runif(0.0, (double) *n)];
        lambda[bs] *= in;
    }
    PutRNGstate();

    /* Compute lambda quantiles */
    R_rsort(lambda, NBOOT);
    ratio = lambda[978] / lambda[51];

    /* Compute the estimates */
    *a = 0.0;
    do {
        /* Find the zero of the equation via Gauss-Newton's step */
        step = chirat(a) - ratio;
        step /= derivata(a, TOL_H);
        *a -= step;
    } while (fabs(step) > TOL_H);
    *a = exp(*a);
    *b = 2.0 * lambda[50] / Rf_qchisq(0.05, *a * 2.0, 1, 0);
    Free(lambda);
}

/**
 * @brief Bootstrap method to compute the hyperparameters based on three methods:
 *   1. Method of moments
 *   2. Maximum likelihood
 *   3. Chi-square approximation
 * @param x Pointer to an array of Poisson distributed data
 * @param n Pointer to the length of the array of data `x`
 * @param a Pointer to the Gamma hyper-parameters `a`
 * @param b Pointer to the Gamma hyper-parameters `b`
 * @param p Pointer to dimension of matrix `a` and `b` (i.e. number of bootstrap simulations)
 * @param met Pointer to the indices of the methods (with values from 1 to 3 as above)
 * @param m Pointer to the length of the array `met`
 */
void hyperbootstrap(int *x, int *n, double *a, double *b, int *p, int *met, int *m) {
    int i, j, k;
    int *xx;
    void (*fun[3])(int *, int *, double *, double *);

    fun[0] = ABmoments;
    fun[1] = ABmle;
    fun[2] = ABChisq;

    xx = Calloc(*n, int);
    GetRNGstate();
    if (xx) {
        for (k = 0; k < *m; k++) { /* For each chosen bootstrap method*/
            if (met[k] > 0 && met[k] < 4) {
                for (j = 0; j < *p; j++) { /* For each bootstrap iteration... */
                    /* Sub-sample with replacement */
                    for (i = 0; i < *n; i++) {
                        xx[i] = x[(int) runif(0.0, (double) *n)];
                    }
                    /* Compute hyper-prameter values*/
                    fun[met[k] - 1](xx, n, &a[*p * k + j], &b[*p * k + j]);
                }
            }
        }
        Free(xx);
    }
    PutRNGstate();
}
