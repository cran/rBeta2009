#include <algorithm>
#include <R.h>
#include <Rdefines.h>
#include "betaAlgs.h"

using namespace std;

/* Beta(shape1, shape2) */
#if __cplusplus
extern "C"
#endif
SEXP newrbeta(SEXP N, SEXP shape1, SEXP shape2) {
    PROTECT(N = AS_INTEGER(N));
    PROTECT(shape1 = AS_NUMERIC(shape1));
    PROTECT(shape2 = AS_NUMERIC(shape2));
    int *Npt = INTEGER(N);
    double *alpha1 = REAL(shape1), *alpha2 = REAL(shape2);
    if (*alpha1 <= 0.0 || *alpha2 <= 0.0) Rf_error("Shape parameters should be positive");
    SEXP out;
    PROTECT(out = NEW_NUMERIC(*Npt));
    double *outpt = REAL(out);
    R_len_t n;
    GetRNGstate();
    if (*alpha1 == 1.0 && *alpha2 == 1.0) { /* Uniform */
        for (n = 0; n < *Npt; ++n)
            *(outpt++) = unif_rand();
    } else if (*alpha1 == 1.0) { /* CDF inversion */
        for (n = 0; n < *Npt; ++n)
            *(outpt++) = (1.0 - power(unif_rand(), 1.0 / *alpha2));
    } else if (*alpha2 == 1.0) { /* CDF inversion */
        for (n = 0; n < *Npt; ++n)
            *(outpt++) = power(unif_rand(), 1.0 / *alpha1);
    } else if (*alpha1 < 1.0 && *alpha2 < 1.0) {
        B00(outpt, Npt, alpha1, alpha2);
    } else if (*alpha1 > 1.0 && *alpha2 > 1.0) {
        if ((fmin(*alpha1, *alpha2) <= 1.2 && fmax(*alpha1, *alpha2) > 4.0)) {
            B4PE(outpt, Npt, alpha1, alpha2);
        } else BPRS(outpt, Npt, alpha1, alpha2);
    } else B01(outpt, Npt, alpha1, alpha2);
    PutRNGstate();
    UNPROTECT(4);
    return out;
}

void shufrule(SEXP shape, int *invorder, int KK) {
    PROTECT(shape = AS_NUMERIC(shape));
    double *shapept = REAL(shape);
    int m = 0, i, j, ind1, ind2;
    /* Shuffling rule */
    double *orig, *s1, *s2;
    orig = new double[KK];
    for (i = 0; i < KK; i++)
        *(orig + i) = *(shapept + i);
    if (*min_element(shapept, shapept + KK) == *max_element(shapept, shapept + KK)) {
        for (i = 0; i < KK; i++)
            *(invorder + i) = i;
        UNPROTECT(1);
        return;
    } else {
        for (i = 0; i < KK; i++)
            if (*(shapept + i) < 1.0) ++m;
    }
    int half = KK / 2;
    sort(shapept, shapept + KK);
    if (*(shapept + KK - 1) < 1.0) {
        if (KK % 2 == 1) {
            s1 = new double[half + 1];
            s2 = new double[half];
            ind1 = 0;
            ind2 = 0;
            for (i = 0; i < KK; i++) {
                if (i % 2 == 0) {
                    *(s1 + (ind1++)) = *(shapept + i);
                } else *(s2 + (ind2++)) = *(shapept + i);
            }
            for (i = 0, ind1 = 0, ind2 = half - 1; i < KK; i++) {
                if (i < half + 1) {
                    *(shapept + i) = *(s1 + (ind1++));
                } else *(shapept + i) = *(s2 + (ind2--));
            }
        } else {
            s1 = new double[half];
            s2 = new double[half - 1];
            ind1 = 0;
            ind2 = 0;
            for (i = 1; i < KK; i++) {
                if (i % 2 == 1) {
                    *(s1 + (ind1++)) = *(shapept + i);
                } else *(s2 + (ind2++)) = *(shapept + i);
            }
            for (i = 1, ind1 = 0, ind2 = half - 2; i < KK; i++) {
                if (i < half + 1) {
                    *(shapept + i) = *(s1 + (ind1++));
                } else *(shapept + i) = *(s2 + (ind2--));
            }
        }
        delete [] s1;
        delete [] s2;
    } else if (*shapept > 1.0) {
        if (*(shapept + KK - 1) > 3.0) {
            swap(*shapept, *(shapept + KK - 3));
            sort(shapept + 1, shapept + KK - 2);
            reverse(shapept, shapept + KK);
        } else {
            swap(*shapept, *(shapept + KK - 2));
            sort(shapept + 1, shapept + KK - 1);
            reverse(shapept, shapept + KK);
        }
    } else if ((2 * m > KK + 1) && (*shapept <= 0.5)) {
        reverse(shapept + m - 1, shapept + KK - 1);
        reverse(shapept + m - 1, shapept + KK);
    } else {
        swap(*shapept, *(shapept + KK - 2));
        sort(shapept + 1, shapept + KK - 1);
        reverse(shapept, shapept + KK);
    }
    /* Inverse ordering */
    for (i = 0; i < KK; i++) {
        for (j = 0; j < KK; j++) {
            if (*(shapept + i) == *(orig + j)) {
                *(invorder + i) = j;
                *(orig + j) = -1.0;
                break;
            }
        }
    }
    delete [] orig;
    UNPROTECT(1);
    return;
}

/* The Dirichlet generating function, using the BETA-M algorithm proposed by Hung et al. (2011). */
#if __cplusplus
extern "C"
#endif
SEXP rdirichlet(SEXP N, SEXP shape) {
    PROTECT(N = AS_INTEGER(N));
    PROTECT(shape = AS_NUMERIC(shape));
    int NN = *INTEGER(N), KK = LENGTH(shape), *invorder, i, j;
    invorder = new int[KK];
    for (i = 0; i < KK; i++)
        if (*(REAL(shape) + i) <= 0.0) Rf_error("Shape parameters should be all positive");
    SEXP out, betas, vecsum, shape1, shapesum, alpha;
    PROTECT(out = Rf_allocMatrix(REALSXP, NN, KK));
    PROTECT(betas = NEW_NUMERIC(NN));
    PROTECT(vecsum = NEW_NUMERIC(NN));
    PROTECT(shape1 = NEW_NUMERIC(1));
    PROTECT(shapesum = NEW_NUMERIC(1));
    PROTECT(alpha = NEW_NUMERIC(KK));
    double *shape1pt = REAL(shape1), *shapesumpt = REAL(shapesum);
    *shapesumpt = 0.0;
    for (i = 0; i < KK; i++) {
        *(REAL(alpha) + i) = *(REAL(shape) + i);
        *shapesumpt += *(REAL(alpha) + i);
    }
    double *outpt = REAL(out), *vecsumpt = REAL(vecsum);
    shufrule(alpha, invorder, KK);
    GetRNGstate();
    *shape1pt = *REAL(alpha);
    *shapesumpt -= *shape1pt;
    betas = newrbeta(N, shape1, shapesum);
    for (j = 0; j < NN; j++) {
        *(outpt + *invorder * NN + j) = *(REAL(betas) + j);
        *(vecsumpt + j) = *(outpt + *invorder * NN + j);
    }
    for (i = 1; i < KK - 1; i++) {
        *shape1pt = *(REAL(alpha) + i);
        *shapesumpt -= *shape1pt;
        betas = newrbeta(N, shape1, shapesum);
        for (j = 0; j < NN; j++) {
            *(outpt + *(invorder + i) * NN + j) = (1.0 - *(vecsumpt + j)) * *(REAL(betas) + j);
            *(vecsumpt + j) += *(outpt + *(invorder + i) * NN + j);
        }
    }
    for (j = 0; j < NN; j++)
        *(outpt + *(invorder + KK - 1) * NN + j) = 1.0 - *(vecsumpt + j);
    PutRNGstate();
    delete [] invorder;
    UNPROTECT(8);
    return out;
}

