#include <R.h>
#include <Rmath.h>
#include "betaAlgs.h"

double power(double a, double b) {
    #if __x86_64__ || __ppc64__
    return(exp(b * log(a)));
    #else
    return(R_pow(a, b));
    #endif
}

void B00(double *out, int *N, double *shape1, double *shape2) {
    int n;
    double RN;
    static double U, V, t, s, r, p, q, c;
    /* Initialization */
    t = (1.0 - *shape1) / (2.0 - *shape1 - *shape2);
    s = (*shape2 - *shape1) * (1.0 - *shape1 - *shape2);
    r = *shape1 * (1.0 - *shape1);
    t = t - ((s * t + 2.0 * r) * t - r) / 2.0 / (s * t + r);
    p = t / *shape1;
    q = (1.0 - t) / *shape2;
    s = power(1.0 - t, *shape2 - 1.0);
    c = power(t, *shape1 - 1.0);
    r = (c - 1.0) / (t - 1.0);
    /* Generation */
    for (n = 0; n < *N; ++n) {
        for (;;) {
            U = unif_rand() * (p + q);
            V = unif_rand();
            if (U <= p) {
                RN = t * power(U / p, 1.0 / *shape1);
                if (s * V < (1.0 - *shape2) * RN + 1.0)
                    break;
                if ((s * V <= (s - 1.0) * RN / t + 1.0) && (s * V <= power(1.0 - RN, *shape2 - 1.0)))
                    break;
            } else {
                RN = 1.0 - (1.0 - t) * power((U - p) / q, 1.0 / *shape2);
                if (c * V < (*shape1 - 1.0) * (RN - 1.0) + 1.0)
                    break;
                if ((c * V <= r * (RN - 1.0) + 1.0) && (c * V <= power(RN, *shape1 - 1.0)))
                    break;
            }
        }
        *(out++) = RN;
    }
}

void B01(double *out, int *N, double *shape1, double *shape2) {
    int n;
    double RN;
    static double U, V, a, b, t, s, r, p ,q, c, d;
    /* Initialization */
    a = fmin(*shape1, *shape2);
    b = fmax(*shape1, *shape2);
    t = (1.0 - a) / (b - a);
    s = power(1.0 - t, b - 2);
    r = a - (a + b - 1.0) * t;
    t = t - (t - s * (1.0 - t) * (1 - r) / b) / (1 - s * r);
    p = t / a;
    q = power(1 - t, b - 1);
    s = fmin(1 - b, (q - 1) / t);
    r = fmax(1 - b, (q - 1) / t);
    q = q * (1 - t) / b;
    c = power(t, a - 1);
    d = (c - 1) / (t - 1);
    /* Generation */
    for (n = 0; n < *N; ++n) {
        for (;;) {
            U = unif_rand() * (p + q);
            V = unif_rand();
            if (U > p) {
                RN = 1.0 - (1.0 - t) * power((U - p) / q, 1.0 / b);
                if (c * V < (a - 1.0) * (RN - 1.0) + 1.0) break;
                if (c * V <= d * (RN - 1.0) + 1.0 && c * V <= power(RN, a - 1.0)) break;
            } else {
                RN = t * power(U / p, 1.0 / a);
                if (V < s * RN + 1.0) break;
                if (V <= r * RN + 1.0 && V <= power(1 - RN, b - 1.0)) break;
            }
        }
        *(out++) = (*shape1 < *shape2) ? RN : (1.0 - RN);
    }
}

bool B4PEaccept(double X, double x3, double V, double P, double Q, double R, double S) {
    double A = log(V);
    if (A > -(X - x3) * (X - x3) * (R + R)) return false;
    if (A <= P * log(X / P) + Q * log((1.0 - X) / Q) + S) return true;
    return false;
}

void B4PE(double *out, int *N, double *shape1, double *shape2) {
    int n;
    double RN;
    static double x1 = 0.0, x2 = 0.0, x3, x4 = 1.0, x5 = 1.0;
    static double f1 = 0.0, f2 = 0.0, f4 = 0.0, f5 = 0.0, pp[10];
    static double P, Q, D, R, S, U, V, W, W2;
    static double lambda1 = 0.0, lambda5 = 0.0;
    /* Initialization */
    P = fmin(*shape1, *shape2) - 1.0;
    Q = fmax(*shape1, *shape2) - 1.0;
    R = P + Q;
    S = R * log(R);
    x3 = P / R;
    if (R > 1.0) {
        D = sqrt(P * Q / (R - 1.0)) / R;
        if (D < x3) {
            x2 = x3 - D;
            x1 = x2 - x2 * (1.0 - x2) / (P - R * x2);
            lambda1 = P / x1 - Q / (1.0 - x1);
            f1 = exp(P * log(x1 / P) + Q * log((1.0 - x1) / Q) + S);
            f2 = exp(P * log(x2 / P) + Q * log((1.0 - x2) / Q) + S);
        }
        if (x3 + D < 1.0) {
            x4 = x3 + D;
            x5 = x4 - x4 * (1.0 - x4) / (P - R * x4);
            lambda5 = Q / (1.0 - x5) - P / x5;
            f4 = exp(P * log(x4 / P) + Q * log((1.0 - x4) / Q) + S);
            f5 = exp(P * log(x5 / P) + Q * log((1.0 - x5) / Q) + S);
        }
    }
    pp[0] = f2 * (x3 - x2);
    pp[1] = f4 * (x4 - x3) + pp[0];
    pp[2] = f1 * (x2 - x1) + pp[1];
    pp[3] = f5 * (x5 - x4) + pp[2];
    pp[4] = (1.0 - f2) * (x3 - x2) + pp[3];
    pp[5] = (1.0 - f4) * (x4 - x3) + pp[4];
    pp[6] = (f2 - f1) * (x2 - x1) / 2 + pp[5];
    pp[7] = (f4 - f5) * (x5 - x4) / 2 + pp[6];
    pp[8] = (lambda1 > 0.0) ? f1 / lambda1 + pp[7] : pp[7];
    pp[9] = (lambda5 > 0.0) ? f5 / lambda5 + pp[8] : pp[8];
    /* Generation */
    for (n = 0; n < *N; ++n) {
        for (;;) {
            U = unif_rand() * pp[9];
            if (U > pp[3]) {
                W = unif_rand();
                if (U > pp[4]) {
                    if (U > pp[5]) {
                        if (U <= pp[7]) {
                            W2 = unif_rand();
                            if (W2 > W) W = W2;
                            if (U > pp[6]) {
                                RN = x5 - W * (x5 - x4);
                                V = f5 + 2 * W * (U - pp[6]) / (x5 - x4);
                                if (V <= f4 * W) break;
                                if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                                continue;
                            } else {
                                RN = x1 + W * (x2 - x1);
                                V = f1 + 2 * W * (U - pp[5]) / (x2 - x1);
                                if (V < f2 * W) break;
                                if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                                continue;
                            }
                        }
                        if (U > pp[8]) {
                            U = (pp[9] - U) / (pp[9] - pp[8]);
                            RN = x5 - log(U) / lambda5;
                            if (RN >= 1.0) continue;
                            if (W <= (lambda5 * (x5 - RN) + 1.0) / U) break;
                            V = W * f5 * U;
                            if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                        } else {
                            U = (pp[8] - U) / (pp[8] - pp[7]);
                            RN = x1 + log(U) / lambda1;
                            if (RN <= 0.0) continue;
                            if (W <= (lambda1 * (RN - x1) + 1.0) / U) break;
                            V = W * f1 * U;
                            if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                            continue;
                        }
                    } else {
                        RN = x3 + W * (x4 - x3);
                        if ((pp[5] - U) / (pp[5] - pp[4]) >= W) break;
                        V = f4 + (U - pp[4]) / (x4 - x3);
                        if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                        continue;
                    }
                } else {
                    RN = x2 + W * (x3 - x2);
                    if ((U - pp[3]) / (pp[4] - pp[3]) <= W) break;
                    V = f2 + (U - pp[3]) / (x3 - x2);
                    if (B4PEaccept(RN, x3, V, P, Q, R, S)) break;
                    continue;
                }
            } else if (U > pp[0]) {
                if (U > pp[1]) {
                    if (U > pp[2]) {
                        RN = x4 + (U - pp[2]) / f5;
                        break;
                    } else {
                        RN = x1 + (U - pp[1]) / f1;
                        break;
                    }
                } else {
                    RN = x3 + (U - pp[0]) / f4;
                    break;
                }
            } else {
                RN = x2 + U / f2;
                break;
            }
        }
        *(out++) = (*shape1 < *shape2) ? RN : (1.0 - RN);
    }
}

double BPRSrescalef(double X, double a, double b, double m) {
    return (power(X / m, a) * power((1.0 - X) / (1.0 - m), b));
}

void BPRS(double *out, int *N, double *shape1, double *shape2) {
    int n;
    double RN;
    static double x1, x2, x4, x5, m, f1, f2, f4, f5, z2, z4, pp[4];
    static double D, U, V, W, Y, a, b, s;
    static double lambda_l, lambda_r, D_l, D_r;
    /* Initialization */
    a = fmin(*shape1, *shape2) - 1.0;
    b = fmax(*shape1, *shape2) - 1.0;
    s = a + b;
    m = a / s;
    if (a > 1.0 || b > 1.0)
        D = sqrt(m * (1.0 - m) / (s - 1.0));
    if (a <= 1.0) {
        D_l = m / 2;
        x2 = D_l;
        x1 = 0.0;
        z2 = 0.0;
        f1 = 0.0;
        lambda_l = 0.0;
    } else {
        D_l = D;
        x2 = m - D_l;
        x1 = x2 - D_l;
        z2 = x2 * (1.0 - (1.0 - x2) / (s * (m - x2)));
        if ((x1 <= 0.0) || (x2 * (s - 6.0) - a + 3.0 > 0.0)) {
            x1 = z2;
            x2 = (x1 + m) / 2.0;
            D_l = m - x2;
        }
        f1 = BPRSrescalef(x1, a, b, m);
        lambda_l = x1 * (1.0 - x1) / (s * (m - x1));
    }
    f2 = BPRSrescalef(x2, a, b, m);
    if (b <= 1.0) {
        D_r = (1.0 - m) / 2.0;
        x4 = 1.0 - D_r;
        x5 = 1.0;
        z4 = 1.0;
        f5 = 0.0;
        lambda_r = 0.0;
    } else {
        D_r = D;
        x4 = m + D_r;
        x5 = x4 + D_r;
        z4 = x4 * (1.0 + (1.0 - x4) / (s * (x4 - m)));
        if ((x5 >= 1.0) || (x4 * (s - 6.0) - a + 3.0 < 0.0)) {
            x5 = z4;
            x4 = (x5 + m) / 2.0;
            D_r = x4 - m;
        }
        f5 = BPRSrescalef(x5, a, b, m);
        lambda_r = x5 * (1.0 - x5) / (s * (x5 - m));
    }
    f4 = BPRSrescalef(x4, a, b, m);
    pp[0] = f2 * (m - x1);
    pp[1] = f4 * (x5 - m) + pp[0];
    pp[2] = f1 * lambda_l + pp[1];
    pp[3] = f5 * lambda_r + pp[2];
    /* Generation */
    for (n = 0; n < *N; ++n) {
        for (;;) {
            U = unif_rand() * pp[3];
            if (U <= pp[0]) { // Step 2
                Y = U / D_l - f2;
                if (Y <= 0.0) {
                    RN = m - U / f2;
                    break;
                }
                if (Y <= f1) {
                    RN = x2 - Y / f1 * D_l;
                    break;
                }
                U = unif_rand();
                W = U * D_l;
                RN = x2 - W;
                V = x2 + W;
                if (Y <= f2 * (RN - z2) / (x2 - z2)) break;
                W = f2 + f2 - Y;
                if (W < 1.0) {
                    if (W <= f2 + (1.0 - f2) * U) {
                        RN = V;
                        break;
                    }
                    if (W <= BPRSrescalef(V, a, b, m)) {
                        RN = V;
                        break;
                    }
                }
            } else if (U <= pp[1]) { // Step 4
                Y = (U - pp[0]) / D_r - f4;
                if (Y <= 0.0) {
                    RN = m + (U - pp[0]) / f4;
                    break;
                }
                if (Y <= f5) {
                    RN = x4 + Y / f5 * D_r;
                    break;
                }
                U = unif_rand();
                W = U * D_r;
                RN = x4 + W;
                V = x4 - W;
                if (Y <= f4 * (z4 - RN) / (z4 - x4)) break;
                W = f4 + f4 - Y;
                if (W < 1.0) {
                    if (W <= f4 + (1.0 - f4) * U) {
                        RN = V;
                        break;
                    }
                    if (W <= BPRSrescalef(V, a, b, m)) {
                        RN = V;
                        break;
                    }
                }
            } else if (U <= pp[2]) { // Step 6
                U = (U - pp[1]) / (pp[2] - pp[1]);
                V = log(U);
                RN = x1 + lambda_l * V;
                if (RN <= 0.0) continue;
                Y = unif_rand() * U;
                if (Y <= 1.0 + V) break;
                Y = Y * f1;
            } else { // Step 7
                U = (U - pp[2]) / (pp[3] - pp[2]);
                V = log(U);
                RN = x5 - lambda_r * V;
                if (RN >= 1.0) continue;
                Y = unif_rand() * U;
                if (Y <= 1.0 + V) break;
                Y = Y * f5;
            }
            if (Y <= BPRSrescalef(RN, a, b, m)) break;
        }
        *(out++) = (*shape1 < *shape2) ? RN : (1.0 - RN);
    }
}
