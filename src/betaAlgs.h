#ifndef _BETAALGS_H
#define _BETAALGS_H

double power(double a, double b);

/* Beta generating algorithms */
void B00(double *out, int *N, double *shape1, double *shape2);
void B01(double *out, int *N, double *shape1, double *shape2);
void B4PE(double *out, int *N, double *shape1, double *shape2);
void BPRS(double *out, int *N, double *shape1, double *shape2);

#endif // _BETAALGS_H
