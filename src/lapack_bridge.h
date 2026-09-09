#ifndef SVGPC_LAPACK_BRIDGE_H
#define SVGPC_LAPACK_BRIDGE_H
// Isolate R's BLAS/LAPACK declarations from Armadillo's complex declarations.
void svgpc_dgeqp3(int *, int *, double *, int *, int *, double *, double *, int *, int *);
void svgpc_dorgqr(int *, int *, int *, double *, int *, double *, double *, int *, int *);
void svgpc_dsyrk(char *, char *, int *, int *, double *, double *, int *, double *, double *,
                 int *);
void svgpc_dsyevr(char *, char *, char *, int *, double *, int *, double *, double *, int *, int *,
                  double *, int *, double *, double *, int *, int *, double *, int *, int *, int *,
                  int *);
#endif
