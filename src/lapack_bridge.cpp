#include "lapack_bridge.h"
#include <R_ext/Lapack.h>
void svgpc_dgeqp3(int *m, int *n, double *a, int *lda, int *piv, double *tau, double *work,
                  int *lwork, int *info) {
  F77_CALL(dgeqp3)(m, n, a, lda, piv, tau, work, lwork, info);
}
void svgpc_dorgqr(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work,
                  int *lwork, int *info) {
  F77_CALL(dorgqr)(m, n, k, a, lda, tau, work, lwork, info);
}
void svgpc_dsyrk(char *u, char *t, int *n, int *k, double *alpha, double *a, int *lda, double *beta,
                 double *c, int *ldc) {
  F77_CALL(dsyrk)(u, t, n, k, alpha, a, lda, beta, c, ldc FCONE FCONE);
}
void svgpc_dsyevr(char *j, char *r, char *u, int *n, double *a, int *lda, double *vl, double *vu,
                  int *il, int *iu, double *tol, int *found, double *w, double *z, int *ldz,
                  int *support, double *work, int *lwork, int *iwork, int *liwork, int *info) {
  F77_CALL(dsyevr)(j, r, u, n, a, lda, vl, vu, il, iu, tol, found, w, z, ldz, support, work, lwork,
                   iwork, liwork, info FCONE FCONE FCONE);
}
