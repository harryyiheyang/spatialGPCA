#ifndef SVGPC_NATIVE_H
#define SVGPC_NATIVE_H
#include "lapack_bridge.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
using arma::mat;
using arma::uword;
using arma::vec;
using Clock = std::chrono::steady_clock;
inline double seconds(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}
const double INF = std::numeric_limits<double>::infinity();

struct QR {
  mat Q;
  int rank;
};
struct Context {
  mat xy, knots, C, Xq, N, V, T, R, cov;
  vec weights, rootw, d, a;
  double rho = 0, energy = 0, b = 0;
  uint64_t p = 0, missing = 0, read_variants = 0;
  bool ready = false;
  std::map<std::string, double> timing;
  int n() const { return weights.n_elem; }
  int rank() const { return d.n_elem; }
  int r() const { return n() - Xq.n_cols; }
};

QR pivot_qr(const arma::mat &input, bool full = false);
mat residual(mat a, const arma::mat &q);
void eigen_top(const arma::mat &input, int count, vec &vals, mat &vectors);
Rcpp::XPtr<Context> context(SEXP x);
Rcpp::List information(Context &c);
void finish(Context &c);
#endif
