// SVGPC native core. Mathematical contract: inst/doc/MATHEMATICS.md.
#include "native.h"
QR pivot_qr(const arma::mat &input, bool full) {
  int nr = input.n_rows, nc = input.n_cols, k = std::min(nr, nc), info, lwork = -1;
  double query;
  if (!nc)
    return {full ? arma::eye<mat>(nr, nr) : mat(nr, 0), 0};
  mat a = input;
  std::vector<int> piv(nc, 0);
  vec tau(k);
  svgpc_dgeqp3(&nr, &nc, a.memptr(), &nr, piv.data(), tau.memptr(), &query, &lwork, &info);
  lwork = std::max(1, (int)query);
  vec work(lwork);
  svgpc_dgeqp3(&nr, &nc, a.memptr(), &nr, piv.data(), tau.memptr(), work.memptr(), &lwork, &info);
  if (info)
    Rcpp::stop("Pivoted QR failed");
  double tol = std::numeric_limits<double>::epsilon() * std::max(nr, nc) * std::abs(a(0, 0));
  int rank = 0;
  for (int j = 0; j < k; ++j)
    if (std::abs(a(j, j)) > tol)
      ++rank;
  int nq = full ? nr : rank;
  if (!nq)
    return {mat(nr, 0), 0};
  mat q(nr, nq, arma::fill::zeros);
  q.cols(0, std::min(nc, nq) - 1) = a.cols(0, std::min(nc, nq) - 1);
  int reflect = std::min(k, nq);
  lwork = -1;
  svgpc_dorgqr(&nr, &nq, &reflect, q.memptr(), &nr, tau.memptr(), &query, &lwork, &info);
  lwork = std::max(1, (int)query);
  work.set_size(lwork);
  svgpc_dorgqr(&nr, &nq, &reflect, q.memptr(), &nr, tau.memptr(), work.memptr(), &lwork, &info);
  if (info)
    Rcpp::stop("QR basis construction failed");
  return {q, rank};
}
mat residual(mat a, const arma::mat &q) {
  if (q.n_cols)
    a -= q * (q.t() * a);
  return a;
}
mat kernel(const arma::mat &x, const arma::mat &k, double rho) {
  mat b(x.n_rows, k.n_rows);
  for (uword j = 0; j < k.n_rows; ++j)
    for (uword i = 0; i < x.n_rows; ++i) {
      double dx = x(i, 0) - k(j, 0), dy = x(i, 1) - k(j, 1);
      b(i, j) = std::exp(-std::hypot(dx, dy) / rho);
    }
  return b;
}
// [[Rcpp::export]]
Rcpp::List cpp_cluster(const arma::mat &xy, int m, int seed) {
  if (xy.n_cols != 2 || !xy.is_finite())
    Rcpp::stop("Clustering requires finite two-dimensional coordinates");
  if (m < 2 || m > (int)xy.n_rows)
    Rcpp::stop("Require 2 <= knots <= locations");
  std::vector<int> order(xy.n_rows);
  std::iota(order.begin(), order.end(), 0);
  std::mt19937 rng(seed);
  std::shuffle(order.begin(), order.end(), rng);
  mat kt(m, 2);
  for (int j = 0; j < m; ++j)
    kt.row(j) = xy.row(order[j]);
  std::vector<int> assignment(xy.n_rows, -1);
  vec counts(m);
  int iterations = 0;
  bool converged = false;
  for (int iteration = 0; iteration < 100; ++iteration) {
    iterations = iteration + 1;
    Rcpp::checkUserInterrupt();
    mat sums(m, 2, arma::fill::zeros);
    counts.zeros();
    int changed = 0;
    for (uword i = 0; i < xy.n_rows; ++i) {
      double best = INF;
      int win = 0;
      for (int j = 0; j < m; ++j) {
        double dx = xy(i, 0) - kt(j, 0), dy = xy(i, 1) - kt(j, 1), dist = dx * dx + dy * dy;
        if (dist < best) {
          best = dist;
          win = j;
        }
      }
      changed += (assignment[i] != win);
      assignment[i] = win;
      sums.row(win) += xy.row(i);
      counts[win] += 1;
    }
    // Empty centres disappear; all assigned locations remain in the analysis.
    arma::uvec keep = arma::find(counts > 0);
    if (keep.n_elem < (uword)m) {
      std::vector<int> remap(m);
      for (uword j = 0; j < keep.n_elem; ++j)
        remap[keep[j]] = j;
      for (int &a : assignment)
        a = remap[a];
      sums = sums.rows(keep);
      counts = counts.elem(keep);
      m = keep.n_elem;
      ++changed;
    }
    kt = sums.each_col() / counts;
    if (!changed) {
      converged = true;
      break;
    }
  }
  // Return membership with the same centres used by the spatial model.
  counts.zeros(m);
  for (uword i = 0; i < xy.n_rows; ++i) {
    double best = INF;
    for (int j = 0; j < m; ++j) {
      double dx = xy(i, 0) - kt(j, 0), dy = xy(i, 1) - kt(j, 1);
      double dist = dx * dx + dy * dy;
      if (dist < best) {
        best = dist;
        assignment[i] = j + 1;
      }
    }
    counts[assignment[i] - 1] += 1;
  }
  arma::uvec keep = arma::find(counts > 0);
  if (keep.n_elem < (uword)m) {
    std::vector<int> remap(m);
    for (uword j = 0; j < keep.n_elem; ++j)
      remap[keep[j]] = j + 1;
    for (int &a : assignment)
      a = remap[a - 1];
    kt = kt.rows(keep);
  }
  if (!converged)
    Rcpp::warning("Clustering reached 100 iterations; inspect the returned centres and membership");
  return Rcpp::List::create(Rcpp::_["centres"] = kt, Rcpp::_["cluster"] = assignment,
                            Rcpp::_["iterations"] = iterations, Rcpp::_["converged"] = converged);
}
vec nearest(const arma::mat &k) {
  vec d(k.n_rows);
  for (uword j = 0; j < k.n_rows; ++j) {
    double v = INF;
    for (uword i = 0; i < k.n_rows; ++i)
      if (i != j)
        v = std::min(v, std::hypot(k(i, 0) - k(j, 0), k(i, 1) - k(j, 1)));
    d[j] = v;
  }
  return d;
}
void eigen_top(const arma::mat &input, int count, vec &vals, mat &vectors) {
  mat a = arma::symmatu(input);
  int n = a.n_rows, lda = n, il = std::max(1, n - count + 1), iu = n, found = 0, info, lwork = -1,
      liwork = -1, iq;
  double vl = 0, vu = 0, tol = 0, query;
  char job = 'V', range = 'I', uplo = 'U';
  vec w(n);
  mat z(n, std::min(count, n));
  std::vector<int> support(2 * n);
  svgpc_dsyevr(&job, &range, &uplo, &n, a.memptr(), &lda, &vl, &vu, &il, &iu, &tol, &found,
               w.memptr(), z.memptr(), &lda, support.data(), &query, &lwork, &iq, &liwork, &info);
  lwork = (int)query;
  liwork = iq;
  vec work(lwork);
  std::vector<int> iwork(liwork);
  svgpc_dsyevr(&job, &range, &uplo, &n, a.memptr(), &lda, &vl, &vu, &il, &iu, &tol, &found,
               w.memptr(), z.memptr(), &lda, support.data(), work.memptr(), &lwork, iwork.data(),
               &liwork, &info);
  if (info)
    Rcpp::stop("Top eigendecomposition failed");
  double cutoff = std::numeric_limits<double>::epsilon() * n * std::max(w[found - 1], 0.0);
  std::vector<uword> keep;
  for (int j = found - 1; j >= 0; --j)
    if (w[j] > cutoff)
      keep.push_back(j);
  arma::uvec ind(keep);
  vals = w.elem(ind);
  vectors = z.cols(ind);
}
// [[Rcpp::export]]
SEXP cpp_prepare(const arma::mat &xy, const arma::mat &pcs, const arma::vec &counts,
                 const arma::mat &knots, double rho) {
  if (xy.n_rows < 3 || xy.n_cols != 2 || pcs.n_rows != xy.n_rows || counts.n_elem != xy.n_rows ||
      !xy.is_finite() || !pcs.is_finite() || !counts.is_finite() || arma::any(counts <= 0))
    Rcpp::stop("Invalid coordinates, covariates or location counts");
  std::unique_ptr<Context> c(new Context);
  c->xy = xy;
  c->weights = counts;
  c->rootw = arma::sqrt(counts);
  c->C = arma::join_rows(arma::ones<mat>(xy.n_rows, 1), pcs);
  auto t = Clock::now();
  c->knots = knots;
  if (c->knots.n_cols != 2 || c->knots.n_rows < 2 || !c->knots.is_finite())
    Rcpp::stop("Invalid supplied knots");
  vec nn = nearest(c->knots);
  if (nn.min() <= 0)
    Rcpp::stop("Duplicate knots");
  c->rho = rho;
  if (!(c->rho > 0) || !std::isfinite(c->rho))
    Rcpp::stop("Range must be positive and finite");
  c->timing["knots_and_spacing"] = seconds(t);
  t = Clock::now();
  auto step = Clock::now();
  mat Z = kernel(xy, c->knots, c->rho);
  c->timing["prepare_detail.01_location_kernel"] = seconds(step);
  step = Clock::now();
  mat Q = kernel(c->knots, c->knots, c->rho);
  c->timing["prepare_detail.02_prior_kernel"] = seconds(step);
  c->timing["kernel_basis_and_Q"] = seconds(t);
  t = Clock::now();
  step = Clock::now();
  mat X = c->C.each_col() % c->rootw;
  c->Xq = pivot_qr(X).Q;
  X.reset();
  c->timing["prepare_detail.03_fixed_effect_QR"] = seconds(step);
  if (c->r() < 1)
    Rcpp::stop("Fixed effects leave no residual degrees of freedom");
  step = Clock::now();
  Z.each_col() %= c->rootw;
  c->timing["prepare_detail.04_weight_basis"] = seconds(step);
  if (c->Xq.n_cols) {
    step = Clock::now();
    mat A = Z.t() * c->Xq;
    c->timing["prepare_detail.05_constraint_crossproduct"] = seconds(step);
    step = Clock::now();
    QR constraint = pivot_qr(A, true);
    A.reset();
    c->timing["prepare_detail.06_constraint_full_QR"] = seconds(step);
    if (constraint.rank >= (int)Q.n_rows)
      Rcpp::stop("Fixed-effect QR removes all spatial directions; increase clusters in svgpc_cluster()");
    step = Clock::now();
    c->N = constraint.Q.cols(constraint.rank, Q.n_rows - 1);
    c->timing["prepare_detail.07_null_basis_copy"] = seconds(step);
    step = Clock::now();
    Z = Z * c->N;
    c->timing["prepare_detail.08_basis_times_null_basis"] = seconds(step);
    step = Clock::now();
    Q = c->N.t() * Q * c->N;
    c->timing["prepare_detail.09_prior_constraint"] = seconds(step);
  } else {
    c->N = arma::eye<mat>(Q.n_rows, Q.n_rows);
  }
  // QR has already restricted the coefficient space: Xq' Z = 0.
  // Do not replace this reduced basis by residualizing the original columns.
  c->timing["weighted_QR_and_constraints"] = seconds(t);
  t = Clock::now();
  step = Clock::now();
  mat S = Z.t() * Z;
  c->timing["prepare_detail.10_basis_Gram"] = seconds(step);
  step = Clock::now();
  if (!arma::chol(Q, arma::symmatu(Q), "lower"))
    Rcpp::stop("Q is not numerically positive definite; no jitter added");
  c->timing["prepare_detail.11_prior_Cholesky"] = seconds(step);
  step = Clock::now();
  S = arma::solve(arma::trimatl(Q), S);
  c->timing["prepare_detail.12_whiten_left_solve"] = seconds(step);
  step = Clock::now();
  S = arma::solve(arma::trimatl(Q), S.t());
  c->timing["prepare_detail.13_whiten_right_solve"] = seconds(step);
  step = Clock::now();
  S = (S + S.t()) * 0.5;
  c->timing["prepare_detail.14_symmetrize"] = seconds(step);
  c->timing["Gram_and_Q_whitening"] = seconds(t);
  t = Clock::now();
  vec d;
  mat U;
  step = Clock::now();
  if (!arma::eig_sym(d, U, S))
    Rcpp::stop("Whitened eigendecomposition failed");
  c->timing["prepare_detail.15_symmetric_eigen"] = seconds(step);
  step = Clock::now();
  S.reset();
  double cutoff = std::numeric_limits<double>::epsilon() * std::max(c->n(), (int)Q.n_rows) *
                  std::max(d.max(), 1.0);
  if (d.min() < -cutoff)
    Rcpp::stop("Material negative whitened eigenvalue");
  arma::uvec keep = arma::find(d > cutoff);
  if (!keep.n_elem)
    Rcpp::stop("No estimable spatial directions");
  c->d = d.elem(keep);
  c->timing["prepare_detail.16_eigen_rank_filter"] = seconds(step);
  step = Clock::now();
  c->V = arma::solve(arma::trimatu(Q.t()), U.cols(keep));
  c->timing["prepare_detail.17_eigenvector_backsolve"] = seconds(step);
  step = Clock::now();
  U.reset();
  Q.reset();
  mat VD = c->V.each_row() / arma::sqrt(c->d).t();
  c->timing["prepare_detail.18_scale_eigenvectors"] = seconds(step);
  step = Clock::now();
  c->T = Z * VD;
  c->timing["prepare_detail.19_construct_design"] = seconds(step);
  step = Clock::now();
  Z.reset();
  VD.reset();
  c->timing["prepare_detail.20_release_workspaces"] = seconds(step);
  c->timing["eigendecomposition_and_design"] = seconds(t);
  t = Clock::now();
  step = Clock::now();
  // Phi = T / sqrt(W); form its metric in row blocks without storing Phi.
  int rank = c->rank();
  c->R.zeros(rank, rank);
  double one = 1;
  char lower = 'L', trans = 'T';
  for (int start = 0; start < c->n(); start += 256) {
    Rcpp::checkUserInterrupt();
    int rows = std::min(256, c->n() - start);
    mat block = c->T.rows(start, start + rows - 1);
    block.each_col() /= c->rootw.subvec(start, start + rows - 1);
    svgpc_dsyrk(&lower, &trans, &rank, &rows, &one, block.memptr(), &rows,
                 &one, c->R.memptr(), &rank);
  }
  c->R = arma::symmatl(c->R);
  c->timing["prepare_detail.21_raw_metric_Gram"] = seconds(step);
  step = Clock::now();
  if (!arma::chol(c->R, c->R))
    Rcpp::stop("Raw F metric is not positive definite");
  c->timing["prepare_detail.22_raw_metric_Cholesky"] = seconds(step);
  c->timing["raw_F_metric_once"] = seconds(t);
  step = Clock::now();
  c->cov.zeros(c->rank(), c->rank());
  c->timing["prepare_detail.23_allocate_covariance"] = seconds(step);
  return Rcpp::XPtr<Context>(c.release(), true);
}
